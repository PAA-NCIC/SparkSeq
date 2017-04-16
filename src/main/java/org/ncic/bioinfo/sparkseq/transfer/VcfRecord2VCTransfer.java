/*
 * Copyright (c) 2017 NCIC, Institute of Computing Technology, Chinese Academy of Sciences
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package org.ncic.bioinfo.sparkseq.transfer;

import htsjdk.tribble.util.ParsingUtils;
import htsjdk.variant.utils.GeneralUtils;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.LazyGenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.AbstractVCFCodec;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.VCFHeaderLineIterable;
import org.ncic.bioinfo.sparkseq.data.basic.VcfRecord;
import org.ncic.bioinfo.sparkseq.data.common.VcfHeaderInfo;
import org.ncic.bioinfo.sparkseq.exceptions.PipelineException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * Author: wbc
 */
public class VcfRecord2VCTransfer {
    private VCFHeader vcfFileHeader;
    private final static int NUM_STANDARD_FIELDS = 8;
    private String name = "Unknown";
    protected String[] infoFieldArray = new String[1000];
    protected String[] infoValueArray = new String[1000];

    // for performance we cache the hashmap of filter encodings for quick lookup
    protected HashMap<String, List<String>> filterHash = new HashMap<String, List<String>>();

    public VcfRecord2VCTransfer(VcfHeaderInfo vcfHeaderInfo) {
        VCFCodec codec = new VCFCodec();
        VCFHeaderLineIterable headerLineIterable = new VCFHeaderLineIterable(vcfHeaderInfo);
        vcfFileHeader = (VCFHeader) codec.readActualHeader(headerLineIterable);
    }

    public VariantContext transfer(VcfRecord vcfRecord) {
        VariantContextBuilder builder = new VariantContextBuilder();
        throw new PipelineException("Not implemented");
        /*builder.source(name);
        String chr = vcfRecord.contigName();
        builder.chr(chr);
        int pos = vcfRecord.position();
        builder.start(pos);
        if (vcfRecord.id().equals(".")) {
            builder.noID();
        } else {
            builder.id(vcfRecord.id());
        }
        String ref = vcfRecord.ref();
        String alts = vcfRecord.alt();
        builder.log10PError(parseQual(vcfRecord.qual()));

        final List<String> filters = parseFilters(vcfRecord.filter());
        if (filters != null) {
            builder.filters(new HashSet<String>(filters));
        }
        final Map<String, Object> attrs = parseInfo(vcfRecord.attributes());
        builder.attributes(attrs);

        if (attrs.containsKey(VCFConstants.END_KEY)) {
            // update stop with the end key if provided
            builder.stop(Integer.valueOf(attrs.get(VCFConstants.END_KEY).toString()));
        } else {
            builder.stop(pos + ref.length() - 1);
        }

        // get our alleles, filters, and setup an attribute map
        final List<Allele> alleles = parseAlleles(ref, alts, lineNo);
        builder.alleles(alleles);

        // do we have genotyping data
        if (parts.length > NUM_STANDARD_FIELDS) {
            final LazyGenotypesContext.LazyParser lazyParser = new AbstractVCFCodec.LazyVCFGenotypesParser(alleles, chr, pos);
            final int nGenotypes = vcfFileHeader.getNGenotypeSamples();
            LazyGenotypesContext lazy = new LazyGenotypesContext(lazyParser, parts[8], nGenotypes);

            // did we resort the sample names?  If so, we need to load the genotype data
            if (!vcfFileHeader.samplesWereAlreadySorted())
                lazy.decode();

            builder.genotypesNoValidation(lazy);
        }*/
    }

    /**
     * parse out the qual value
     *
     * @return return a double
     */
    protected static Double parseQual(double qual) {
        // if we're the VCF 4 missing char, return immediately
        if (qual == 0)
            return VariantContext.NO_LOG10_PERROR;

        Double val = Double.valueOf(qual);

        // check to see if they encoded the missing qual score in VCF 3 style, with either the -1 or -1.0.  check for val < 0 to save some CPU cycles
        if ((val < 0) && (Math.abs(val - VCFConstants.MISSING_QUALITY_v3_DOUBLE) < VCFConstants.VCF_ENCODING_EPSILON))
            return VariantContext.NO_LOG10_PERROR;

        // scale and return the value
        return val / -10.0;
    }

    protected List<String> parseFilters(final String filterString) {
        // null for unfiltered
        if (filterString.equals(VCFConstants.UNFILTERED))
            return null;

        if (filterString.equals(VCFConstants.PASSES_FILTERS_v4))
            return Collections.emptyList();

        // do we have the filter string cached?
        if (filterHash.containsKey(filterString))
            return filterHash.get(filterString);

        // empty set for passes filters
        final List<String> fFields = new LinkedList<String>();
        // otherwise we have to parse and cache the value
        if (!filterString.contains(VCFConstants.FILTER_CODE_SEPARATOR))
            fFields.add(filterString);
        else
            fFields.addAll(Arrays.asList(filterString.split(VCFConstants.FILTER_CODE_SEPARATOR)));

        filterHash.put(filterString, Collections.unmodifiableList(fFields));

        return fFields;
    }

    /**
     * parse out the info fields
     *
     * @param infoField the fields
     * @return a mapping of keys to objects
     */
    private Map<String, Object> parseInfo(String infoField) {
        Map<String, Object> attributes = new HashMap<String, Object>();

        if (!infoField.equals(VCFConstants.EMPTY_INFO_FIELD)) {
            int infoFieldSplitSize = ParsingUtils.split(infoField, infoFieldArray, VCFConstants.INFO_FIELD_SEPARATOR_CHAR, false);
            for (int i = 0; i < infoFieldSplitSize; i++) {
                String key;
                Object value;

                int eqI = infoFieldArray[i].indexOf("=");
                if (eqI != -1) {
                    key = infoFieldArray[i].substring(0, eqI);
                    String valueString = infoFieldArray[i].substring(eqI + 1);

                    // split on the INFO field separator
                    int infoValueSplitSize = ParsingUtils.split(valueString, infoValueArray, VCFConstants.INFO_FIELD_ARRAY_SEPARATOR_CHAR, false);
                    if (infoValueSplitSize == 1) {
                        value = infoValueArray[0];
                        final VCFInfoHeaderLine headerLine = vcfFileHeader.getInfoHeaderLine(key);
                        if (headerLine != null && headerLine.getType() == VCFHeaderLineType.Flag && value.equals("0")) {
                            // deal with the case where a flag field has =0, such as DB=0, by skipping the add
                            continue;
                        }
                    } else {
                        ArrayList<String> valueList = new ArrayList<String>(infoValueSplitSize);
                        for (int j = 0; j < infoValueSplitSize; j++)
                            valueList.add(infoValueArray[j]);
                        value = valueList;
                    }
                } else {
                    key = infoFieldArray[i];
                    final VCFInfoHeaderLine headerLine = vcfFileHeader.getInfoHeaderLine(key);
                    if (headerLine != null && headerLine.getType() != VCFHeaderLineType.Flag) {
                        value = VCFConstants.MISSING_VALUE_v4;
                    } else {
                        value = true;
                    }
                }

                // this line ensures that key/value pairs that look like key=; are parsed correctly as MISSING
                if ("".equals(value)) value = VCFConstants.MISSING_VALUE_v4;

                attributes.put(key, value);
            }
        }

        return attributes;
    }
}
