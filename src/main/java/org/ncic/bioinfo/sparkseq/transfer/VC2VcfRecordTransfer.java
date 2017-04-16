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
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.LazyGenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.IntGenotypeFieldAccessors;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.ncic.bioinfo.sparkseq.data.basic.VcfRecord;
import org.ncic.bioinfo.sparkseq.data.common.RefContigInfo;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

/**
 * Author: wbc
 */
public class VC2VcfRecordTransfer {
    private static final String QUAL_FORMAT_STRING = "%.2f";
    private static final String QUAL_FORMAT_EXTENSION_TO_TRIM = ".00";

    private final IntGenotypeFieldAccessors GENOTYPE_FIELD_ACCESSORS = new IntGenotypeFieldAccessors();

    private VCFHeader header;
    private RefContigInfo refContigInfo;

    public VC2VcfRecordTransfer(VCFHeader header, RefContigInfo refContigInfo) {
        this.header = header;
        this.refContigInfo = refContigInfo;
    }

    public VcfRecord transfer(VariantContext variantContext) {
        return transfer(variantContext, true);
    }

    public VcfRecord transfer(VariantContext variantContext, boolean withGenotype) {
        String contigName = variantContext.getChr();
        int contigID= refContigInfo.getId(contigName);
        int position = variantContext.getStart();
        String id = variantContext.getID();
        String ref = variantContext.getReference().getDisplayString();

        StringBuilder altBuilder = new StringBuilder();
        if (variantContext.isVariant()) {
            Allele altAllele = variantContext.getAlternateAllele(0);
            String altStr = altAllele.getDisplayString();
            altBuilder.append(altStr);

            for (int i = 1; i < variantContext.getAlternateAlleles().size(); i++) {
                altAllele = variantContext.getAlternateAllele(i);
                altStr = altAllele.getDisplayString();
                altBuilder.append(",");
                altBuilder.append(altStr);
            }
        } else {
            altBuilder.append(VCFConstants.EMPTY_ALTERNATE_ALLELE_FIELD);
        }
        String alt = altBuilder.toString();

        double qual = (!variantContext.hasLog10PError())?  0 : variantContext.getPhredScaledQual();

        String filter = getFilterString(variantContext);

        StringBuilder rest = new StringBuilder();
        // INFO
        final Map<String, String> infoFields = new TreeMap<String, String>();
        for (final Map.Entry<String, Object> field : variantContext.getAttributes().entrySet()) {
            final String outputValue = formatVCFField(field.getValue());
            if (outputValue != null) infoFields.put(field.getKey(), outputValue);
        }
        writeInfoString(infoFields, rest);

        // FORMAT
        if(withGenotype) {
            final GenotypesContext gc = variantContext.getGenotypes();
            if (gc.isLazyWithData() && ((LazyGenotypesContext) gc).getUnparsedGenotypeData() instanceof String) {
                rest.append(VCFConstants.FIELD_SEPARATOR);
                rest.append(((LazyGenotypesContext) gc).getUnparsedGenotypeData().toString());
            } else {
                final List<String> genotypeAttributeKeys = variantContext.calcVCFGenotypeKeys(this.header);
                if (!genotypeAttributeKeys.isEmpty()) {
                    final String genotypeFormatString = ParsingUtils.join(VCFConstants.GENOTYPE_FIELD_SEPARATOR, genotypeAttributeKeys);

                    rest.append(VCFConstants.FIELD_SEPARATOR);
                    rest.append(genotypeFormatString);

                    final Map<Allele, String> alleleStrings = buildAlleleStrings(variantContext);
                    addGenotypeData(variantContext, alleleStrings, genotypeAttributeKeys, rest);
                }
            }
        }

        return new VcfRecord(contigName, contigID, position, id, ref, alt, qual, filter, rest.toString());
    }

    private String getFilterString(final VariantContext vc) {
        if (vc.isFiltered()) {
            return ParsingUtils.join(";", ParsingUtils.sortList(vc.getFilters()));
        } else if (vc.filtersWereApplied()) return VCFConstants.PASSES_FILTERS_v4;
        else return VCFConstants.UNFILTERED;
    }

    String formatVCFField(final Object val) {
        final String result;
        if (val == null)
            result = VCFConstants.MISSING_VALUE_v4;
        else if (val instanceof Double)
            result = formatVCFDouble((Double) val);
        else if (val instanceof Boolean)
            result = (Boolean) val ? "" : null; // empty string for true, null for false
        else if (val instanceof List) {
            result = formatVCFField(((List) val).toArray());
        } else if (val.getClass().isArray()) {
            final int length = Array.getLength(val);
            if (length == 0)
                return formatVCFField(null);
            final StringBuilder sb = new StringBuilder(formatVCFField(Array.get(val, 0)));
            for (int i = 1; i < length; i++) {
                sb.append(",");
                sb.append(formatVCFField(Array.get(val, i)));
            }
            result = sb.toString();
        } else
            result = val.toString();

        return result;
    }

    public static String formatVCFDouble(final double d) {
        final String format;
        if (d < 1) {
            if (d < 0.01) {
                if (Math.abs(d) >= 1e-20)
                    format = "%.3e";
                else {
                    // return a zero format
                    return "0.00";
                }
            } else {
                format = "%.3f";
            }
        } else {
            format = "%.2f";
        }

        return String.format(format, d);
    }

    private void writeInfoString(final Map<String, String> infoFields, final StringBuilder builder) {
        if (infoFields.isEmpty()) {
            builder.append(VCFConstants.EMPTY_INFO_FIELD);
            return;
        }

        boolean isFirst = true;
        for (final Map.Entry<String, String> entry : infoFields.entrySet()) {
            if (isFirst) isFirst = false;
            else builder.append(VCFConstants.INFO_FIELD_SEPARATOR);

            builder.append(entry.getKey());

            if (!entry.getValue().equals("")) {
                final VCFInfoHeaderLine metaData = this.header.getInfoHeaderLine(entry.getKey());
                if (metaData == null || metaData.getCountType() != VCFHeaderLineCount.INTEGER || metaData.getCount() != 0) {
                    builder.append("=");
                    builder.append(entry.getValue());
                }
            }
        }
    }

    public Map<Allele, String> buildAlleleStrings(final VariantContext vc) {
        final Map<Allele, String> alleleMap = new HashMap<Allele, String>(vc.getAlleles().size()+1);
        alleleMap.put(Allele.NO_CALL, VCFConstants.EMPTY_ALLELE); // convenience for lookup

        final List<Allele> alleles = vc.getAlleles();
        for ( int i = 0; i < alleles.size(); i++ ) {
            alleleMap.put(alleles.get(i), String.valueOf(i));
        }

        return alleleMap;
    }

    public void addGenotypeData(final VariantContext vc, final Map<Allele, String> alleleMap, final List<String> genotypeFormatKeys, final StringBuilder builder) {
        final int ploidy = vc.getMaxPloidy(2);

        for (final String sample : this.header.getGenotypeSamples()) {
            builder.append(VCFConstants.FIELD_SEPARATOR);

            Genotype g = vc.getGenotype(sample);
            if (g == null) g = GenotypeBuilder.createMissing(sample, ploidy);

            final List<String> attrs = new ArrayList<String>(genotypeFormatKeys.size());
            for (final String field : genotypeFormatKeys) {
                if (field.equals(VCFConstants.GENOTYPE_KEY)) {
                    if ( ! g.isAvailable()) {
                        throw new IllegalStateException("GTs cannot be missing for some samples if they are available for others in the record");
                    }

                    writeAllele(g.getAllele(0), alleleMap, builder);
                    for (int i = 1; i < g.getPloidy(); i++) {
                        builder.append(g.isPhased() ? VCFConstants.PHASED : VCFConstants.UNPHASED);
                        writeAllele(g.getAllele(i), alleleMap, builder);
                    }
                    continue;

                } else {
                    final String outputValue;
                    if ( field.equals(VCFConstants.GENOTYPE_FILTER_KEY ) ) {
                        outputValue = g.isFiltered() ? g.getFilters() : VCFConstants.PASSES_FILTERS_v4;
                    } else {
                        final IntGenotypeFieldAccessors.Accessor accessor = GENOTYPE_FIELD_ACCESSORS.getAccessor(field);
                        if ( accessor != null ) {
                            final int[] intValues = accessor.getValues(g);
                            if ( intValues == null )
                                outputValue = VCFConstants.MISSING_VALUE_v4;
                            else if ( intValues.length == 1 ) // fast path
                                outputValue = Integer.toString(intValues[0]);
                            else {
                                final StringBuilder sb = new StringBuilder();
                                sb.append(intValues[0]);
                                for ( int i = 1; i < intValues.length; i++) {
                                    sb.append(",");
                                    sb.append(intValues[i]);
                                }
                                outputValue = sb.toString();
                            }
                        } else {
                            Object val = g.hasExtendedAttribute(field) ? g.getExtendedAttribute(field) : VCFConstants.MISSING_VALUE_v4;

                            final VCFFormatHeaderLine metaData = this.header.getFormatHeaderLine(field);
                            if ( metaData != null ) {
                                final int numInFormatField = metaData.getCount(vc);
                                if ( numInFormatField > 1 && val.equals(VCFConstants.MISSING_VALUE_v4) ) {
                                    // If we have a missing field but multiple values are expected, we need to construct a new string with all fields.
                                    // For example, if Number=2, the string has to be ".,."
                                    final StringBuilder sb = new StringBuilder(VCFConstants.MISSING_VALUE_v4);
                                    for ( int i = 1; i < numInFormatField; i++ ) {
                                        sb.append(",");
                                        sb.append(VCFConstants.MISSING_VALUE_v4);
                                    }
                                    val = sb.toString();
                                }
                            }

                            // assume that if key is absent, then the given string encoding suffices
                            outputValue = formatVCFField(val);
                        }
                    }

                    if ( outputValue != null )
                        attrs.add(outputValue);
                }
            }

            // strip off trailing missing values
            for (int i = attrs.size() - 1; i >= 0; i--) {
                if (isMissingValue(attrs.get(i))) attrs.remove(i);
                else break;
            }

            for (int i = 0; i < attrs.size(); i++) {
                if ( i > 0 || genotypeFormatKeys.contains(VCFConstants.GENOTYPE_KEY)) {
                    builder.append(VCFConstants.GENOTYPE_FIELD_SEPARATOR);
                }
                builder.append(attrs.get(i));
            }
        }
    }

    private void writeAllele(final Allele allele, final Map<Allele, String> alleleMap, final StringBuilder builder) {
        final String encoding = alleleMap.get(allele);
        if ( encoding == null )
            throw new RuntimeException("Allele " + allele + " is not an allele in the variant context");
        builder.append(encoding);
    }

    static int countOccurrences(final char c, final String s) {
        int count = 0;
        for (int i = 0; i < s.length(); i++) {
            count += s.charAt(i) == c ? 1 : 0;
        }
        return count;
    }

    static boolean isMissingValue(final String s) {
        // we need to deal with the case that it's a list of missing values
        return (countOccurrences(VCFConstants.MISSING_VALUE_v4.charAt(0), s) + countOccurrences(',', s) == s.length());
    }
}
