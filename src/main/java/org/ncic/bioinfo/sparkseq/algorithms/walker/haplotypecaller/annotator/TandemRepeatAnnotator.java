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
package org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator;

import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.variantcontext.VariantContext;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GATKVariantContextUtils;
import org.ncic.bioinfo.sparkseq.algorithms.data.basic.Pair;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefMetaDataTracker;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.ReferenceContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.AlignmentContext;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.PerReadAlleleLikelihoodMap;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.interfaces.AnnotatorCompatible;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.interfaces.InfoFieldAnnotation;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.interfaces.StandardAnnotation;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Tandem repeat unit composition and counts per allele
 *
 * <p>This annotation tags variants that fall within tandem repeat sets. It also provides the composition of the tandem repeat units and the number of times they are repeated for each allele (including the REF allele).</p>
 *
 * <p>A tandem repeat unit is composed of one or more nucleotides that are repeated multiple times in series. Repetitive sequences are difficult to map to the reference because they are associated with multiple alignment possibilities. Knowing the number of repeat units in a set of tandem repeats tells you the number of different positions the tandem repeat can be placed in. The observation of many tandem repeat units multiplies the number of possible representations that can be made of the region.
 *
 * <h3>Caveats</h3>
 * <ul>
 *     <li>This annotation is currently not compatible with HaplotypeCaller.</li>
 * </ul>
 *
 */
public class TandemRepeatAnnotator extends InfoFieldAnnotation implements StandardAnnotation {
    private static final String STR_PRESENT = "STR";
    private static final String REPEAT_UNIT_KEY = "RU";
    private static final String REPEATS_PER_ALLELE_KEY = "RPA";
    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {
        if ( !vc.isIndel())
            return null;

        Pair<List<Integer>,byte[]> result = GATKVariantContextUtils.getNumTandemRepeatUnits(vc, ref.getForwardBases());
        if (result == null)
            return null;

        byte[] repeatUnit = result.second;
        List<Integer> numUnits = result.first;

        Map<String, Object> map = new HashMap<String, Object>();
        map.put(STR_PRESENT,true);
        map.put(REPEAT_UNIT_KEY,new String(repeatUnit));
        map.put(REPEATS_PER_ALLELE_KEY, numUnits);

        return map;
    }

    protected static final String[] keyNames = {STR_PRESENT, REPEAT_UNIT_KEY,REPEATS_PER_ALLELE_KEY };
    protected static final VCFInfoHeaderLine[] descriptions = {
            new VCFInfoHeaderLine(STR_PRESENT, 0, VCFHeaderLineType.Flag, "Variant is a short tandem repeat"),
            new VCFInfoHeaderLine(REPEAT_UNIT_KEY, 1, VCFHeaderLineType.String, "Tandem repeat unit (bases)"),
            new VCFInfoHeaderLine(REPEATS_PER_ALLELE_KEY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Number of times tandem repeat unit is repeated, for each allele (including reference)") };

    public List<String> getKeyNames() {
        return Arrays.asList(keyNames);
    }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(descriptions); }

}
