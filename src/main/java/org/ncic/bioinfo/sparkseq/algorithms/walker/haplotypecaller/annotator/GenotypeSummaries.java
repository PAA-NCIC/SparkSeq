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

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefMetaDataTracker;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.ReferenceContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.AlignmentContext;
import org.ncic.bioinfo.sparkseq.algorithms.utils.MathUtils;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.PerReadAlleleLikelihoodMap;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.interfaces.ActiveRegionBasedAnnotation;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.interfaces.AnnotatorCompatible;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.interfaces.InfoFieldAnnotation;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Author: wbc
 */
public class GenotypeSummaries extends InfoFieldAnnotation implements ActiveRegionBasedAnnotation {

    public final static String CCC = "CCC";
    public final static String NCC = "NCC";
    public final static String HWP = "HWP";
    public final static String GQ_MEAN = "GQ_MEAN";
    public final static String GQ_STDDEV = "GQ_STDDEV";

    @Override
    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap ) {
        if ( ! vc.hasGenotypes() )
            return null;

        final Map<String,Object> returnMap = new HashMap<>();
        returnMap.put(NCC, vc.getNoCallCount());

        final MathUtils.RunningAverage average = new MathUtils.RunningAverage();
        for( final Genotype g : vc.getGenotypes() ) {
            if( g.hasGQ() ) {
                average.add(g.getGQ());
            }
        }
        if( average.observationCount() > 0L ) {
            returnMap.put(GQ_MEAN, String.format("%.2f", average.mean()));
            if( average.observationCount() > 1L ) {
                returnMap.put(GQ_STDDEV, String.format("%.2f", average.stddev()));
            }
        }

        return returnMap;
    }

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(CCC, NCC, HWP, GQ_MEAN, GQ_STDDEV);
    }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(
                new VCFInfoHeaderLine(CCC, 1, VCFHeaderLineType.Integer, "Number of called chromosomes"),
                new VCFInfoHeaderLine(NCC, 1, VCFHeaderLineType.Integer, "Number of no-called samples"),
                new VCFInfoHeaderLine(HWP, 1, VCFHeaderLineType.Float, "P value from test of Hardy Weinberg Equilibrium"),
                new VCFInfoHeaderLine(GQ_MEAN, 1, VCFHeaderLineType.Float, "Mean of all GQ values"),
                new VCFInfoHeaderLine(GQ_STDDEV, 1, VCFHeaderLineType.Float, "Standard deviation of all GQ values")
        );
    }
}
