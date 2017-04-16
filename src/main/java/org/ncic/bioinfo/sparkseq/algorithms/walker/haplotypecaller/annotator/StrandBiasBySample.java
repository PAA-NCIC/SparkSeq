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
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefMetaDataTracker;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.ReferenceContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.AlignmentContext;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.PerReadAlleleLikelihoodMap;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.interfaces.AnnotatorCompatible;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.interfaces.GenotypeAnnotation;

import java.util.Collections;
import java.util.List;

/**
 * Author: wbc
 */
public class StrandBiasBySample extends GenotypeAnnotation {

    public final static String STRAND_BIAS_BY_SAMPLE_KEY_NAME = "SB";

    @Override
    public void annotate(final RefMetaDataTracker tracker,
                         final AnnotatorCompatible walker,
                         final ReferenceContext ref,
                         final AlignmentContext stratifiedContext,
                         final VariantContext vc,
                         final Genotype g,
                         final GenotypeBuilder gb,
                         final PerReadAlleleLikelihoodMap alleleLikelihoodMap) {
        if ( ! isAppropriateInput(alleleLikelihoodMap, g) )
            return;

        final int[][] table = FisherStrand.getContingencyTable(Collections.singletonMap(g.getSampleName(), alleleLikelihoodMap), vc, 0);

        gb.attribute(STRAND_BIAS_BY_SAMPLE_KEY_NAME, FisherStrand.getContingencyArray(table));
    }

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(STRAND_BIAS_BY_SAMPLE_KEY_NAME); }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() { return Collections.singletonList(new VCFFormatHeaderLine(getKeyNames().get(0), 4, VCFHeaderLineType.Integer, "Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.")); }

    private boolean isAppropriateInput(final PerReadAlleleLikelihoodMap map, final Genotype g) {
        return ! (map == null || g == null || !g.isCalled());
    }
}
