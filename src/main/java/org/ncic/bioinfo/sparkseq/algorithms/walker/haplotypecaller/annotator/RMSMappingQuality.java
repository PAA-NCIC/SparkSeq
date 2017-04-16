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

import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import htsjdk.variant.variantcontext.VariantContext;
import org.ncic.bioinfo.sparkseq.algorithms.utils.MathUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.QualityUtils;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefMetaDataTracker;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.ReferenceContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.AlignmentContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.GATKSAMRecord;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.PileupElement;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.PerReadAlleleLikelihoodMap;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.interfaces.ActiveRegionBasedAnnotation;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.interfaces.AnnotatorCompatible;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.interfaces.InfoFieldAnnotation;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.interfaces.StandardAnnotation;

import java.util.*;

/**
 * Root Mean Square of the mapping quality of reads across all samples.
 *
 * <p>This annotation provides an estimation of the overall mapping quality of reads supporting a variant call, averaged over all samples in a cohort.</p>
 *
 * <h3>Statistical notes</h3>
 * <p>The root mean square is equivalent to the mean of the mapping qualities plus the standard deviation of the mapping qualities.</p>
 *
 * <h3>Related annotations</h3>
 * <ul>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_MappingQualityRankSumTest.php">MappingQualityRankSumTest</a></b> compares the mapping quality of reads supporting the REF and ALT alleles.</li>
 * </ul>
 *
 */
public class RMSMappingQuality extends InfoFieldAnnotation implements StandardAnnotation, ActiveRegionBasedAnnotation {

    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap ) {

        final List<Integer> qualities = new ArrayList<>();
        if ( stratifiedContexts != null ) {
            if ( stratifiedContexts.size() == 0 )
                return null;

            for ( final Map.Entry<String, AlignmentContext> sample : stratifiedContexts.entrySet() ) {
                final AlignmentContext context = sample.getValue();
                for ( final PileupElement p : context.getBasePileup() )
                    fillMappingQualitiesFromPileup(p.getRead().getMappingQuality(), qualities);
            }
        }
        else if (perReadAlleleLikelihoodMap != null) {
            if ( perReadAlleleLikelihoodMap.size() == 0 )
                return null;

            for ( final PerReadAlleleLikelihoodMap perReadLikelihoods : perReadAlleleLikelihoodMap.values() ) {
                for ( final GATKSAMRecord read : perReadLikelihoods.getStoredElements() )
                    fillMappingQualitiesFromPileup(read.getMappingQuality(), qualities);
            }
        }
        else
            return null;

        final double rms = MathUtils.rms(qualities);
        return Collections.singletonMap(getKeyNames().get(0), (Object)String.format("%.2f", rms));
    }

    private static void fillMappingQualitiesFromPileup(final int mq, final List<Integer> qualities) {
        if ( mq != QualityUtils.MAPPING_QUALITY_UNAVAILABLE ) {
            qualities.add(mq);
        }
    }

    public List<String> getKeyNames() { return Arrays.asList(VCFConstants.RMS_MAPPING_QUALITY_KEY); }

    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(VCFStandardHeaderLines.getInfoLine(getKeyNames().get(0)));
    }
}