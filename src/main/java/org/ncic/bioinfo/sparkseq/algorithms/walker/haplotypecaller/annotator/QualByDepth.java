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

import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GATKVariantContextUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.MathUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.RandomGenerator;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefMetaDataTracker;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.ReferenceContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.AlignmentContext;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.PerReadAlleleLikelihoodMap;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.interfaces.ActiveRegionBasedAnnotation;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.interfaces.AnnotatorCompatible;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.interfaces.InfoFieldAnnotation;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.interfaces.StandardAnnotation;

import java.util.*;

/**
 * Variant confidence normalized by unfiltered depth of variant samples
 *
 * <p>This annotation puts the variant confidence QUAL score in perspective by normalizing for the amount of coverage available. Because each read contributes a little to the QUAL score, variants in regions with deep coverage can have artificially inflated QUAL scores, giving the impression that the call is supported by more evidence than it really is. To compensate for this, we normalize the variant confidence by depth, which gives us a more objective picture of how well supported the call is.</p>
 *
 * <h3>Statistical notes</h3>
 * <p>The calculation only takes into account coverage from samples genotyped as having the variant allele(s). This removes the influence of any homozygous-reference samples that might be present in the same cohort, which would otherwise penalize the call unfairly.</p>
 *
 * <h3>Caveats</h3>
 * <p>This annotation can only be calculated for sites for which at least one sample was genotyped as carrying a variant allele.</p>
 *
 * <h3>Related annotations</h3>
 * <ul>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_Coverage.php">Coverage</a></b> gives the filtered depth of coverage for each sample and the unfiltered depth across all samples.</li>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_DepthPerAlleleBySample.php">DepthPerAlleleBySample</a></b> calculates depth of coverage for each allele per sample (AD).</li>
 * </ul>
 */
public class QualByDepth extends InfoFieldAnnotation implements StandardAnnotation, ActiveRegionBasedAnnotation {
//    private final static Logger logger = Logger.getLogger(QualByDepth.class);

    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap ) {
        if ( !vc.hasLog10PError() )
            return null;

        final GenotypesContext genotypes = vc.getGenotypes();
        if ( genotypes == null || genotypes.size() == 0 )
            return null;

        int standardDepth = 0;
        int ADrestrictedDepth = 0;

        for ( final Genotype genotype : genotypes ) {

            // we care only about variant calls with likelihoods
            if ( !genotype.isHet() && !genotype.isHomVar() )
                continue;

            // if we have the AD values for this sample, let's make sure that the variant depth is greater than 1!
            // TODO -- If we like how this is working and want to apply it to a situation other than the single sample HC pipeline,
            // TODO --  then we will need to modify the annotateContext() - and related - routines in the VariantAnnotatorEngine
            // TODO --  so that genotype-level annotations are run first (to generate AD on the samples) and then the site-level
            // TODO --  annotations must come afterwards (so that QD can use the AD).
            if ( genotype.hasAD() ) {
                final int[] AD = genotype.getAD();
                final int totalADdepth = (int) MathUtils.sum(AD);
                if ( totalADdepth - AD[0] > 1 )
                    ADrestrictedDepth += totalADdepth;
                standardDepth += totalADdepth;
                continue;
            }

            if (stratifiedContexts!= null && !stratifiedContexts.isEmpty()) {
                final AlignmentContext context = stratifiedContexts.get(genotype.getSampleName());
                if ( context == null )
                    continue;
                standardDepth += context.getBasePileup().depthOfCoverage();

            } else if (perReadAlleleLikelihoodMap != null) {
                final PerReadAlleleLikelihoodMap perReadAlleleLikelihoods = perReadAlleleLikelihoodMap.get(genotype.getSampleName());
                if (perReadAlleleLikelihoods == null || perReadAlleleLikelihoods.isEmpty())
                    continue;

                standardDepth += perReadAlleleLikelihoods.getNumberOfStoredElements();
            } else if ( genotype.hasDP() ) {
                standardDepth += genotype.getDP();
            }
        }

        // if the AD-restricted depth is a usable value (i.e. not zero), then we should use that one going forward
        if ( ADrestrictedDepth > 0 )
            standardDepth = ADrestrictedDepth;

        if ( standardDepth == 0 )
            return null;

        final double altAlleleLength = GATKVariantContextUtils.getMeanAltAlleleLength(vc);
        // Hack: when refContext == null then we know we are coming from the HaplotypeCaller and do not want to do a
        //  full length-based normalization (because the indel length problem is present only in the UnifiedGenotyper)
        double QD = -10.0 * vc.getLog10PError() / ((double)standardDepth * indelNormalizationFactor(altAlleleLength, ref != null));

        // Hack: see note in the fixTooHighQD method below
        QD = fixTooHighQD(QD);

        final Map<String, Object> map = new HashMap<>();
        map.put(getKeyNames().get(0), String.format("%.2f", QD));
        return map;
    }

    /**
     * Generate the indel normalization factor.
     *
     * @param altAlleleLength  the average alternate allele length for the call
     * @param increaseNormalizationAsLengthIncreases should we apply a normalization factor based on the allele length?
     * @return a possitive double
     */
    private double indelNormalizationFactor(final double altAlleleLength, final boolean increaseNormalizationAsLengthIncreases) {
        return ( increaseNormalizationAsLengthIncreases ? Math.max(altAlleleLength / 3.0, 1.0) : 1.0);
    }

    /**
     * The haplotype caller generates very high quality scores when multiple events are on the
     * same haplotype.  This causes some very good variants to have unusually high QD values,
     * and VQSR will filter these out.  This code looks at the QD value, and if it is above
     * threshold we map it down to the mean high QD value, with some jittering
     *
     * // TODO -- remove me when HaplotypeCaller bubble caller is live
     *
     * @param QD the raw QD score
     * @return a QD value
     */
    private double fixTooHighQD(final double QD) {
        if ( QD < MAX_QD_BEFORE_FIXING ) {
            return QD;
        } else {
            return IDEAL_HIGH_QD + RandomGenerator.getRandomGenerator().nextGaussian() * JITTER_SIGMA;
        }
    }

    private final static double MAX_QD_BEFORE_FIXING = 35;
    private final static double IDEAL_HIGH_QD = 30;
    private final static double JITTER_SIGMA = 3;

    public List<String> getKeyNames() { return Arrays.asList("QD"); }

    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFInfoHeaderLine(getKeyNames().get(0), 1, VCFHeaderLineType.Float, "Variant Confidence/Quality by Depth"));
    }


}
