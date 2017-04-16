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

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import org.ncic.bioinfo.sparkseq.algorithms.engine.Walker;
import org.ncic.bioinfo.sparkseq.algorithms.utils.MathUtils;
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
 * Likelihood-based test for the inbreeding among samples
 * <p>
 * <p>This annotation estimates whether there is evidence of inbreeding in a population. The higher the score, the higher the chance that there is inbreeding.</p>
 * <p>
 * <h3>Statistical notes</h3>
 * <p>The calculation is a continuous generalization of the Hardy-Weinberg test for disequilibrium that works well with limited coverage per sample. The output is a Phred-scaled p-value derived from running the HW test for disequilibrium with PL values. See the <a href="http://www.broadinstitute.org/gatk/guide/article?id=4732">method document on statistical tests</a> for a more detailed explanation of this statistical test.</p>
 * <p>
 * <h3>Caveats</h3>
 * <h4>Note that the Inbreeding Coefficient can only be calculated for cohorts containing at least 10 founder samples.</h4>
 *
 * Author: wbc
 */
public class InbreedingCoeff extends InfoFieldAnnotation implements StandardAnnotation, ActiveRegionBasedAnnotation {

    private static final int MIN_SAMPLES = 10;
    private static final String INBREEDING_COEFFICIENT_KEY_NAME = "InbreedingCoeff";
    private Set<String> founderIds;
    private int sampleCount;

    @Override
    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {
        //If available, get the founder IDs and cache them. the IC will only be computed on founders then.
        if (founderIds == null && walker != null) {
            founderIds = new HashSet<String>();
            founderIds.add(((Walker) walker).getSampleName());
        }
        return makeCoeffAnnotation(vc);
    }

    protected double calculateIC(final VariantContext vc, final GenotypesContext genotypes) {

        final boolean doMultiallelicMapping = !vc.isBiallelic();

        int idxAA = 0, idxAB = 1, idxBB = 2;

        double refCount = 0.0;
        double hetCount = 0.0;
        double homCount = 0.0;
        sampleCount = 0; // number of samples that have likelihoods

        for (final Genotype g : genotypes) {
            if (g.isCalled() && g.hasLikelihoods() && g.getPloidy() == 2)  // only work for diploid samples
                sampleCount++;
            else
                continue;
            final double[] normalizedLikelihoods = MathUtils.normalizeFromLog10(g.getLikelihoods().getAsVector());
            if (doMultiallelicMapping) {
                if (g.isHetNonRef()) {
                    //all likelihoods go to homCount
                    homCount += 1;
                    continue;
                }

                //get alternate allele for each sample
                final Allele a1 = g.getAllele(0);
                final Allele a2 = g.getAllele(1);
                if (a2.isNonReference()) {
                    final int[] idxVector = vc.getGLIndecesOfAlternateAllele(a2);
                    idxAA = idxVector[0];
                    idxAB = idxVector[1];
                    idxBB = idxVector[2];
                }
                //I expect hets to be reference first, but there are no guarantees (e.g. phasing)
                else if (a1.isNonReference()) {
                    final int[] idxVector = vc.getGLIndecesOfAlternateAllele(a1);
                    idxAA = idxVector[0];
                    idxAB = idxVector[1];
                    idxBB = idxVector[2];
                }
            }

            refCount += normalizedLikelihoods[idxAA];
            hetCount += normalizedLikelihoods[idxAB];
            homCount += normalizedLikelihoods[idxBB];
        }

        final double p = (2.0 * refCount + hetCount) / (2.0 * (refCount + hetCount + homCount)); // expected reference allele frequency
        final double q = 1.0 - p; // expected alternative allele frequency
        final double F = 1.0 - (hetCount / (2.0 * p * q * (double) sampleCount)); // inbreeding coefficient

        return F;
    }

    protected Map<String, Object> makeCoeffAnnotation(final VariantContext vc) {
        final GenotypesContext genotypes = (founderIds == null || founderIds.isEmpty()) ? vc.getGenotypes() : vc.getGenotypes(founderIds);
        if (genotypes == null || genotypes.size() < MIN_SAMPLES || !vc.isVariant())
            return null;
        double F = calculateIC(vc, genotypes);
        if (sampleCount < MIN_SAMPLES)
            return null;
        return Collections.singletonMap(getKeyNames().get(0), (Object) String.format("%.4f", F));
    }

    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(INBREEDING_COEFFICIENT_KEY_NAME);
    }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Collections.singletonList(new VCFInfoHeaderLine(INBREEDING_COEFFICIENT_KEY_NAME, 1, VCFHeaderLineType.Float, "Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation"));
    }
}