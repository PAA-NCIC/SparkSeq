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
package org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.genotyper;

import htsjdk.variant.variantcontext.VariantContext;
import org.ncic.bioinfo.sparkseq.algorithms.utils.AlleleListUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.MathUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.SampleListUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.haplotype.Haplotype;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.AlleleList;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.readlikelihood.PairHMMLikelihoodCalculationEngine;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.readlikelihood.ReadLikelihoods;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Author: wbc
 */
public class HaplotypeLDCalculator {
    private final List<Haplotype> haplotypes;
    private final ReadLikelihoods<Haplotype> readLikelihoods;
    private List<Map<Haplotype, Double>> haplotypeLikelihoodsPerSample = null;

    // linear contigency table with table[0] == [0][0], table[1] = [0][1], table[2] = [1][0], table[3] = [1][1]
    private final double[] table = new double[4];

    public HaplotypeLDCalculator(final List<Haplotype> haplotypes, final ReadLikelihoods<Haplotype> haplotypeReadMap) {
        this.haplotypes = haplotypes;
        this.readLikelihoods = haplotypeReadMap;
    }

    /**
     * Construct the cached list of summed haplotype likelihoods per sample if it
     * hasn't already been computed.  This data structure is lazy created but only
     * needs to be made once when we make 1 merge decision as the data doesn't change
     * no matter how many calls to computeProbOfBeingPhased
     */
    private void buildHaplotypeLikelihoodsPerSampleIfNecessary() {
        if (haplotypeLikelihoodsPerSample == null) {
            // do the lazy computation
            final Set<String> samples = new LinkedHashSet<>(readLikelihoods.samples());
            haplotypeLikelihoodsPerSample = new LinkedList<>();
            for (final String sample : samples) {
                final Map<Haplotype, Double> map = new HashMap<>(haplotypes.size());
                for (final Haplotype h : haplotypes) {
                    // count up the co-occurrences of the events for the R^2 calculation
                    final double haplotypeLikelihood =
                            PairHMMLikelihoodCalculationEngine.computeDiploidHaplotypeLikelihoods(
                                    sample, readLikelihoods, Collections.singletonList(h), false)[0][0];
                    map.put(h, haplotypeLikelihood);
                }
                haplotypeLikelihoodsPerSample.add(map);
            }
        }
    }

    /**
     * Compute the likelihood based probability that that haplotypes for first and second are only x11 and x22
     * <p>
     * As opposed to the hypothesis that all four haplotypes (x11, x12, x21, and x22) exist in the population
     *
     * @param first  a non-null VariantContext
     * @param second a non-null VariantContext
     * @return the probability that only x11 and x22 exist among the samples
     */
    protected double computeProbOfBeingPhased(final VariantContext first, final VariantContext second) {
        buildHaplotypeLikelihoodsPerSampleIfNecessary();

        Arrays.fill(table, Double.NEGATIVE_INFINITY);

        for (final Map<Haplotype, Double> entry : haplotypeLikelihoodsPerSample) {
            for (final Map.Entry<Haplotype, Double> haplotypeLikelihood : entry.entrySet()) {
                final Haplotype h = haplotypeLikelihood.getKey();
                // count up the co-occurrences of the events for the R^2 calculation
                final VariantContext thisHapVC = h.getEventMap().get(first.getStart());
                final VariantContext nextHapVC = h.getEventMap().get(second.getStart()); // TODO -- add function to take a VC
                final int i = thisHapVC == null ? 0 : 1;
                final int j = nextHapVC == null ? 0 : 1;
                final int index = 2 * i + j;
                table[index] = MathUtils.approximateLog10SumLog10(table[index], haplotypeLikelihood.getValue());
            }
        }

        return pPhased(table);
    }

    /**
     * Compute probability that two variants are in phase with each other and that no
     * compound hets exist in the population.
     * <p>
     * Implemented as a likelihood ratio test of the hypothesis:
     * <p>
     * x11 and x22 are the only haplotypes in the populations
     * <p>
     * vs.
     * <p>
     * all four haplotype combinations (x11, x12, x21, and x22) all exist in the population.
     * <p>
     * Now, since we have to have both variants in the population, we exclude the x11 & x11 state.  So the
     * p of having just x11 and x22 is P(x11 & x22) + p(x22 & x22).
     * <p>
     * Alternatively, we might have any configuration that gives us both 1 and 2 alts, which are:
     * <p>
     * - P(x11 & x12 & x21) -- we have hom-ref and both hets
     * - P(x22 & x12 & x21) -- we have hom-alt and both hets
     * - P(x22 & x12) -- one haplotype is 22 and the other is het 12
     * - P(x22 & x21) -- one haplotype is 22 and the other is het 21
     * <p>
     * The probability is just p11_22 / (p11_22 + p hets)
     *
     * @param table linear contigency table with table[0] == [0][0], table[1] = [0][1], table[2] = [1][0], table[3] = [1][1]
     *              doesn't have to be normalized as this function does the normalization internally
     * @return the real space probability that the data is phased
     */
    protected double pPhased(double[] table) {
        final double[] normTable = MathUtils.normalizeFromLog10(table, true);

        final double x11 = normTable[0], x12 = normTable[1], x21 = normTable[2], x22 = normTable[3];

        // probability that we are only x11 && x22
        final double p11_22 = MathUtils.approximateLog10SumLog10(x11 + x22, x22 + x22);

        // probability of having any of the other pairs
        final double p11_12_21 = MathUtils.approximateLog10SumLog10(x11 + x12, x11 + x21, x12 + x21);
        final double p22_12_21 = MathUtils.approximateLog10SumLog10(x22 + x12, x22 + x21, x12 + x21);
        final double p22_12 = x22 + x12;
        final double p22_21 = x22 + x21;
        final double pOthers = MathUtils.approximateLog10SumLog10(new double[]{p11_12_21, p22_12_21, p22_12, p22_21});

        // probability of being phases is the ratio of p11_22 / pOthers which in log space is just a substraction
        final double log10phased = p11_22 - (MathUtils.approximateLog10SumLog10(p11_22, pOthers));

        return Math.pow(10.0, log10phased);
    }

    protected double pPhasedTest(final double x11, final double x12, final double x21, final double x22) {
        return pPhased(new double[]{x11, x12, x21, x22});
    }
}
