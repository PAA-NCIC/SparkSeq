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
package org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.afcalculate;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.ncic.bioinfo.sparkseq.algorithms.utils.MathUtils;
import org.ncic.bioinfo.sparkseq.algorithms.data.basic.Pair;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Map;

/**
 * Author: wbc
 */
class OriginalDiploidExactAFCalculator extends DiploidExactAFCalculator {
    protected OriginalDiploidExactAFCalculator() {
    }

    @Override
    protected AFCalculationResult computeLog10PNonRef(final VariantContext vc,
                                                      @SuppressWarnings("unused")
                                                      final int defaultPloidy,
                                                      final double[] log10AlleleFrequencyPriors,
                                                      final StateTracker stateTracker) {
        final double[] log10AlleleFrequencyLikelihoods = new double[log10AlleleFrequencyPriors.length];
        final double[] log10AlleleFrequencyPosteriors  = new double[log10AlleleFrequencyPriors.length];
        final Pair<Integer, Integer> result = linearExact(vc, log10AlleleFrequencyPriors, log10AlleleFrequencyLikelihoods, log10AlleleFrequencyPosteriors);
        final int lastK = result.getFirst();
        final int mleK = result.getSecond();

        final double log10LikelihoodAFGt0 = lastK == 0 ? MathUtils.LOG10_P_OF_ZERO : MathUtils.log10sumLog10(log10AlleleFrequencyLikelihoods, 1, lastK+1);
        final double[] log10Likelihoods = new double[]{log10AlleleFrequencyLikelihoods[0], log10LikelihoodAFGt0};
        final double[] log10Priors = new double[]{log10AlleleFrequencyPriors[0], MathUtils.log10sumLog10(log10AlleleFrequencyPriors, 1)};
        final double[] log10Posteriors = MathUtils.vectorSum(log10Likelihoods, log10Priors);

        final double log10PRef = log10Posteriors[1] > log10Posteriors[0] ? MathUtils.LOG10_P_OF_ZERO : 0.0;
        final Map<Allele, Double> log10pRefByAllele = Collections.singletonMap(vc.getAlternateAllele(0), log10PRef);

        return new AFCalculationResult(new int[]{mleK}, 0, vc.getAlleles(),
                MathUtils.normalizeFromLog10(log10Likelihoods, true),
                MathUtils.normalizeFromLog10(log10Priors, true),
                log10pRefByAllele);
    }

    /**
     * A simple data structure that holds the current, prev, and prev->prev likelihoods vectors
     * for the exact model calculation
     */
    private final static class ExactACCache {
        double[] kMinus2, kMinus1, kMinus0;

        private static double[] create(int n) {
            return new double[n];
        }

        public ExactACCache(int n) {
            kMinus2 = create(n);
            kMinus1 = create(n);
            kMinus0 = create(n);
        }

        final public void rotate() {
            double[] tmp = kMinus2;
            kMinus2 = kMinus1;
            kMinus1 = kMinus0;
            kMinus0 = tmp;
        }

        final public double[] getkMinus2() {
            return kMinus2;
        }

        final public double[] getkMinus1() {
            return kMinus1;
        }

        final public double[] getkMinus0() {
            return kMinus0;
        }
    }

    private Pair<Integer, Integer> linearExact(final VariantContext vc,
                                               double[] log10AlleleFrequencyPriors,
                                               double[] log10AlleleFrequencyLikelihoods,
                                               double[] log10AlleleFrequencyPosteriors) {
        final ArrayList<double[]> genotypeLikelihoods = getGLs(vc.getGenotypes(), true);
        final int numSamples = genotypeLikelihoods.size()-1;
        final int numChr = 2*numSamples;

        final ExactACCache logY = new ExactACCache(numSamples+1);
        logY.getkMinus0()[0] = 0.0; // the zero case

        double maxLog10L = Double.NEGATIVE_INFINITY;
        boolean done = false;
        int lastK = -1, mleK = -1;

        for (int k=0; k <= numChr && ! done; k++ ) {
            final double[] kMinus0 = logY.getkMinus0();

            if ( k == 0 ) { // special case for k = 0
                for ( int j=1; j <= numSamples; j++ ) {
                    kMinus0[j] = kMinus0[j-1] + genotypeLikelihoods.get(j)[0];
                }
            } else { // k > 0
                final double[] kMinus1 = logY.getkMinus1();
                final double[] kMinus2 = logY.getkMinus2();

                for ( int j=1; j <= numSamples; j++ ) {
                    final double[] gl = genotypeLikelihoods.get(j);
                    final double logDenominator = MathUtils.Log10Cache.get(2*j) + MathUtils.Log10Cache.get(2*j-1);

                    double aa = Double.NEGATIVE_INFINITY;
                    double ab = Double.NEGATIVE_INFINITY;
                    if (k < 2*j-1)
                        aa = MathUtils.Log10Cache.get(2*j-k) + MathUtils.Log10Cache.get(2*j-k-1) + kMinus0[j-1] + gl[0];

                    if (k < 2*j)
                        ab = MathUtils.Log10Cache.get(2*k) + MathUtils.Log10Cache.get(2*j-k)+ kMinus1[j-1] + gl[1];

                    double log10Max;
                    if (k > 1) {
                        final double bb = MathUtils.Log10Cache.get(k) + MathUtils.Log10Cache.get(k-1) + kMinus2[j-1] + gl[2];
                        log10Max = MathUtils.approximateLog10SumLog10(aa, ab, bb);
                    } else {
                        // we know we aren't considering the BB case, so we can use an optimized log10 function
                        log10Max = MathUtils.approximateLog10SumLog10(aa, ab);
                    }

                    // finally, update the L(j,k) value
                    kMinus0[j] = log10Max - logDenominator;
                }
            }

            // update the posteriors vector
            final double log10LofK = kMinus0[numSamples];
            log10AlleleFrequencyLikelihoods[k] = log10LofK;
            log10AlleleFrequencyPosteriors[k] = log10LofK + log10AlleleFrequencyPriors[k];

            // can we abort early?
            lastK = k;
            if ( log10LofK > maxLog10L ) {
                maxLog10L = log10LofK;
                mleK = k;
            }

            if ( log10LofK < maxLog10L - StateTracker.MAX_LOG10_ERROR_TO_STOP_EARLY ) {
                //if ( DEBUG ) System.out.printf("  *** breaking early k=%d log10L=%.2f maxLog10L=%.2f%n", k, log10LofK, maxLog10L);
                done = true;
            }

            logY.rotate();
        }

        return new Pair<>(lastK, mleK);
    }
}
