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
package org.ncic.bioinfo.sparkseq.algorithms.walker.mutect;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.BinomialDistribution;
import org.apache.commons.math.distribution.BinomialDistributionImpl;

/**
 * Author: wbc
 */
public class TumorPowerCalculator extends AbstractPowerCalculator{
    private double constantContamination;
    private boolean enableSmoothing;

    public TumorPowerCalculator(double constantEps, double constantLodThreshold, double constantContamination) {
        this(constantEps, constantLodThreshold, constantContamination, true);
    }

    public TumorPowerCalculator(double constantEps, double constantLodThreshold, double constantContamination, boolean enableSmoothing) {
        this.constantEps = constantEps;
        this.constantLodThreshold = constantLodThreshold;
        this.constantContamination = constantContamination;
        this.enableSmoothing = enableSmoothing;
    }

    public double cachingPowerCalculation(int n, double delta) throws MathException {
        PowerCacheKey key = new PowerCacheKey(n, delta);
        Double power = cache.get(key);
        if (power == null) {
            power = calculatePower(n, constantEps, constantLodThreshold, delta, constantContamination, enableSmoothing);
            cache.put(key, power);
        }
        return power;
    }




    protected static double calculateTumorLod(int depth, int alts, double eps, double contam) {
        double f = (double) alts / (double) depth;
        return (calculateLogLikelihood(depth, alts, eps, f) - calculateLogLikelihood(depth, alts, eps, Math.min(f,contam)));
    }

    protected static double calculatePower(int depth, double eps, double lodThreshold, double delta, double contam, boolean enableSmoothing) throws MathException {
        if (depth==0) return 0;

        // calculate the probability of each configuration
        double p_alt_given_e_delta = delta*(1d-eps) + (1d-delta)*eps;
        BinomialDistribution binom = new BinomialDistributionImpl(depth, p_alt_given_e_delta);
        double[] p = new double[depth+1];
        for(int i=0; i<p.length; i++) {
            p[i] = binom.probability(i);
        }

        // calculate the LOD scores
        double[] lod = new double[depth+1];
        for(int i=0; i<lod.length; i++) {
            lod[i] = calculateTumorLod(depth, i, eps, contam);
        }

        int k = -1;
        for(int i=0; i<lod.length; i++) {
            if (lod[i] >= lodThreshold) {
                k = i;
                break;
            }
        }

        // if no depth meets the lod score, the power is zero
        if (k == -1) {
            return 0;
        }

        double power = 0;

        // here we correct for the fact that the exact lod threshold is likely somewhere between
        // the k and k-1 bin, so we prorate the power from that bin
        // the k and k-1 bin, so we prorate the power from that bin
        // if k==0, it must be that lodThreshold == lod[k] so we don't have to make this correction
        if ( enableSmoothing && k > 0 ) {
            double x = 1d - (lodThreshold - lod[k-1]) / (lod[k] - lod[k-1]);
            power = x*p[k-1];
        }

        for(int i=k; i<p.length; i++) {
            power += p[i];
        }

        return(power);
    }

}
