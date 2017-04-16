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

import org.ncic.bioinfo.sparkseq.algorithms.utils.MathUtils;

import java.util.Arrays;

/**
 * Author: wbc
 */
public class HeterozygosityAFPriorProvider extends AFPriorProvider {

    private final double heterozygosity;
    private final double log10Heterozygosity;

    /**
     * Construct a new provider given the heterozygosity value.
     * @param heterozygosity must be a valid heterozygosity between larger than 0 and smaller than 1.
     * @throws IllegalArgumentException if {@code heterozygosity} is not valid one in the interval (0,1).
     */
    public HeterozygosityAFPriorProvider(final double heterozygosity) {
        if (heterozygosity <= 0)
            throw new IllegalArgumentException("the heterozygosity must be greater than 0");
        if (heterozygosity >= 1)
            throw new IllegalArgumentException("the heterozygosity must be less than 1");
        if (Double.isNaN(heterozygosity))
            throw new IllegalArgumentException("the heterozygosity cannot be a NaN");
        this.heterozygosity = heterozygosity;
        this.log10Heterozygosity = Math.log10(heterozygosity);
    }

    @Override
    protected double[] buildPriors(final int totalPloidy) {
        final double[] result = new double [totalPloidy + 1];
        Arrays.fill(result, log10Heterozygosity);
        result[0] = Double.NEGATIVE_INFINITY;
        MathUtils.Log10Cache.ensureCacheContains(totalPloidy);
        for (int i = 1; i <= totalPloidy; i++)
            result[i] -= MathUtils.Log10Cache.get(i);
        final double log10Sum = MathUtils.approximateLog10SumLog10(result);
        if (log10Sum >= 0)
            throw new IllegalArgumentException("heterosygosity " + heterozygosity + " is too large of total ploidy " + totalPloidy);
        result[0] = MathUtils.log10OneMinusPow10(log10Sum);
        return result;
    }
}
