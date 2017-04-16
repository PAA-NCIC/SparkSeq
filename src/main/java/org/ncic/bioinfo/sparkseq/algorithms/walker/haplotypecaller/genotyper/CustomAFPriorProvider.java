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

import java.util.List;

/**
 * Author: wbc
 */
public class CustomAFPriorProvider extends AFPriorProvider {

    private final double[] priors;

    /**
     *
     * @param priorValues must exactly the number of genomes in the samples (the total ploidy).
     */
    public CustomAFPriorProvider(final List<Double> priorValues) {
        if (priorValues == null)
            throw new IllegalArgumentException("the input prior values cannot be null");
        priors = new double[priorValues.size() + 1];

        int i = 1;
        double sum = 0;
        for (double value : priorValues) {
            if (value <= 0 || value >= 1)
                throw new IllegalArgumentException("the AF prior value "+ value + " is out of the valid interval (0,1)");
            if (Double.isNaN(value))
                throw new IllegalArgumentException("NaN is not a valid prior AF value");
            priors[i++] = Math.log10(value);
            sum += value;
        }
        if (sum >= 1)
            throw new IllegalArgumentException("the AF prior value sum must be less than 1: " + sum);
        priors[0] = Math.log10(1 - sum);
    }

    @Override
    protected double[] buildPriors(final int totalPloidy) {
        if (totalPloidy != priors.length - 1)
            throw new IllegalStateException("requesting an invalid prior total ploidy " + totalPloidy + " != " + (priors.length - 1));
        return priors;
    }
}
