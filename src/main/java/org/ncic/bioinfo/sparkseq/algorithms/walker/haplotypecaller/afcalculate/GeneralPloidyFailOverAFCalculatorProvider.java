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

import org.apache.log4j.Logger;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.argcollection.GenotypeCalculationArgumentCollection;

/**
 * Author: wbc
 */
public class GeneralPloidyFailOverAFCalculatorProvider extends AFCalculatorProvider {

    private final AFCalculator preferred;
    private final AFCalculatorImplementation preferredImplementation;
    private final AFCalculator failOver;

    /**
     * Creates a new AF calculator provider given the genotyping arguments and logger reference.
     * @param genotypeArgs genotyping parameter collection.
     * @param logger where the AF calculator logging messages go. If {@code null}, logging message will not be emitted.
     *
     * @throws IllegalStateException if {@code genotypeArgs} is {@code null}.
     */
    public GeneralPloidyFailOverAFCalculatorProvider(final GenotypeCalculationArgumentCollection genotypeArgs, final Logger logger) {
        if (genotypeArgs == null)
            throw new IllegalArgumentException("genotype arguments object cannot be null");
        preferredImplementation = AFCalculatorImplementation.bestValue(genotypeArgs.samplePloidy,genotypeArgs.MAX_ALTERNATE_ALLELES, null);
        preferred = preferredImplementation.newInstance();
        preferred.setLogger(logger);
        failOver = AFCalculatorImplementation.EXACT_GENERAL_PLOIDY.newInstance();
        failOver.setLogger(logger);
    }

    /**
     * {@inheritDoc}
     * @param ploidy {@inheritDoc}
     * @param maximumAlternativeAlleles {@inheritDoc}
     * @return {@inheritDoc}
     */
    @Override
    public AFCalculator getInstance(final int ploidy, final int maximumAlternativeAlleles) {
        return preferredImplementation.usableForParams(ploidy,maximumAlternativeAlleles) ? preferred : failOver;
    }

    /**
     * Creates a AF calculator provider that complies to the contract of this class and that make sures that is
     * safe to use in multi-thread walker runs.
     *
     * <p>
     *    Whether this is part of a multi-thread run is determined by looking into the engine configuration in {@code toolkit}.
     * </p>
     *
     * @param genotypeArgs the genotyping arguments.
     * @param logger where to output logging messages for instantiated AF calculators.
     * @return never {@code null}.
     */
    public static AFCalculatorProvider createThreadSafeProvider(final GenotypeCalculationArgumentCollection genotypeArgs,
                                                                final Logger logger) {
        if (genotypeArgs == null)
            throw new IllegalArgumentException("genotype arguments object cannot be null");
        return new GeneralPloidyFailOverAFCalculatorProvider(genotypeArgs,logger);
    }
}
