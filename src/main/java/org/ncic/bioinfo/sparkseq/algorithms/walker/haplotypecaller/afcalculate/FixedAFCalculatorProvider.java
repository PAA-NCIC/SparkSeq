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

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.log4j.Logger;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.argcollection.GenotypeCalculationArgumentCollection;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.argcollection.StandardCallerArgumentCollection;

/**
 * Author: wbc
 */
public class FixedAFCalculatorProvider extends AFCalculatorProvider {

    private final AFCalculator singleton;

    private final boolean verifyRequests;

    private final int maximumAltAlleleCount;

    private final int ploidy;

    /**
     * Constructs a fixed AF Calculator provider.
     *
     * @param configuration  the called configuration. This is the source of the fixed ploidy and maximum number of
     *                       supported alleles.
     * @param logger         the logger to be used by the AF calculator instance.
     * @param verifyRequests whether this provider will verify that each request for the AF calculator meets the
     *                       initial parameter values (ploidy, sample-count and maximum number of alleles.
     * @throws NullPointerException if {@code configuration} is {@code null}, or it contains invalid values for
     *                              sample ploidy and maximum number of alternative alleles, or {@code sampleCount} is less than 0.
     */
    public FixedAFCalculatorProvider(final StandardCallerArgumentCollection configuration,
                                     final Logger logger, final boolean verifyRequests) {
        this(configuration.requestedAlleleFrequencyCalculationModel, configuration.genotypeArgs, logger, verifyRequests);

    }


    /**
     * Constructs a fixed AF Calculator provider.
     *
     * @param configuration  the called configuration. This is the source of the fixed ploidy and maximum number of
     *                       supported alleles.
     * @param logger         the logger to be used by the AF calculator instance.
     * @param verifyRequests whether this provider will verify that each request for the AF calculator meets the
     *                       initial parameter values (ploidy, sample-count and maximum number of alleles.
     * @throws IllegalArgumentException if {@code configuration} is {@code null}, or it contains invalid values for
     *                                  sample ploidy and maximum number of alternative alleles, or {@code sampleCount} is less than 0.
     */
    public FixedAFCalculatorProvider(final GenotypeCalculationArgumentCollection configuration,
                                     final Logger logger, final boolean verifyRequests) {
        this(null, configuration, logger, verifyRequests);
    }

    /**
     * Constructs a fixed AF Calculator provider.
     *
     * @param preferred      preferred implementation.
     * @param configuration  the called configuration. This is the source of the fixed ploidy and maximum number of
     *                       supported alleles.
     * @param logger         the logger to be used by the AF calculator instance.
     * @param verifyRequests whether this provider will verify that each request for the AF calculator meets the
     *                       initial parameter values (ploidy, sample-count and maximum number of alleles.
     * @throws IllegalArgumentException if {@code configuration} is {@code null}, or it contains invalid values for
     *                                  sample ploidy and maximum number of alternative alleles, or {@code sampleCount} is less than 0.
     */
    public FixedAFCalculatorProvider(final AFCalculatorImplementation preferred, final GenotypeCalculationArgumentCollection configuration,
                                     final Logger logger, final boolean verifyRequests) {

        if (configuration == null)
            throw new IllegalArgumentException("null configuration");
        if (configuration == null)
            throw new IllegalArgumentException("null configuration genotype arguments");
        if (configuration.samplePloidy < 1)
            throw new IllegalArgumentException("invalid sample ploidy " + configuration.samplePloidy);
        if (configuration.MAX_ALTERNATE_ALLELES < 0)
            throw new IllegalArgumentException("invalid maximum number of alleles " + (configuration.MAX_ALTERNATE_ALLELES + 1));

        ploidy = configuration.samplePloidy;
        maximumAltAlleleCount = configuration.MAX_ALTERNATE_ALLELES;
        singleton = AFCalculatorImplementation.bestValue(ploidy, maximumAltAlleleCount, preferred).newInstance();
        singleton.setLogger(logger);
        this.verifyRequests = verifyRequests;
    }

    @Override
    public AFCalculator getInstance(final VariantContext vc, final int defaultPloidy, final int maximumAlleleCount) {
        if (verifyRequests)
            // supers implementation will call eventually one of the other methods, so no need to verify anything here.
            return super.getInstance(vc, defaultPloidy, maximumAlleleCount);
        return singleton;
    }

    @Override
    public AFCalculator getInstance(final int ploidy, final int maxAltAlleleCount) {
        if (verifyRequests) {
            if (this.ploidy != AFCalculatorImplementation.UNBOUND_PLOIDY && ploidy != this.ploidy)
                throw new IllegalStateException("non-supported ploidy");
            if (maximumAltAlleleCount != AFCalculatorImplementation.UNBOUND_ALTERNATIVE_ALLELE_COUNT && maxAltAlleleCount > maximumAltAlleleCount)
                throw new IllegalStateException("non-supported alleleCount");
        }
        return singleton;
    }

    /**
     * Creates a fixed AF calculator provider that is thread safe given the engine configuration.
     *
     * @param config the caller configuration.
     * @param logger reference to the logger for the AF calculator to dump messages to.
     * @return never {@code null}
     * @throws IllegalArgumentException if any of the input argument is {@code null} or contain invalid configuration
     *                                  like zero-samples, zero or negative ploidy or negative-zero maximum number of alleles.
     */
    public static AFCalculatorProvider createThreadSafeProvider(final StandardCallerArgumentCollection config,
                                                                final Logger logger) {

        return new FixedAFCalculatorProvider(config, logger, false);
    }

    /**
     * Creates a fixed AF calculator provider that is thread safe given the engine configuration.
     *
     * @param config         the caller configuration.
     * @param logger         reference to the logger for the AF calculator to dump messages to.
     * @param verifyRequests whether each request should be check for a compatible ploidy and max-alt-allele count.
     *                       A non-compliant request would result in an exception when requesting the AF calculator instance.
     * @return never {@code null}.
     * @throws IllegalArgumentException if any of the input argument is {@code null} or contain invalid configuration
     *                                  like zero-samples, zero or negative ploidy or negative-zero maximum number of alleles.
     */
    @SuppressWarnings("unused")
    public static AFCalculatorProvider createThreadSafeProvider(
            final StandardCallerArgumentCollection config,
            final Logger logger, final boolean verifyRequests) {
        return new FixedAFCalculatorProvider(config, logger, verifyRequests);
    }

}
