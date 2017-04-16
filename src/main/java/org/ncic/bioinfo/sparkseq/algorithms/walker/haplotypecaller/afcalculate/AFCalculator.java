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
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.log4j.Logger;

import java.util.List;

/**
 * Author: wbc
 */
public abstract class AFCalculator implements Cloneable {
    private final static Logger defaultLogger = Logger.getLogger(AFCalculator.class);


    protected Logger logger = defaultLogger;

    private StateTracker stateTracker;

    /**
     * Create a new AFCalc object capable of calculating the prob. that alleles are
     * segregating among many samples.
     *
     * <p>
     *    Restrictions in ploidy and number of alternative alleles that a instance can handle will be determined
     *    by its implementation class {@link AFCalculatorImplementation}
     * </p>
     */
    protected AFCalculator() {
    }

    /**
     * Use this logger instead of the default logger
     *
     * @param logger
     */
    public void setLogger(Logger logger) {
        this.logger = logger;
    }

    /**
     * Compute the probability of the alleles segregating given the genotype likelihoods of the samples in vc
     *
     * @param vc the VariantContext holding the alleles and sample information.  The VariantContext
     *           must have at least 1 alternative allele
     * @param log10AlleleFrequencyPriors a prior vector nSamples x 2 in length indicating the Pr(AF = i)
     * @return result (for programming convenience)
     */
    public AFCalculationResult getLog10PNonRef(final VariantContext vc, final int defaultPloidy, final int maximumAlternativeAlleles, final double[] log10AlleleFrequencyPriors) {
        if ( vc == null ) throw new IllegalArgumentException("VariantContext cannot be null");
        if ( vc.getNAlleles() == 1 ) throw new IllegalArgumentException("VariantContext has only a single reference allele, but getLog10PNonRef requires at least one at all " + vc);
        if ( log10AlleleFrequencyPriors == null ) throw new IllegalArgumentException("priors vector cannot be null");

        // reset the result, so we can store our new result there
        final StateTracker stateTracker = getStateTracker(true,maximumAlternativeAlleles);

        final VariantContext vcWorking = reduceScope(vc,defaultPloidy, maximumAlternativeAlleles);

        final AFCalculationResult result = computeLog10PNonRef(vcWorking, defaultPloidy, log10AlleleFrequencyPriors, stateTracker);
        return result;
    }

    /**
     * Convert the final state of the state tracker into our result as an AFCalculationResult
     *
     * Assumes that stateTracker has been updated accordingly
     *
     * @param vcWorking the VariantContext we actually used as input to the calc model (after reduction)
     * @param log10AlleleFrequencyPriors the priors by AC vector
     * @return a AFCalculationResult describing the result of this calculation
     */
    protected AFCalculationResult getResultFromFinalState(final VariantContext vcWorking, final double[] log10AlleleFrequencyPriors, final StateTracker stateTracker) {
        stateTracker.setAllelesUsedInGenotyping(vcWorking.getAlleles());
        return stateTracker.toAFCalculationResult(log10AlleleFrequencyPriors);
    }

    // ---------------------------------------------------------------------------
    //
    // Abstract methods that should be implemented by concrete implementations
    // to actually calculate the AF
    //
    // ---------------------------------------------------------------------------

    /**
     * Look at VC and perhaps return a new one of reduced complexity, if that's necessary
     *
     * Used before the call to computeLog10PNonRef to simply the calculation job at hand,
     * if vc exceeds bounds.  For example, if VC has 100 alt alleles this function
     * may decide to only genotype the best 2 of them.
     *
     * @param vc the initial VC provided by the caller to this AFcalculation
     * @return a potentially simpler VC that's more tractable to genotype
     */
    protected abstract VariantContext reduceScope(final VariantContext vc, final int defaultPloidy, final int maximumAlternativeAlleles);

    /**
     * Actually carry out the log10PNonRef calculation on vc, storing results in results
     *
     * @param vc                                variant context with alleles and genotype likelihoods,
     *                                          must have at least one alt allele
     * @param log10AlleleFrequencyPriors        priors
     * @return a AFCalcResult object describing the results of this calculation
     */
    protected abstract AFCalculationResult computeLog10PNonRef(final VariantContext vc, final int defaultPloidy,
                                                               final double[] log10AlleleFrequencyPriors, final StateTracker stateTracker);

    /**
     * Subset VC to the just allelesToUse, updating genotype likelihoods
     *
     * Must be overridden by concrete subclasses
     *
     * @param vc                                variant context with alleles and genotype likelihoods
     * @param defaultPloidy                     default ploidy to assume in case {@code vc} does not indicate it for a sample.
     * @param allelesToUse                      alleles to subset
     * @param assignGenotypes
     * @return GenotypesContext object
     */
    public abstract GenotypesContext subsetAlleles(final VariantContext vc,
                                                   final int defaultPloidy,
                                                   final List<Allele> allelesToUse,
                                                   final boolean assignGenotypes);

    // ---------------------------------------------------------------------------
    //
    // accessors
    //
    // ---------------------------------------------------------------------------


    /**
     * Retrieves the state tracker.
     *
     * <p>
     *     The tracker will be reset if so requested or if it needs to be resized due to an increase in the
     *     maximum number of alleles is must be able to handle.
     * </p>
     *
     * @param reset make sure the tracker is reset.
     * @param maximumAlternativeAlleleCount the maximum alternative allele count it must be able to handle. Has no effect if
     *                                     the current tracker is able to handle that number.
     *
     * @return never {@code null}
     */
    protected StateTracker getStateTracker(final boolean reset, final int maximumAlternativeAlleleCount) {
        if (stateTracker == null)
            stateTracker = new StateTracker(maximumAlternativeAlleleCount);
        else if (reset)
            stateTracker.reset(maximumAlternativeAlleleCount);
        else
            stateTracker.ensureMaximumAlleleCapacity(maximumAlternativeAlleleCount);
        return stateTracker;
    }

    /**
     * Used by testing code.
     *
     * Please don't use this method in production.
     *
     * @deprecated
     */
    @Deprecated
    protected int getAltAlleleCountOfMAP(final int allele) {
        return getStateTracker(false,allele + 1).getAlleleCountsOfMAP()[allele];
    }


}
