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

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;

/**
 * Instantiates Exact AF calculators given the required ploidy specs.
 *
 * <p>Unless you know better implementations of this class might return the same instance several times
 * and so the client code might need to make sure that there are no collisions or race conditions.</p>
 *
 * Author: wbc
 */
public abstract class AFCalculatorProvider {

    /**
     * Returns a AF calculator capable to handle a particular variant-context.
     * @param variantContext the target context build.
     * @param defaultPloidy the assumed ploidy in case that there is no a GT call present to determine it.
     * @return never {@code null}
     */
    public AFCalculator getInstance(final VariantContext variantContext, final int defaultPloidy, final int maximumAltAlleles) {
        if (variantContext == null)
            throw new IllegalArgumentException("variant context cannot be null");

        final int sampleCount = variantContext.getNSamples();
        if  (sampleCount == 0)
            return getInstance(defaultPloidy,maximumAltAlleles);

        final GenotypesContext genotypes = variantContext.getGenotypes();

        final Genotype firstGenotype = genotypes.get(0);
        int ploidy = firstGenotype.getPloidy();
        if (ploidy <= 0) ploidy = defaultPloidy;
        for (int i = 1 ; i < sampleCount; i++) {
            final Genotype genotype = genotypes.get(i);
            final int declaredPloidy = genotype.getPloidy();
            final int actualPloidy = declaredPloidy <= 0 ? defaultPloidy : declaredPloidy;
            if (actualPloidy != ploidy) {
                ploidy = AFCalculatorImplementation.UNBOUND_PLOIDY;
                break;
            }
        }
        return getInstance(ploidy,Math.min(variantContext.getNAlleles() - 1, maximumAltAlleles));
    }

    /**
     * Returns a AF calculator given the required homogeneous ploidy and allele count (including the reference).
     * @param ploidy the required ploidy.
     * @param maximumAltAlleles the allele count.
     * @return never {@code null}
     */
    public abstract AFCalculator getInstance(final int ploidy, final int maximumAltAlleles);

}