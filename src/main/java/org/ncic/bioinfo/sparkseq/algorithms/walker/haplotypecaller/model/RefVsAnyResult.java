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
package org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.model;

import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.HomoSapiensConstants;

/**
 * Author: wbc
 */
public final class RefVsAnyResult {
    /**
     * The genotype likelihoods for ref/ref ref/non-ref non-ref/non-ref
     */
    public final double[] genotypeLikelihoods;

    /**
     * AD field value for ref / non-ref
     */
    final int[] AD_Ref_Any = new int[2];

    /**
     * @return Get the DP (sum of AD values)
     */
    protected int getDP() { return AD_Ref_Any[0] + AD_Ref_Any[1]; }

    /**
     * Cap the het and hom var likelihood values by the hom ref likelihood.
     */
    protected void capByHomRefLikelihood() {
        final int likelihoodCount = genotypeLikelihoods.length;
        for (int i = 1; i < likelihoodCount; i++)
            genotypeLikelihoods[i] = Math.min(genotypeLikelihoods[0],genotypeLikelihoods[i]);
    }

    /**
     * Creates a new ref-vs-alt result assuming 3 as the number of genotype likelihoods (human ploidy.
     */
    @Deprecated
    public RefVsAnyResult() {
        genotypeLikelihoods =
                new double[(HomoSapiensConstants.DEFAULT_PLOIDY * (HomoSapiensConstants.DEFAULT_PLOIDY + 1)) >> 1];
    }

    /**
     * Creates a new ref-vs-alt result indicating the genotype likelihood vector capacity.
     * @param likelihoodCapacity the required capacity of the likelihood array, should match the possible number of
     *                           genotypes given the number of alleles (always 2), ploidy (arbitrary) less the genotyping
     *                           model non-sense genotype count if applies.
     * @throws IllegalArgumentException if {@code likelihoodCapacity} is negative.
     */
    public RefVsAnyResult(final int likelihoodCapacity) {
        if (likelihoodCapacity < 0)
            throw new IllegalArgumentException("likelihood capacity is negative");
        genotypeLikelihoods = new double[likelihoodCapacity];
    }
}
