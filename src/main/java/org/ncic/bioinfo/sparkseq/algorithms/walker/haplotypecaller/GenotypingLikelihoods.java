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
package org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;

import java.util.List;

/**
 * Author: wbc
 */
public class GenotypingLikelihoods<A extends Allele> implements SampleList, AlleleList<A> {

    private final GenotypeLikelihoods[] likelihoods;

    private final PloidyModel ploidyModel;

    private final AlleleList<A> alleles;

    /**
     * Creates a new genotyping-likelihoods collection given the genotype alleles, the sample ploidy model and the
     *   likelihoods themselves.
     * <p>
     * Notice that this constructor does not check whether the likelihood array lengths corresponds to the sample plodies and
     * number of alleles.
     * </p>
     *
     * @param alleles the genotyping alleles.
     * @param ploidyModel the ploidy model.
     * @param likelihoods the actual genotype likelihoods, one element per sample.
     *
     * @throws IllegalArgumentException if any argument is {@code null}, or the number of samples in {@code ploidyModel}
     *  does not correspond with the number of likelihoods arrays in {@code likelihoods}
     */
    public GenotypingLikelihoods(final AlleleList<A> alleles, final PloidyModel ploidyModel,
                          final List<GenotypeLikelihoods> likelihoods) {
        if (alleles == null)
            throw new IllegalArgumentException("allele list cannot be null");
        if (ploidyModel == null)
            throw new IllegalArgumentException("the ploidy model cannot be null");
        if (likelihoods == null)
            throw new IllegalArgumentException("the likelihood collection cannot be null");
        if (ploidyModel.sampleCount() !=  likelihoods.size())
            throw new IllegalArgumentException("there must be exactly one likelihood set for each sample");

        this.likelihoods = likelihoods.toArray(new GenotypeLikelihoods[likelihoods.size()]);
        for (final GenotypeLikelihoods likelihood : this.likelihoods)
            if (likelihood == null)
                throw new IllegalArgumentException("no genotype likelihood is allowed to be null");

        this.alleles = alleles;
        this.ploidyModel = ploidyModel;
    }

    @Override
    public int sampleCount() {
        return ploidyModel.sampleCount();
    }

    @Override
    public int sampleIndex(final String sample) {
        return ploidyModel.sampleIndex(sample);
    }

    @Override
    public String sampleAt(final int sampleIndex) {
        return ploidyModel.sampleAt(sampleIndex);
    }

    /**
     * Returns the ploidy of the sample given its index in the collection.
     *
     * @param sampleIndex the query sample index.
     *
     * @throws IllegalArgumentException if {@code sampleIndex} is not a valid index for this collection:
     *   [0,{@link #sampleCount()).
     *
     * @return 0 or greater.
     */
    public int samplePloidy(final int sampleIndex) {
        return ploidyModel.samplePloidy(sampleIndex);
    }

    /**
     * Returns the genotype-likelihoods of the sample given its index in the collection.
     *
     * @param sampleIndex the query sample index.
     *
     * @throws IllegalArgumentException if {@code sampleIndex} is not a valid index for this collection:
     *   [0,{@link #sampleCount()).
     *
     * @return never {@code null}.
     */
    public GenotypeLikelihoods sampleLikelihoods(final int sampleIndex) {
        return likelihoods[sampleIndex];
    }

    @Override
    public int alleleCount() {
        return alleles.alleleCount();
    }

    @Override
    public int alleleIndex(final A allele) {
        return alleles.alleleIndex(allele);
    }

    @Override
    public A alleleAt(final int index) {
        return alleles.alleleAt(index);
    }
}
