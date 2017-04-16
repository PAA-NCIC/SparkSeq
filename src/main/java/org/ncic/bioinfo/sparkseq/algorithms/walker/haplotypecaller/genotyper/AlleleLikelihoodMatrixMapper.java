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

import htsjdk.variant.variantcontext.Allele;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.GATKSAMRecord;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.AlleleListPermutation;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.readlikelihood.ReadLikelihoods;

import java.util.List;

/**
 * Author: wbc
 */
public abstract class AlleleLikelihoodMatrixMapper<A extends Allele> {

    public abstract ReadLikelihoods.Matrix<A> map(final ReadLikelihoods.Matrix<A> original);

    /**
     * Instantiates a new mapper given an allele-list permutation.
     * @param permutation the requested permutation.
     * @param <A> the allele type.
     *
     * @throws IllegalArgumentException if {@code permutation} is {@code null}.
     *
     * @return never {@code null}.
     */
    public static <A extends Allele> AlleleLikelihoodMatrixMapper<A> newInstance(final AlleleListPermutation<A> permutation) {
        if (permutation == null)
            throw new IllegalArgumentException("the permutation must not be null");
        if (permutation.isNonPermuted())
            return asIs();
        else
            return general(permutation);
    }

    /**
     * Returns trivial mapper that just maps to the original matrix without changes.
     *
     * @param <A> the allele type.
     * @return never {@code null}.
     */
    @SuppressWarnings("unchecked")
    private static <A extends Allele> AlleleLikelihoodMatrixMapper<A> asIs() {
        return AS_IS;
    }

    @SuppressWarnings("unchecked")
    private static final AlleleLikelihoodMatrixMapper AS_IS = new AlleleLikelihoodMatrixMapper<Allele>() {
        @Override
        public ReadLikelihoods.Matrix<Allele> map(final ReadLikelihoods.Matrix original) {
            return original;
        }
    };

    /**
     * Constructs a new mapper instance that work with general permutation without making any assumption.
     * @param permutation the permutation to apply to requested matrices wrappers.
     * @param <A> allele type.
     * @return never {@code null}.
     */
    private static <A extends Allele> AlleleLikelihoodMatrixMapper<A> general(final AlleleListPermutation<A> permutation) {
        return new AlleleLikelihoodMatrixMapper<A>() {
            @Override
            public ReadLikelihoods.Matrix<A> map(final ReadLikelihoods.Matrix<A> original) {
                return new ReadLikelihoods.Matrix<A>() {

                    @Override
                    public List<GATKSAMRecord> reads() {
                        return original.reads();
                    }

                    @Override
                    public List<A> alleles() {
                        return permutation.toList();
                    }

                    @Override
                    public void set(final int alleleIndex, final int readIndex, final double value) {
                        original.set(permutation.fromIndex(alleleIndex),readIndex,value);
                    }

                    @Override
                    public double get(final int alleleIndex, final int readIndex) {
                        return original.get(permutation.fromIndex(alleleIndex),readIndex);
                    }

                    @Override
                    public int alleleIndex(final A allele) {
                        return permutation.alleleIndex(allele);
                    }

                    @Override
                    public int readIndex(GATKSAMRecord read) {
                        return original.readIndex(read);
                    }

                    @Override
                    public int alleleCount() {
                        return permutation.toSize();
                    }

                    @Override
                    public int readCount() {
                        return original.readCount();
                    }

                    @Override
                    public A alleleAt(final int alleleIndex) {
                        return original.alleleAt(permutation.fromIndex(alleleIndex));
                    }

                    @Override
                    public GATKSAMRecord readAt(final int readIndex) {
                        return original.readAt(readIndex);
                    }

                    @Override
                    public void copyAlleleLikelihoods(final int alleleIndex, final double[] dest, final int offset) {
                        original.copyAlleleLikelihoods(permutation.fromIndex(alleleIndex),dest,offset);
                    }
                };
            }
        };
    }
}
