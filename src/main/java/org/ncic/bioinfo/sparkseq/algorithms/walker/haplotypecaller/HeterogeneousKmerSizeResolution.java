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

/**
 * Author: wbc
 */
public enum HeterogeneousKmerSizeResolution {

    /**
     * Combine haplotypes using a haplotype graph with the largest kmerSize amongst the ones that generated some haplotype.
     */
    COMBO_MAX,

    /**
     * Combine haplotypes using a haplotype graph with the largest kmerSize amongst the ones that generated some haplotype.
     */
    COMBO_MIN,

    /**
     * Take just the haplotypes from largest kmersize that generated any.
     */
    MAX_ONLY,

    /**
     * Take just the haplotypes from the smallest kmerSize that generated any.
     */
    @SuppressWarnings("unused")
    MIN_ONLY;

    /**
     * Indicates whether we should use the maximum kmerSize for the haplotypeGraph or not.
     *
     * @return true if we need to use the maximum, false otherwise.
     */
    public boolean useMaximum() {
        switch (this) {
            case COMBO_MAX: return true;
            case MAX_ONLY: return true;
            default: return false;
        }
    }

    /**
     * Indicates whether we should use the minimum kmerSize for the haplotypeGraph or not.
     *
     * @return true if we need to use the minimum, false otherwise.
     */
    @SuppressWarnings("unused")
    public boolean useMinimum() {
        return ! useMaximum();
    }

    /**
     * Tell whether this policy combines kmer-sizes or not.
     * @return true iff it does.
     */
    public boolean combinesKmerSizes() {
        switch (this) {
            case COMBO_MAX: return true;
            case COMBO_MIN: return true;
            default: return false;
        }

    }
}
