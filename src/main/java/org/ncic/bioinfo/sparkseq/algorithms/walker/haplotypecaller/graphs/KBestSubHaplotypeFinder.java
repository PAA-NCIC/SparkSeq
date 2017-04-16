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
package org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.graphs;

import org.ncic.bioinfo.sparkseq.algorithms.data.basic.Pair;

import java.util.Set;

/**
 * Common interface for K-Best sub-haplotype finders.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
interface KBestSubHaplotypeFinder {

    /**
     * Return an unique id for this sub-haplotype finder to be used when outputting diagrams.
     *
     * @return never {@code null}.
     */
    public String id();

    /**
     * Returns a label with human readable representation of this finder.
     *
     * <p>This is used when generating a diagram to illustrate the search space and costs</p>
     *
     * @return never {@code null}.
     */
    public String label();

    /**
     * Returns the set of subfinder from this finder together with a label for the connection with the current finder.
     *
     * <p>The label is used when generating a diagram to illustrate the search space and costs</p>
     */
    public Set<Pair<? extends KBestSubHaplotypeFinder,String>> subFinderLabels();

    /**
     * Returns the total number of possible sub-haplotypes.
     * @return 0 or greater.
     */
    public abstract int getCount();

    /**
     * Return the k-best sub-haplotype solution.
     *
     *
     * @param k the requested solution rank.
     * @throws IllegalArgumentException if {@code k} is outside bounds [0 .. {@link #getCount()}).
     *
     * @return never {@code null}.
     */
    public KBestHaplotype getKBest(int k);

    /**
     * Checks whether the top vertex for this finder is a reference haplotype vertex.
     *
     * @return {@code true} iff the top vertex for this finder is a reference vertex.
     */
    public boolean isReference();

    /**
     * Calculate the score of a sequence of bases.
     *
     * @param bases array containing the query base sequence.
     * @param offset first position of the query base sequence in {@code bases} .
     * @param length length of the query base sequence.
     * @return {@link Double#NaN} if there is no score for this sequence, otherwise a valid score value.
     */
    public double score(byte[] bases, int offset, int length);
}
