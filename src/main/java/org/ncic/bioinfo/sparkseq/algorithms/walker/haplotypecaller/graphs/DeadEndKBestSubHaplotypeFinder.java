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

import java.util.Collections;
import java.util.Set;

/**
 * Represents a trivial k-best sub haplotype finder with no solutions.
 *
 * <p>To be used at vertices that do not have any valid path to the requested sink vertices</p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
final class DeadEndKBestSubHaplotypeFinder implements KBestSubHaplotypeFinder {

    /**
     * Sole instance of this class.
     */
    public static DeadEndKBestSubHaplotypeFinder INSTANCE = new DeadEndKBestSubHaplotypeFinder();

    /**
     * Prevents instantiation of more than one instance; please use {@link #INSTANCE}.
     */
    protected DeadEndKBestSubHaplotypeFinder() {
    }

    @Override
    public String id() {
        return "<DEAD>";
    }

    @Override
    public String label() {
        return "&lt;DEAD&gt;";
    }

    @Override
    public Set<Pair<? extends KBestSubHaplotypeFinder, String>> subFinderLabels() {
        return Collections.emptySet();
    }

    @Override
    public int getCount() {
        return 0;
    }

    @Override
    public KBestHaplotype getKBest(int k) {
        if (k < 0)
            throw new IllegalArgumentException("k cannot be negative");
        else
            throw new IllegalArgumentException("k cannot be equal or greater to the haplotype count");
    }

    @Override
    public boolean isReference() {
        return false;
    }

    @Override
    public double score(final byte[] bases, final int offset, final int length) {
        if (bases == null) throw new IllegalArgumentException("bases cannot be null");
        if (offset < 0) throw new IllegalArgumentException("the offset cannot be negative");
        if (length < 0) throw new IllegalArgumentException("the length cannot be negative");
        if (offset + length > bases.length) throw new IllegalArgumentException("the offset and length go beyond the array size");
        return Double.NaN;
    }
}
