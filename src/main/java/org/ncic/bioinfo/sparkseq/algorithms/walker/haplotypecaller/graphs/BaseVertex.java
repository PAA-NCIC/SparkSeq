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

import java.util.Arrays;

/**
 * A graph vertex that holds some sequence information
 *
 * @author: depristo
 * @since 03/2013
 */
public class BaseVertex {
    /** placeholder to store additional information for debugging purposes */
    String additionalInfo = "";
    final byte[] sequence;
    private final static int UNASSIGNED_HASHCODE = -1;
    int cachedHashCode = UNASSIGNED_HASHCODE;

    /**
     * Create a new sequence vertex with sequence
     *
     * This code doesn't copy sequence for efficiency reasons, so sequence should absolutely not be modified
     * in any way after passing this sequence to the BaseVertex
     *
     * @param sequence a non-null, non-empty sequence of bases contained in this vertex
     */
    public BaseVertex(final byte[] sequence) {
        if ( sequence == null ) throw new IllegalArgumentException("Sequence cannot be null");
        this.sequence = sequence;
    }

    /**
     * Does this vertex have an empty sequence?
     *
     * That is, is it a dummy node that's only present for structural reasons but doesn't actually
     * contribute to the sequence of the graph?
     *
     * @return true if sequence is empty, false otherwise
     */
    public boolean isEmpty() {
        return length() == 0;
    }

    /**
     * Get the length of this sequence
     * @return a positive integer >= 1
     */
    public int length() {
        return sequence.length;
    }

    /**
     * For testing purposes only -- low performance
     * @param sequence the sequence as a string
     */
    protected BaseVertex(final String sequence) {
        this(sequence.getBytes());
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        BaseVertex that = (BaseVertex) o;

        if (!Arrays.equals(sequence, that.sequence)) return false;

        return true;
    }

    /**
     * Are b and this equal according to their base sequences?
     *
     * @param b the vertex to compare ourselves to
     * @return true if b and this have the same sequence, regardless of other attributes that might differentiate them
     */
    public boolean seqEquals(final BaseVertex b) {
        return Arrays.equals(this.getSequence(), b.getSequence());
    }

    /**
     * necessary to override here so that graph.containsVertex() works the same way as vertex.equals() as one might expect
     * @return
     */
    @Override
    public int hashCode() {
        if ( cachedHashCode == UNASSIGNED_HASHCODE ) {
            cachedHashCode = Arrays.hashCode(sequence);
        }
        return cachedHashCode;
    }

    @Override
    public String toString() {
        return getSequenceString();
    }

    /**
     * Get the sequence of bases contained in this vertex
     *
     * Do not modify these bytes in any way!
     *
     * @return a non-null pointer to the bases contained in this vertex
     */
    public byte[] getSequence() {
        return sequence;
    }

    /**
     * Get a string representation of the bases in this vertex
     * @return a non-null String
     */
    public String getSequenceString() {
        return new String(sequence);
    }

    /**
     * Get the sequence unique to this vertex
     *
     * This function may not return the entire sequence stored in the vertex, as kmer graphs
     * really only provide 1 base of additional sequence (the last base of the kmer).
     *
     * The base implementation simply returns the sequence.
     *
     * @param source is this vertex a source vertex (i.e., no in nodes) in the graph
     * @return a byte[] of the sequence added by this vertex to the overall sequence
     */
    public byte[] getAdditionalSequence(final boolean source) {
        return getSequence();
    }

    /**
     * Set additional debugging information for this vertex
     * @param info the new info value.
     */
    public void setAdditionalInfo(final String info) {
        if ( info == null ) throw new IllegalArgumentException("info cannot be null");
        additionalInfo = info;
    }

    /**
     * @return the additional information for display about this vertex
     */
    public String additionalInfo() { return additionalInfo; }

    /**
     * Checks whether the vertex sequence is ambiguous or not.
     *
     * <p>
     *     Ambiguity may come about as a result of either:
     *     <ul>
     *        <li>by construction as the generating sequence (read or haplotype) had ambiguous bases</li>
     *        <li>or because this vertex is the result of merging two or more vertices with some variation upstream
     *        no more than kmerSize bases away </li>
     *     </ul>
     * </p>
     *
     * @return {@code true} iff so.
     */
    public boolean hasAmbiguousSequence() {
        for (final byte base : sequence)
            switch (Character.toUpperCase(base)) {
                case 'A' :
                case 'T' :
                case 'G' :
                case 'C' :
                    continue;
                default :
                    return true;
            }
        return false;
    }
}
