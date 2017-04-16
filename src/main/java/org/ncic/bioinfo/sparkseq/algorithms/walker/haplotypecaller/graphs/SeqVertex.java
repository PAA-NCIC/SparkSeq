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
import java.util.concurrent.atomic.AtomicInteger;

/**
 * A graph vertex containing a sequence of bases and a unique ID that
 * allows multiple distinct nodes in the graph to have the same sequence.
 *
 * This is essential when thinking about representing the actual sequence of a haplotype
 * in a graph.  There can be many parts of the sequence that have the same sequence, but
 * are distinct elements in the graph because they have a different position in the graph.  For example:
 *
 * A -> C -> G -> A -> T
 *
 * The two As are not the same, because they occur with different connections.  In a kmer graph equals()
 * is based on the sequence itself, as each distinct kmer can only be represented once.  But the transformation
 * of the kmer graph into a graph of base sequences, without their kmer prefixes, means that nodes that
 * where once unique including their prefix can become equal after shedding the prefix.  So we need to
 * use some mechanism -- here a unique ID per node -- to separate nodes that have the same sequence
 * but are distinct elements of the graph.
 *
 * @author: depristo
 * @since 03/2013
 */
public final class SeqVertex extends BaseVertex {
    // Note that using an AtomicInteger is critical to allow multi-threaded HaplotypeCaller
    private static final AtomicInteger idCounter = new AtomicInteger(0);
    private int id = idCounter.getAndIncrement();

    /**
     * Create a new SeqVertex with sequence and the next available id
     * @param sequence our base sequence
     */
    public SeqVertex(final byte[] sequence) {
        super(sequence);
    }

    /**
     * Create a new SeqVertex having bases of sequence.getBytes()
     * @param sequence the string representation of our bases
     */
    public SeqVertex(final String sequence) {
        super(sequence);
    }

    /**
     * Create a copy of toCopy
     * @param toCopy a SeqVertex to copy into this newly allocated one
     */
    public SeqVertex(final SeqVertex toCopy) {
        super(toCopy.sequence);
        this.id = toCopy.id;
    }

    /**
     * Get the unique ID for this SeqVertex
     * @return a positive integer >= 0
     */
    public int getId() {
        return id;
    }

    @Override
    public String toString() {
        return "SeqVertex_id_" + id + "_seq_" + getSequenceString();
    }

    /**
     * Two SeqVertex are equal only if their ids are equal
     * @param o
     * @return
     */
    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        SeqVertex seqVertex = (SeqVertex) o;
        if (id != seqVertex.id) return false;

        // note that we don't test for super equality here because the ids are unique
        //if (!super.equals(o)) return false;

        return true;
    }

    @Override
    public int hashCode() {
        return id;
    }

    /**
     * Return a new SeqVertex derived from this one but not including the suffix bases
     *
     * @param suffix the suffix bases to remove from this vertex
     * @return a newly allocated SeqVertex with appropriate prefix, or null if suffix removes all bases from this node
     */
    public SeqVertex withoutSuffix(final byte[] suffix) {
        final int prefixSize = sequence.length - suffix.length;
        return prefixSize > 0 ? new SeqVertex(Arrays.copyOf(sequence, prefixSize)) : null;
    }

    /**
     * Return a new SeqVertex derived from this one but not including prefix or suffix bases
     *
     * @param prefix the previx bases to remove
     * @param suffix the suffix bases to remove from this vertex
     * @return a newly allocated SeqVertex
     */
    public SeqVertex withoutPrefixAndSuffix(final byte[] prefix, final byte[] suffix) {
        final int start = prefix.length;
        final int length = sequence.length - suffix.length - prefix.length;
        final int stop = start + length;
        return length > 0 ? new SeqVertex(Arrays.copyOfRange(sequence, start, stop)) : null;
    }
}
