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
package org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.readthreading;

import org.ncic.bioinfo.sparkseq.algorithms.utils.Utils;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.graphs.DeBruijnVertex;

import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * A DeBruijnVertex that supports multiple copies of the same kmer
 *
 * Author: wbc
 */
public final class MultiDeBruijnVertex extends DeBruijnVertex {
    private final static boolean KEEP_TRACK_OF_READS = false;

    // Note that using an AtomicInteger is critical to allow multi-threaded HaplotypeCaller
    private static final AtomicInteger idCounter = new AtomicInteger(0);
    private int id = idCounter.getAndIncrement();

    private final List<String> reads = new LinkedList<String>();

    /**
     * Create a new MultiDeBruijnVertex with kmer sequence
     * @param sequence the kmer sequence
     */
    MultiDeBruijnVertex(byte[] sequence) {
        super(sequence);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        MultiDeBruijnVertex that = (MultiDeBruijnVertex) o;

        return id == that.id;
    }

    @Override
    public String toString() {
        return "MultiDeBruijnVertex_id_" + id + "_seq_" + getSequenceString();
    }

    /**
     * Add name information to this vertex for debugging
     *
     * This information will be captured as a list of strings, and displayed in DOT if this
     * graph is written out to disk
     *
     * This functionality is only enabled when KEEP_TRACK_OF_READS is true
     *
     * @param name a non-null string
     */
    protected void addRead(final String name) {
        if ( name == null ) throw new IllegalArgumentException("name cannot be null");
        if ( KEEP_TRACK_OF_READS ) reads.add(name);
    }

    @Override
    public int hashCode() { return id; }

    @Override
    public String additionalInfo() {
        return super.additionalInfo() + (KEEP_TRACK_OF_READS ? (! reads.contains("ref") ? "__" + Utils.join(",", reads) : "") : "");
    }

     int getId() {
        return id;
    }
}
