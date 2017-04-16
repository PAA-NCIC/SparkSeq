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

import org.ncic.bioinfo.sparkseq.algorithms.utils.Utils;
import org.ncic.bioinfo.sparkseq.algorithms.data.basic.Pair;

import java.util.Collections;
import java.util.Set;

/**
 * Trivial k-best sub-haplotype finder where the source and sink vertex are the same one.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
class EmptyPathHaplotypeFinderNode implements KBestSubHaplotypeFinder {

    /**
     * Caches the only solution returned by this finder.
     */
    private final KBestHaplotype singleHaplotypePath;

    /**
     * Constructs a new empty k-best haplotype finder.
     *
     * @param graph the search graph.
     * @param vertex the source and sink vertex of the only solution returned by this finder.
     */
    public EmptyPathHaplotypeFinderNode(final SeqGraph graph, final SeqVertex vertex) {
        singleHaplotypePath = new MyBestHaplotypePath(graph,vertex);
    }

    @Override
    public String id() {
        return "v" + singleHaplotypePath.head().getId();
    }

    @Override
    public String label() {
        return singleHaplotypePath.head().getSequenceString();
    }

    @Override
    public Set<Pair<? extends KBestSubHaplotypeFinder, String>> subFinderLabels() {
        return Collections.emptySet();
    }

    @Override
    public int getCount() {
        return 1;
    }

    @Override
    public KBestHaplotype getKBest(int k) {
        if (k < 0)
            throw new IllegalArgumentException("k cannot be negative");
        if (k > 0)
            throw new IllegalArgumentException("k cannot greater than the possible haplotype count");
        return singleHaplotypePath;
    }

    @Override
    public boolean isReference() {
        return singleHaplotypePath.isReference();
    }

    @Override
    public double score(final byte[] bases, final int offset, final int length) {
        if (bases == null) throw new IllegalArgumentException("bases cannot be null");
        if (offset < 0) throw new IllegalArgumentException("the offset cannot be negative");
        if (length < 0) throw new IllegalArgumentException("the length cannot be negative");
        if (offset + length > bases.length) throw new IllegalArgumentException("the offset and length go beyond the array size");
        final byte[] vertexBases = singleHaplotypePath.head().getSequence();
        if (length != vertexBases.length)
            return Double.NaN;
        else
            return Utils.equalRange(bases, offset, vertexBases, 0, length)? 0 : Double.NaN;
    }

    /**
     * Custom extension of {@link KBestHaplotype} that implements the single solution behaviour.
     */
    private class MyBestHaplotypePath extends KBestHaplotype {

        /**
         * The solution's only vertex.
         */
        private final SeqVertex vertex;

        /**
         * The search graph.
         */
        private final SeqGraph graph;

        /**
         * Whether the vertex is a reference vertex.
         *
         * <p>Initialize lazily.</p>
         */
        private Boolean isReference;

        /**
         * Constructs a new empty k-best haplotype solution.
         *
         * @param graph the search graph.
         * @param vertex the source and sink vertex of the only solution returned by the outer finder.
         */
        public MyBestHaplotypePath(final SeqGraph graph, final SeqVertex vertex) {
            this.vertex = vertex;
            this.graph = graph;
        }

        @Override
        public SeqGraph graph() {
            return graph;
        }

        @Override
        public double score() {
            return 0;
        }

        @Override
        public int rank() {
            return 0;
        }

        @Override
        protected SeqVertex head() {
            return vertex;
        }

        @Override
        protected KBestHaplotype tail() {
            return null;
        }

        @Override
        public boolean isReference() {
            return (isReference != null) ? isReference: (isReference = graph.isReferenceNode(vertex));
        }
    }
}