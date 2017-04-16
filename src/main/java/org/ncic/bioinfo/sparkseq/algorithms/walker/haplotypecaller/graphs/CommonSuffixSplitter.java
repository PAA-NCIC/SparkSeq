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

import java.util.*;

/**
 * Split a collection of middle nodes in a graph into their shared prefix and suffix values
 *
 * This code performs the following transformation.  Suppose I have a set of vertices V, such
 * that each vertex is composed of sequence such that
 *
 * Vi = prefix + seq_i + suffix
 *
 * where prefix and suffix are shared sequences across all vertices V.  This replaces each
 * Vi with three nodes prefix, seq_i, and suffix connected in a simple chain.
 *
 * This operation can be performed in a very general case, without too much worry about the incoming
 * and outgoing edge structure of each Vi.  The partner algorithm SharedSequenceMerger can
 * put these pieces back together in a smart way that maximizes the sharing of nodes
 * while respecting complex connectivity.
 *
 * User: depristo
 * Date: 3/22/13
 * Time: 8:31 AM
 */
public class CommonSuffixSplitter {
    /**
     * Create a new graph that contains the vertices in toMerge with their shared suffix and prefix
     * sequences extracted out.
     *
     */
    public CommonSuffixSplitter() {}

    /**
     * Simple single-function interface to split and then update a graph
     *
     * @param graph the graph containing the vertices in toMerge
     * @param v The bottom node whose incoming vertices we'd like to split
     * @return true if some useful splitting was done, false otherwise
     */
    public boolean split(final SeqGraph graph, final SeqVertex v) {
        if ( graph == null ) throw new IllegalArgumentException("graph cannot be null");
        if ( v == null ) throw new IllegalArgumentException("v cannot be null");
        if ( ! graph.vertexSet().contains(v) ) throw new IllegalArgumentException("graph doesn't contain vertex v " + v);

        final Collection<SeqVertex> toSplit = graph.incomingVerticesOf(v);
        if ( toSplit.size() < 2 )
            // Can only split at least 2 vertices
            return false;
        else if ( ! safeToSplit(graph, v, toSplit) ) {
            return false;
        } else {
            final SeqVertex suffixVTemplate = commonSuffix(toSplit);
            if ( suffixVTemplate.isEmpty() ) {
                return false;
            } else if ( wouldEliminateRefSource(graph, suffixVTemplate, toSplit) ) {
                return false;
            } else if ( allVerticesAreTheCommonSuffix(suffixVTemplate, toSplit) ) {
                return false;
            } else {
                final List<BaseEdge> edgesToRemove = new LinkedList<BaseEdge>();

//                graph.printGraph(new File("split.pre_" + v.getSequenceString() + "." + counter + ".dot"), 0);
                for ( final SeqVertex mid : toSplit ) {
                    // create my own copy of the suffix
                    final SeqVertex suffixV = new SeqVertex(suffixVTemplate.getSequence());
                    graph.addVertex(suffixV);
                    final SeqVertex prefixV = mid.withoutSuffix(suffixV.getSequence());
                    final BaseEdge out = graph.outgoingEdgeOf(mid);

                    final SeqVertex incomingTarget;
                    if ( prefixV == null ) {
                        // this node is entirely explained by suffix
                        incomingTarget = suffixV;
                    } else {
                        incomingTarget = prefixV;
                        graph.addVertex(prefixV);
                        graph.addEdge(prefixV, suffixV, new BaseEdge(out.isRef(), 1));
                        edgesToRemove.add(out);
                    }

                    graph.addEdge(suffixV, graph.getEdgeTarget(out), out.copy());

                    for ( final BaseEdge in : graph.incomingEdgesOf(mid) ) {
                        graph.addEdge(graph.getEdgeSource(in), incomingTarget, in.copy());
                        edgesToRemove.add(in);
                    }
                }

                graph.removeAllVertices(toSplit);
                graph.removeAllEdges(edgesToRemove);
//                graph.printGraph(new File("split.post_" + v.getSequenceString() + "." + counter++ + ".dot"), 0);

                return true;
            }
        }
    }

    /**
     * Would factoring out this suffix result in elimating the reference source vertex?
     * @param graph the graph
     * @param commonSuffix the common suffix of all toSplits
     * @param toSplits the list of vertices we're are trying to split
     * @return true if toSplit contains the reference source and this ref source has all and only the bases of commonSuffix
     */
    private boolean wouldEliminateRefSource(final SeqGraph graph, final SeqVertex commonSuffix, final Collection<SeqVertex> toSplits) {
        for ( final SeqVertex toSplit : toSplits ) {
            if ( graph.isRefSource(toSplit) )
                return toSplit.length() == commonSuffix.length();
        }
        return false;
    }

//    private static int counter = 0;

    /**
     * Would all vertices that we'd split just result in the common suffix?
     *
     * That is, suppose we have prefix nodes ABC and ABC.  After splitting all of the vertices would
     * just be ABC again, and we'd enter into an infinite loop.
     *
     * @param commonSuffix the common suffix of all vertices in toSplits
     * @param toSplits the collection of vertices we want to split
     * @return true if all of the vertices are equal to the common suffix
     */
    private boolean allVerticesAreTheCommonSuffix(final SeqVertex commonSuffix, final Collection<SeqVertex> toSplits) {
        for ( final SeqVertex toSplit : toSplits ) {
            if ( toSplit.length() != commonSuffix.length() )
                return false;
        }

        return true;
    }

    /**
     * Can we safely split up the vertices in toMerge?
     *
     * @param graph a graph
     * @param bot a vertex whose incoming vertices we want to split
     * @param toMerge the set of vertices we'd be splitting up
     * @return true if we can safely split up toMerge
     */
    private boolean safeToSplit(final SeqGraph graph, final SeqVertex bot, final Collection<SeqVertex> toMerge) {
        final Set<SeqVertex> outgoingOfBot = new HashSet<SeqVertex>(graph.outgoingVerticesOf(bot));
        for ( final SeqVertex m : toMerge ) {
            final Set<BaseEdge> outs = graph.outgoingEdgesOf(m);
            if ( m == bot || outs.size() != 1 || ! graph.outgoingVerticesOf(m).contains(bot) )
                // m == bot => don't allow self cycles in the graph
                return false;
            if ( outgoingOfBot.contains(m) )
                // forbid cycles from bottom -> mid
                return false;
        }

        return true;
    }

    /**
     * Return the longest suffix of bases shared among all provided vertices
     *
     * For example, if the vertices have sequences AC, CC, and ATC, this would return
     * a single C.  However, for ACC and TCC this would return CC.  And for AC and TG this
     * would return null;
     *
     * @param middleVertices a non-empty set of vertices
     * @return a single vertex that contains the common suffix of all middle vertices
     */
    protected static SeqVertex commonSuffix(final Collection<SeqVertex> middleVertices) {
        final List<byte[]> kmers = GraphUtils.getKmers(middleVertices);
        final int min = GraphUtils.minKmerLength(kmers);
        final int suffixLen = GraphUtils.compSuffixLen(kmers, min);
        final byte[] kmer = kmers.get(0);
        final byte[] suffix = Arrays.copyOfRange(kmer, kmer.length - suffixLen, kmer.length);
        return new SeqVertex(suffix);
    }
}