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

import java.util.*;

/**
 * Split a collection of middle nodes in a graph into their shared prefix and suffix values
 *
 * This code performs the following transformation.  Suppose I have a set of vertices V, such
 * that each vertex is composed of sequence such that
 *
 * Vi = prefix + seq_i + suffix
 *
 * where prefix and suffix are shared sequences across all vertices V
 *
 * This algorithm creates a new SeqGraph with the following configuration
 *
 * prefix -> has outgoing edges to all seq_i
 * suffix -> has incoming edges for all seq_i
 *
 * There are a few special cases that must be handled.  First, Vi could be simply
 * == to the prefix or the suffix.  These generate direct connections between
 * the prefix and suffix nodes, and they are handled internally by the algorithm.
 *
 * Note that for convenience, we will always create newTop and newBottom nodes, but
 * these may be empty node (i.e., they contain no sequence).  That allows them to be
 * trivially merged, if desired, when the graph is incorporated into an overall
 * graph.
 *
 * The product of this operation is a SeqGraph that contains the split.  There's a
 * function to merge reconnect this graph into the graph that contains the middle nodes
 *
 * The process guarentees a few things about the output:
 *
 * -- Preserves the paths and weights among all vertices
 *
 * It produces a graph that has some unusual properties
 *
 * -- May add nodes with no sequence (isEmpty() == true) to preserve connectivity among the graph
 * -- May introduce edges with no multiplicity to preserve paths through the graph
 *
 * The overall workflow of using this class is simple:
 *
 * find vertices V in graph that you want to split out
 * s = new SharedVertexSequenceSplitter(graph, V)
 * s.updateGraph(graph)
 *
 * to update the graph with the modifications created by this splitter
 *
 * User: depristo
 * Date: 3/22/13
 * Time: 8:31 AM
 */
public class SharedVertexSequenceSplitter {
    final private SeqGraph outer;
    final protected SeqVertex prefixV, suffixV;
    final protected Collection<SeqVertex> toSplits;

    // updated in split routine
    protected SeqGraph splitGraph = null;
    protected Collection<SeqVertex> newMiddles = null;
    protected List<BaseEdge> edgesToRemove = null;

    /**
     * Create a new graph that contains the vertices in toSplitsArg with their shared suffix and prefix
     * sequences extracted out.
     *
     * @param graph the graph containing the vertices in toSplitsArg
     * @param toSplitsArg a collection of vertices to split.  Must be contained within graph, and have only connections
     *                    from a single shared top and/or bottom node
     */
    public SharedVertexSequenceSplitter(final SeqGraph graph, final Collection<SeqVertex> toSplitsArg) {
        if ( graph == null ) throw new IllegalArgumentException("graph cannot be null");
        if ( toSplitsArg == null ) throw new IllegalArgumentException("toSplitsArg cannot be null");
        if ( toSplitsArg.size() < 2 ) throw new IllegalArgumentException("Can only split at least 2 vertices but only got " + toSplitsArg);
        if ( ! graph.vertexSet().containsAll(toSplitsArg) ) throw new IllegalArgumentException("graph doesn't contain all of the vertices to split");

        this.outer = graph;
        this.toSplits = toSplitsArg;

        // all of the edges point to the same sink, so it's time to merge
        final Pair<SeqVertex, SeqVertex> prefixAndSuffix = commonPrefixAndSuffixOfVertices(toSplits);
        prefixV = prefixAndSuffix.getFirst();
        suffixV = prefixAndSuffix.getSecond();
    }

    /**
     * Given sequencing that are all equal, does this splitter make those into prefix or suffix nodes?
     * @return true if we merge equal nodes into prefix nodes or suffix nodes
     */
    protected static boolean prefersPrefixMerging() {
        return true;
    }

    /**
     * Simple single-function interface to split and then update a graph
     *
     * @see #updateGraph(SeqVertex, SeqVertex) for a full description of top and bottom
     *
     * @param top the top vertex, may be null
     * @param bottom the bottom vertex, may be null
     * @return true if some useful splitting was done, false otherwise
     */
    public boolean splitAndUpdate(final SeqVertex top, final SeqVertex bottom) {
        split();
        updateGraph(top, bottom);
        return true;
    }

    /**
     * Does either the common suffix or prefix have at least minCommonSequence bases in it?
     * @param minCommonSequence a minimum length of the common sequence, must be >= 0
     * @return true if either suffix or prefix length >= minCommonSequence
     */
    public boolean meetsMinMergableSequenceForEitherPrefixOrSuffix(final int minCommonSequence) {
        return meetsMinMergableSequenceForPrefix(minCommonSequence) || meetsMinMergableSequenceForSuffix(minCommonSequence);
    }

    /**
     * Does the common prefix have at least minCommonSequence bases in it?
     * @param minCommonSequence a minimum length of the common sequence, must be >= 0
     * @return true if prefix length >= minCommonSequence
     */
    public boolean meetsMinMergableSequenceForPrefix(final int minCommonSequence) {
        return prefixV.length() >= minCommonSequence;
    }

    /**
     * Does the common suffix have at least minCommonSequence bases in it?
     * @param minCommonSequence a minimum length of the common sequence, must be >= 0
     * @return true if suffix length >= minCommonSequence
     */
    public boolean meetsMinMergableSequenceForSuffix(final int minCommonSequence) {
        return suffixV.length() >= minCommonSequence;
    }

    /**
     * Actually do the splitting up of the vertices
     *
     * Must be called before calling updateGraph
     */
    public void split() {
        splitGraph = new SeqGraph(outer.getKmerSize());
        newMiddles = new LinkedList<SeqVertex>();
        edgesToRemove = new LinkedList<BaseEdge>();

        splitGraph.addVertices(prefixV, suffixV);

        for ( final SeqVertex mid : toSplits ) {
            final BaseEdge toMid = processEdgeToRemove(mid, outer.incomingEdgeOf(mid));
            final BaseEdge fromMid = processEdgeToRemove(mid, outer.outgoingEdgeOf(mid));

            final SeqVertex remaining = mid.withoutPrefixAndSuffix(prefixV.getSequence(), suffixV.getSequence());
            if ( remaining != null ) {
                // there's some sequence prefix + seq + suffix, so add the node and make edges
                splitGraph.addVertex(remaining);
                newMiddles.add(remaining);
                // update edge from top -> middle to be top -> without suffix
                splitGraph.addEdge(prefixV, remaining, toMid);
                splitGraph.addEdge(remaining, suffixV, fromMid);
            } else {
                // prefix + suffix completely explain this node
                splitGraph.addOrUpdateEdge(prefixV, suffixV, toMid.copy().add(fromMid));
            }
        }
    }

    /**
     * Update graph outer, replacing the previous middle vertices that were split out with the new
     * graph structure of the split, linking this subgraph into the graph at top and bot (the
     * vertex connecting the middle nodes and the vertex outgoing of all middle node)
     *
     * @param top an optional top node that must have outgoing edges to all split vertices.  If null, this subgraph
     *            will be added without any incoming edges
     * @param bot an optional bottom node that must have incoming edges to all split vertices.  If null, this subgraph
     *            will be added without any outgoing edges to the rest of the graph
     */
    public void updateGraph(final SeqVertex top, final SeqVertex bot) {
        if ( ! outer.vertexSet().containsAll(toSplits) ) throw new IllegalArgumentException("graph doesn't contain all of the original vertices to split");
        if ( top == null && bot == null ) throw new IllegalArgumentException("Cannot update graph without at least one top or bot vertex, but both were null");
        if ( top != null && ! outer.containsVertex(top) ) throw new IllegalArgumentException("top " + top + " not found in graph " + outer);
        if ( bot != null && ! outer.containsVertex(bot) ) throw new IllegalArgumentException("bot " + bot + " not found in graph " + outer);
        if ( splitGraph == null ) throw new IllegalStateException("Cannot call updateGraph until split() has been called");

        outer.removeAllVertices(toSplits);
        outer.removeAllEdges(edgesToRemove);

        outer.addVertices(newMiddles);

        final boolean hasPrefixSuffixEdge = splitGraph.getEdge(prefixV, suffixV) != null;
        final boolean hasOnlyPrefixSuffixEdges = hasPrefixSuffixEdge && splitGraph.outDegreeOf(prefixV) == 1;
        final boolean needPrefixNode = ! prefixV.isEmpty() || (top == null && ! hasOnlyPrefixSuffixEdges);
        final boolean needSuffixNode = ! suffixV.isEmpty() || (bot == null && ! hasOnlyPrefixSuffixEdges);

        // if prefix / suffix are needed, keep them
        final SeqVertex topForConnect = needPrefixNode ? prefixV : top;
        final SeqVertex botForConnect = needSuffixNode ? suffixV : bot;

        if ( needPrefixNode ) {
            outer.addVertex(prefixV);
            if ( top != null ) outer.addEdge(top, prefixV, BaseEdge.orRef(splitGraph.outgoingEdgesOf(prefixV), 1));
        }

        if ( needSuffixNode ) {
            outer.addVertex(suffixV);
            if ( bot != null ) outer.addEdge(suffixV, bot, BaseEdge.orRef(splitGraph.incomingEdgesOf(suffixV), 1));
        }

        if ( topForConnect != null ) {
            for ( final BaseEdge e : splitGraph.outgoingEdgesOf(prefixV) ) {
                final SeqVertex target = splitGraph.getEdgeTarget(e);

                if ( target == suffixV ) { // going straight from prefix -> suffix
                    if ( botForConnect != null )
                        outer.addEdge(topForConnect, botForConnect, e);
                } else {
                    outer.addEdge(topForConnect, target, e);
                }
            }
        }

        if ( botForConnect != null ) {
            for ( final BaseEdge e : splitGraph.incomingEdgesOf(suffixV) ) {
                outer.addEdge(splitGraph.getEdgeSource(e), botForConnect, e);
            }
        }
    }

    /**
     * Return the longest suffix of bases shared among all provided vertices
     *
     * For example, if the vertices have sequences AC, CC, and ATC, this would return
     * a single C.  However, for ACC and TCC this would return CC.  And for AC and TG this
     * would return null;
     *
     * @param middleVertices a non-empty set of vertices
     * @return
     */
    protected static Pair<SeqVertex, SeqVertex> commonPrefixAndSuffixOfVertices(final Collection<SeqVertex> middleVertices) {
        final List<byte[]> kmers = new ArrayList<byte[]>(middleVertices.size());

        int min = Integer.MAX_VALUE;
        for ( final SeqVertex v : middleVertices ) {
            kmers.add(v.getSequence());
            min = Math.min(min, v.getSequence().length);
        }

        final int prefixLen = GraphUtils.compPrefixLen(kmers, min);
        final int suffixLen = GraphUtils.compSuffixLen(kmers, min - prefixLen);

        final byte[] kmer = kmers.get(0);
        final byte[] prefix = Arrays.copyOfRange(kmer, 0, prefixLen);
        final byte[] suffix = Arrays.copyOfRange(kmer, kmer.length - suffixLen, kmer.length);
        return new Pair<SeqVertex, SeqVertex>(new SeqVertex(prefix), new SeqVertex(suffix));
    }

    /**
     * Helper function that returns an edge that we should use for splitting
     *
     * If e is null, creates a new 0 multiplicity edge, set to ref is any edges to V are ref
     * If e is not null, returns a new copy of e, and schedules e for removal
     *
     * @param e a non-null edge
     * @return a non-null edge
     */
    private BaseEdge processEdgeToRemove(final SeqVertex v, final BaseEdge e) {
        if ( e == null ) {
            // there's no edge, so we return a newly allocated one and don't schedule e for removal
            // the weight must be 0 to preserve sum through the diamond
            return new BaseEdge(outer.isReferenceNode(v), 0);
        } else {
            // schedule edge for removal, and return a freshly allocated one for our graph to use
            edgesToRemove.add(e);
            return e.copy();
        }
    }
}