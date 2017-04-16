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
 /**
 * Prune all chains from this graph where all edges in the path have multiplicity <= pruneFactor
 *
 * Unlike pruneGraph, this function will remove only linear chains in the graph where all edges have weight <= pruneFactor.
 *
 * For A -[1]> B -[1]> C -[1]> D would be removed with pruneFactor 1
 * but A -[1]> B -[2]> C -[1]> D would not be because the linear chain includes an edge with weight >= 2
 *
 * User: depristo
 * Date: 5/2/13
 * Time: 10:38 AM
 */
public class LowWeightChainPruner<V extends BaseVertex, E extends BaseEdge> {
    private final int pruneFactor;

    public LowWeightChainPruner(int pruneFactor) {
        if ( pruneFactor < 0 ) throw new IllegalArgumentException("pruneFactor must be >= 0 but got " + pruneFactor);
        this.pruneFactor = pruneFactor;
    }

    /**
     * Prune graph
     * @param graph the graph to prune
     */
    public void pruneLowWeightChains(final BaseGraph<V,E> graph) {
        if ( graph == null ) throw new IllegalArgumentException("Graph cannot be null");

        if ( pruneFactor > 0 ) {
            final Set<E> edgesToKeep = new LinkedHashSet<>();

            for ( final Path<V,E> linearChain : getLinearChains(graph) ) {
                if( mustBeKept(linearChain, pruneFactor) ) {
                    // we must keep edges in any path that contains a reference edge or an edge with weight > pruneFactor
                    edgesToKeep.addAll(linearChain.getEdges());
                }
            }

            // we want to remove all edges not in the keep set
            final Set<E> edgesToRemove = new HashSet<>(graph.edgeSet());
            edgesToRemove.removeAll(edgesToKeep);
            graph.removeAllEdges(edgesToRemove);

            graph.removeSingletonOrphanVertices();
        }
    }

    /**
     * Traverse the edges in the path and determine if any are either ref edges or have weight above or equal to
     * the pruning factor and should therefore not be pruned away.
     *
     * @param path the path in question
     * @param pruneFactor the integer pruning factor
     * @return true if any edge in the path must be kept
     */
    private boolean mustBeKept(final Path<V, E> path, final int pruneFactor) {
        for ( final E edge : path.getEdges() ) {
            if ( edge.getPruningMultiplicity() >= pruneFactor || edge.isRef() )
                return true;
        }
        return false;
    }

    /**
     * Get all of the linear chains in graph
     *
     * A linear chain is a series of vertices that start from either a source of a vertex with
     * out-degree > 1 and extend through all vertices accessible via an outgoing edge from this
     * vertex that have in == 1 and out degree of 0 or 1.
     *
     * @param graph the graph
     * @return a non-null collection of paths in graph
     */
    protected final Collection<Path<V,E>> getLinearChains(final BaseGraph<V,E> graph) {
        final Set<V> chainStarts = new LinkedHashSet<>();

        for ( final V v : graph.vertexSet() ) {
            // we want a list of all chain start vertices.  These are all vertices with out
            // degree > 1, or all source vertices.
            final int outDegree = graph.outDegreeOf(v);
            final int inDegree = graph.inDegreeOf(v);
            if ( outDegree > 1 || inDegree > 1 || (inDegree == 0 && outDegree > 0)) // don't add isolated vertices
                chainStarts.add(v);
        }

        // must be after since we can add duplicate starts in the above finding algorithm
        final List<Path<V, E>> linearChains = new LinkedList<>();
        for ( final V chainStart : chainStarts ) {
            for ( final E outEdge : graph.outgoingEdgesOf(chainStart) ) {
                // these chains are composed of the starts + their next vertices
                linearChains.add(extendLinearChain(new Path<>(new Path<>(chainStart, graph), outEdge)));
            }
        }

        return linearChains;
    }

    /**
     * Extend path while the last vertex has in and out degrees of 1 or 0
     * @param path the path to extend
     * @return a fully extended linear path
     */
    protected final Path<V,E> extendLinearChain(final Path<V, E> path) {
        final V last = path.getLastVertex();
        final Set<E> outEdges = path.getGraph().outgoingEdgesOf(last);

        final int outDegree = outEdges.size();
        final int inDegree = path.getGraph().inDegreeOf(last);

        if ( outDegree != 1 || inDegree > 1 ) {
            // out next vertex has multiple outgoing edges, so we are done with the linear path
            return path;
        } else {
            final V next = path.getGraph().getEdgeTarget(outEdges.iterator().next());
            if ( path.containsVertex(next) ) {
                // we are done if the path contains a cycle
                return path;
            } else {
                // we now know that last has outdegree == 1, so we keep extending the chain
                return extendLinearChain(new Path<>(path, outEdges.iterator().next()));
            }
        }
    }
}
