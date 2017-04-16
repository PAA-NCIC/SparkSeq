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

import java.util.Stack;

/**
 * A collection of route building methods.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class RouteFinder {


    /**
     * Completes a path backwards in the graph that would explain the sequence if bytes ending in the vertex provided.
     *
     * @param graph the graph to build the path upon.
     * @param sequence contains the sequence to backtrack.
     * @param start inclusive start position of the sequence to backtrack.
     * @param end exclusive end position of the sequence to backtrack.
     * @param vertex final vertex of the resulting path.
     * @return {@code null} if there is not such path, otherwise a path such that vertex is the last vertex of it
     * and its sequence is squence[start to end] + v.getSuffix();
     */
    private static <V extends BaseVertex, E extends BaseEdge> Route<V,E> extendRouteBackwards(final BaseGraph<V, E> graph,
                                                                                             final byte[] sequence,
                                                                                             final int start,
                                                                                             final int end,
                                                                                             final V vertex) {
        final Route<V,E> emptyPath = new Route<>(vertex,graph);
        if (end <= start) // trivial case.
            return emptyPath;
        final int kmerSize = graph.getKmerSize();
        final Stack<Pair<Route<V,E>,Integer>> stack = new Stack<>();
        stack.ensureCapacity(end - start + 1);
        stack.push(new Pair<>(emptyPath,end));
        while (!stack.isEmpty()) {
            final Pair<Route<V,E>,Integer> next = stack.pop();
            final Route<V,E> nextRoute = next.getFirst();
            final int nextEnd = next.getSecond();
            if (nextEnd <= start) {
                return nextRoute.splicePrefix(kmerSize - 1); // gotcha!!!
            }
            final V nextFirstVertex = nextRoute.getFirstVertex();
            if (graph.isSource(nextFirstVertex)) {
                final byte[] fullFirstVertexSequence = nextFirstVertex.getSequence();
                if (nextEnd - start != fullFirstVertexSequence.length - 1) {
                    continue; // you need to have the right length to accept a source vertex.
                }
                boolean mismatchFound = false;
                for (int i = 0; i < fullFirstVertexSequence.length - 1; i++) {
                    if (fullFirstVertexSequence[i] != sequence[i + start]) {
                        mismatchFound = true;
                        break;
                    }
                }
                if (!mismatchFound)
                    return nextRoute;
            } else {
                final Integer newNextEnd = nextEnd - 1;
                for (final E edge : graph.incomingEdgesOf(nextFirstVertex)) {
                    final V prevVertex = graph.getEdgeSource(edge);
                    final byte[] prevSequence = prevVertex.getSequence();
                    final byte prevByte = prevSequence[prevSequence.length - 1];
                    if (prevByte == sequence[newNextEnd]) {
                        stack.push(new Pair<>(new Route<>(edge,nextRoute),newNextEnd));
                    }
                }
            }
        }
        return null;
    }

    /**
     * Completes a path forward in the graph that would explain the sequence if bytes starting by the prefix provided.
     *
     * @param sequence missing sequence we want to
     * @param start inclusive first position in {@code sequence} that starts the extension
     * @param end exclusive position after the last of bases to be added to the extension.
     * @param prefix the seed prefix of the path.
     * @return {@code null} if there is not such path, otherwise a path such that vertex is the last vertex of it
     * and its sequence is prefix.getBases() + sequence[start to end];
     */
    private static  <V extends BaseVertex, E extends BaseEdge> Route<V,E> extendRouteForwards(
            final BaseGraph<V, E> graph, final byte[] sequence, final int start, final int end,
            final Route<V, E> prefix) {
        if (end <= start) // trivial case.
            return prefix;

        final Stack<Pair<Route<V,E>,Integer>> stack = new Stack<>();
        stack.ensureCapacity(end - start + 1);
        stack.push(new Pair<>(prefix,start));
        while (!stack.isEmpty()) {
            final Pair<Route<V,E>,Integer> next = stack.pop();
            final Route<V,E> nextRoute = next.getFirst();
            final int nextStart = next.getSecond();
            if (end <= nextStart)
                return nextRoute; // gotcha!!!
            final V lastVertex = nextRoute.getLastVertex();
            final Integer newNextStart = nextStart + 1;
                for (final E edge : graph.outgoingEdgesOf(lastVertex)) {
                    final V nextVertex = graph.getEdgeTarget(edge);
                    final byte[] nextSequence = nextVertex.getSequence();
                    final byte nextByte = nextSequence[nextSequence.length - 1];
                    if (nextByte == sequence[nextStart]) {
                        stack.push(new Pair<>(new Route<>(nextRoute,edge),newNextStart));
                    }
                }
        }
        return null;
    }

    /**
     * Construct a new route object give a sequence using unique kmer mappings.
     *
     * @param sequence base sequence.
     * @return {@code null} if there is no way such route on the graph or the start kmer is not unique.
     */
    @SuppressWarnings("unchecked")
    public static <V extends BaseVertex, E extends BaseEdge> Route<V,E> findRoute(final BaseGraph<V,E> graph,
                                                                                  final byte[] sequence) {
        if (graph == null)
            throw new NullPointerException();
        if (!(graph instanceof KmerSearchableGraph))
            throw new IllegalArgumentException("the input graph must implement " + KmerSearchableGraph.class.getName());

        final int kmerSize = graph.getKmerSize();
        final KmerSequence haplotypeKmers = new KmerSequence(sequence,kmerSize);

        if (haplotypeKmers.kmerSize() != graph.getKmerSize())
            throw new IllegalArgumentException("incompatible kmer sizes " + graph.getKmerSize() + " != " + haplotypeKmers.kmerSize());

        V vertex = null;
        int i;
        for (i = 0; i < haplotypeKmers.size(); i++)
            if ((vertex = ((KmerSearchableGraph<V,E>)graph).findKmer(haplotypeKmers.get(i))) != null)
                break;
        if (vertex == null)
            return null;
        if (!graph.containsVertex(vertex))
            throw new IllegalStateException("vertex does not belong to graph.");
        Route<V,E> result = i == 0 ? new Route<>(vertex,graph) :
                extendRouteBackwards(graph, sequence, 0, i + kmerSize - 1, vertex);
        if (result == null)
            return null;
        result = extendRouteForwards(graph, sequence, i + kmerSize, sequence.length, result);
        return result;
    }
}
