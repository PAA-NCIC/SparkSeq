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

import java.util.List;
import java.util.ListIterator;

/**
 * Represents a route or path through a graph.
 * <p>
 *     In contrast with a {@link Path}, a route keeps track of the
 * path taken at furcations in order to speed up some path comparisions like the
 * one implemented by {@link #isSuffix}.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class Route<V extends BaseVertex, E extends BaseEdge> extends Path<V,E> {

    protected final Route<V,E> previousRouteWithLastVertexThatIsForkOrJoin;
    protected final boolean lastVertexIsForkOrJoin;

    /**
     * Create a zero length route with a start in a particular vertex:
     *
     * @param initialVertex the first vertex of the route.
     * @param graph the new route's graph.
     *
     * @throws IllegalArgumentException if {@code initialVertex} or {@code graph} are {@code null}.
     *      or if {@code initialVertex} does not belong to {@code graph}.
     */
    public Route(final V initialVertex, final BaseGraph<V, E> graph) {
        super(initialVertex, graph);
        previousRouteWithLastVertexThatIsForkOrJoin = null;
        lastVertexIsForkOrJoin = graph.inDegreeOf(initialVertex) > 1;
    }

    @Override
    public boolean equals(final Object other) {
        if (other == null) return false;
        if (other == this) return true;
        if (! (other instanceof Route)) return false;
        @SuppressWarnings("unchecked")
        final Route<V,E> otherRoute = (Route<V,E>) other;
        return otherRoute.length() == this.length() && isSuffix(otherRoute);
    }

    /**
     * Extends a route into a new instance.
     *
     * @param prefix the route to extend.
     * @param nextVertex the vertex to extend the route to.
     *
     * @throws IllegalArgumentException if {@code prefix} is {@code null} or {@code nextVertex} is {@code null}
     *   or {@code nextVertex} does not belong to {@code prefix}'s graph or there is no edge that in the graph
     *   that would connect {@code prefix}'s last vertex with {@code nextVertex} directly.
     */
    public Route(final Route<V,E> prefix, final V nextVertex) {
        this(prefix,resolveSuffixEdge(prefix,nextVertex));
    }


    /**
     * Extends a route into a new instance.
     *
     * @param prevVertex the vertex to extend the route to.
     * @param suffix the route to extend.
     *
     * @throws IllegalArgumentException if {@code suffix} is {@code null} or {@code prevVertex} is {@code null}
     *   or {@code prevVertex} does not belong to {@code suffix}'s graph or there is no edge that in the graph
     *   that would connect {@code suffix}'s first vertex with {@code prevVertex} directly.
     */
    public Route(final V prevVertex, final Route<V,E> suffix) {
        this(resolvePrefixEdge(prevVertex, suffix),suffix);
    }

    /**
     * Resolves the prefix edge as required by {@link Route(V,Route)}.
     */
    private static <V extends BaseVertex,E extends BaseEdge>  E resolvePrefixEdge(final V prevVertex, final Route<V, E> suffix) {
        if (prevVertex == null) throw new NullPointerException();
        if (!suffix.getGraph().containsVertex(prevVertex)) throw new IllegalArgumentException();
        final E result = suffix.getGraph().getEdge(prevVertex,suffix.getFirstVertex());
        if (result == null)
            throw new IllegalArgumentException("there is no such edge in the graph");
        return result;
    }

    /**
     * Resolves the suffix edge as required by {@link Route(Route,V)}
     */
    private static <V extends BaseVertex,E extends BaseEdge>  E resolveSuffixEdge(final Route<V,E> prefix, final V nextVertex) {
        if (nextVertex == null) throw new IllegalArgumentException();
        if (!prefix.getGraph().containsVertex(nextVertex)) throw new IllegalArgumentException();
        final E result = prefix.getGraph().getEdge(prefix.getLastVertex(),nextVertex);
        if (result == null)
            throw new IllegalArgumentException("there is no such edge in the graph");
        return result;
    }

    /**
     * Extends a route by prefixing an edge.
     *
     * @param initialEdge the extending edge.
     * @param suffix the original path.
     *
     * @throws IllegalArgumentException if {@code suffix} or {@code initialEdge} are {@code null}, or {@code initialEdge} is
     * not part of {@code suffix}'s graph, or {@code initialEdge} does not have as a target the first vertex in {@code suffix}.
     */
    public Route(final E initialEdge, final Route<V,E> suffix) {
        super(initialEdge,suffix);
        final V firstVertex = getFirstVertex();
        if(suffix.length() == 0) {
            lastVertexIsForkOrJoin = suffix.lastVertexIsForkOrJoin || graph.outDegreeOf(firstVertex) > 1;
            previousRouteWithLastVertexThatIsForkOrJoin = graph.inDegreeOf(firstVertex) > 1 ? new Route<>(firstVertex,graph) : null;
        } else {
            lastVertexIsForkOrJoin = suffix.lastVertexIsForkOrJoin;
            if (suffix.previousRouteWithLastVertexThatIsForkOrJoin != null)
                previousRouteWithLastVertexThatIsForkOrJoin = new Route<>(initialEdge,suffix.previousRouteWithLastVertexThatIsForkOrJoin);
            else
                previousRouteWithLastVertexThatIsForkOrJoin = graph.outDegreeOf(firstVertex) > 1 ?
                        new Route<>(new Route<>(firstVertex,graph),edgesInOrder.get(0)) :
                            graph.inDegreeOf(firstVertex) > 1 ? new Route<>(firstVertex,graph) : null;
        }
    }

    /**
     * Create copy of an existing route.
     * @param route the route to copy
     *
     * @throws NullPointerException if {@code route} is {@code null}.
     */
    protected Route(final Route<V, E> route) {
        super(route);
        lastVertexIsForkOrJoin = route.lastVertexIsForkOrJoin;
        previousRouteWithLastVertexThatIsForkOrJoin = route.previousRouteWithLastVertexThatIsForkOrJoin;
    }

    /**
     * Create a new Route extending another one with an edge
     *
     * @param route the route to extend.
     * @param edge the edge to extend the route with.
     *
     * @throws IllegalArgumentException if {@code route} or {@code edge} are {@code null}, or {@code edge} is
     * not part of {@code route}'s graph, or {@code edge} does not have as a source the last vertex in {@code route}.
     */
    public Route(final Route<V, E> route, final E edge) {
        super(route, edge);
        lastVertexIsForkOrJoin = graph.outDegreeOf(route.lastVertex) > 1 || graph.inDegreeOf(lastVertex) > 1;
        previousRouteWithLastVertexThatIsForkOrJoin = route.lastVertexIsForkOrJoin ? route : route.previousRouteWithLastVertexThatIsForkOrJoin;
    }

    @Override
    public boolean isSuffix(final Path<V,E> other) {
        if (other == this)
            return true;
        else if (other == null)
            throw new IllegalArgumentException("other path must not be null");
        else if (getGraph() != other.getGraph())
            throw new IllegalArgumentException("other path must be part of the same graph");
        else if (other instanceof Route)
            return isRouteSuffix((Route<V,E>)other);
        else
            return super.isSuffix(other);
    }

    @Override
    public String toString() {
        return super.toString().replace("Path{", "Route{");
    }

    /**
     * Faster version when comparing with a route.
     */
    protected boolean isRouteSuffix(final Route<V,E> other) {
        if (other.getGraph() != this.getGraph())
            throw new IllegalArgumentException("you cannot compare routes on different graphs");
        else if (lastVertex != other.lastVertex)  // obvious case.
            return false;
        else if (this.previousRouteWithLastVertexThatIsForkOrJoin == null
                && other.previousRouteWithLastVertexThatIsForkOrJoin != null) // I am shorter or different path for sure.
            return false;
        else if (this.edgesInOrder.size() < other.edgesInOrder.size())  // I am shorter regardless of path, no way Jose!
            return false;
        else if (this.previousRouteWithLastVertexThatIsForkOrJoin == null || other.previousRouteWithLastVertexThatIsForkOrJoin == null) {
            final ListIterator<E> myEdges = edgesInOrder.listIterator(edgesInOrder.size());
            final ListIterator<E> otherEdges = other.edgesInOrder.listIterator(other.edgesInOrder.size());
            while (otherEdges.hasPrevious())
                if (myEdges.previous() != otherEdges.previous())
                    return false;
            return true;
        } else
            return (other.previousRouteWithLastVertexThatIsForkOrJoin == this.previousRouteWithLastVertexThatIsForkOrJoin)
                || (previousRouteWithLastVertexThatIsForkOrJoin.lastVertex == other.previousRouteWithLastVertexThatIsForkOrJoin.lastVertex
              && previousRouteWithLastVertexThatIsForkOrJoin.isRouteSuffix(other.previousRouteWithLastVertexThatIsForkOrJoin));
    }

    /**
     * Checks whether the last vertex in the route is a fork or a joining vertex.
     * @return {@code true} iff so.
     */
    public boolean lastVertexIsForkOrJoin() {
        return lastVertexIsForkOrJoin;
    }

    /**
     * Returns the longest prefix route that has as a last vertex a join or furcation vertex.
     *
     * @return never {@code null}.
     */
    public Route<V,E> getPrefixRouteWithLastVertexThatIsForkOrJoin() {
        return previousRouteWithLastVertexThatIsForkOrJoin;
    }



    /**
     * Splice out the first few vertices of the route.
     *
     * @param length how many vertices to splice out
     * @return a new route without those spliced vertices.
     *
     * @throws IllegalArgumentException if {@code length} is equal to the route's length or greater or if it is negative.
     * Notice that non-vertex route are no legal routes.
     */
    public Route<V,E> splicePrefix(final int length) {
        if (length == 0)
            return this;
        if (length >= length())
            throw new IllegalArgumentException("prefix slicing to long");
        if (length < 0)
            throw new IllegalArgumentException("prefix cannot be negative");

        final List<E> resultEdges = getEdges().subList(length,length());
        Route<V,E> result = new Route<>(graph.getEdgeSource(resultEdges.get(0)),graph);
        for (final E edge : resultEdges)
            result = new Route<>(result,edge);
        return result;
    }
}
