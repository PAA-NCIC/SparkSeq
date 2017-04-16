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

import htsjdk.samtools.Cigar;
import org.apache.commons.lang.ArrayUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.CigarUtils;

import java.util.*;

/**
 * A path thought a BaseGraph
 *
 * class to keep track of paths
 *
 * User: depristo
 * Date: 3/19/13
 * Time: 2:34 PM
 *
 */
public class Path<T extends BaseVertex, E extends BaseEdge> {

    // the last vertex seen in the path
    protected final T lastVertex;

    // the list of edges comprising the path
    private Set<E> edgesAsSet = null;
    protected final ArrayList<E> edgesInOrder;

    // the scores for the path
    protected final int totalScore;

    // the graph from which this path originated
    protected final BaseGraph<T, E> graph;

    /**
     * Create a new Path containing no edges and starting at initialVertex
     * @param initialVertex the starting vertex of the path
     * @param graph the graph this path will follow through
     */
    public Path(final T initialVertex, final BaseGraph<T, E> graph) {
        if ( initialVertex == null ) throw new IllegalArgumentException("initialVertex cannot be null");
        if ( graph == null ) throw new IllegalArgumentException("graph cannot be null");
        if ( ! graph.containsVertex(initialVertex) ) throw new IllegalArgumentException("Vertex " + initialVertex + " must be part of graph " + graph);

        lastVertex = initialVertex;
        edgesInOrder = new ArrayList<>(0);
        totalScore = 0;
        this.graph = graph;
    }

    /**
     * Convenience constructor for testing that creates a path through vertices in graph
     */
    protected static <T extends BaseVertex, E extends BaseEdge> Path<T,E> makePath(final List<T> vertices, final BaseGraph<T, E> graph) {
        Path<T,E> path = new Path<T,E>(vertices.get(0), graph);
        for ( int i = 1; i < vertices.size(); i++ )
            path = new Path<T,E>(path, graph.getEdge(path.lastVertex, vertices.get(i)));
        return path;
    }

    /**
     * Create a new path with the same field values.
     *
     * @param p the template path.
     *
     * @throws NullPointerException if {@code p} is {@code null}.
     */
    protected Path(final Path<T,E> p) {
        this.edgesInOrder = p.edgesInOrder;
        this.lastVertex = p.lastVertex;
        this.edgesAsSet = p.edgesAsSet;
        this.totalScore = p.totalScore;
        this.graph = p.graph;
    }

    /**
     * Create a new Path extending p with edge
     *
     * @param p the path to extend.
     * @param edge the edge to extend path with.
     *
     * @throws IllegalArgumentException if {@code p} or {@code edge} are {@code null}, or {@code edge} is
     * not part of {@code p}'s graph, or {@code edge} does not have as a source the last vertex in {@code p}.
     */
    public Path(final Path<T,E> p, final E edge) {
        if ( p == null ) throw new IllegalArgumentException("Path cannot be null");
        if ( edge == null ) throw new IllegalArgumentException("Edge cannot be null");
        if ( ! p.graph.containsEdge(edge) ) throw new IllegalArgumentException("Graph must contain edge " + edge + " but it doesn't");
        if ( ! p.graph.getEdgeSource(edge).equals(p.lastVertex) ) { throw new IllegalStateException("Edges added to path must be contiguous."); }

        graph = p.graph;
        lastVertex = p.graph.getEdgeTarget(edge);
        edgesInOrder = new ArrayList<>(p.length() + 1);
        edgesInOrder.addAll(p.edgesInOrder);
        edgesInOrder.add(edge);
        totalScore = p.totalScore + edge.getMultiplicity();
    }

    /**
     * Length of the path in edges.
     *
     * @return {@code 0} or greater.
     */
    public int length() {
        return edgesInOrder.size();
    }

    /**
     * Prepend a path with an edge.
     *
     * @param edge the extending edge.
     * @param p the original path.
     *
     * @throws IllegalArgumentException if {@code p} or {@code edge} are {@code null}, or {@code edge} is
     * not part of {@code p}'s graph, or {@code edge} does not have as a target the first vertex in {@code p}.
     */
    public Path(final E edge, final Path<T,E> p) {
        if ( p == null ) throw new IllegalArgumentException("Path cannot be null");
        if ( edge == null ) throw new IllegalArgumentException("Edge cannot be null");
        if ( ! p.graph.containsEdge(edge) ) throw new IllegalArgumentException("Graph must contain edge " + edge + " but it doesn't");
        if ( ! p.graph.getEdgeTarget(edge).equals(p.getFirstVertex())) { throw new IllegalStateException("Edges added to path must be contiguous."); }
        graph = p.graph;
        lastVertex = p.lastVertex;
        edgesInOrder = new ArrayList<>(p.length() + 1);
        edgesInOrder.add(edge);
        edgesInOrder.addAll(p.getEdges());
        totalScore = p.totalScore + edge.getMultiplicity();
    }

    /**
     * Get the collection of edges leaving the last vertex of this path
     * @return a non-null collection
     */
    public Collection<E> getOutgoingEdgesOfLastVertex() {
        return getGraph().outgoingEdgesOf(getLastVertex());
    }

    /**
     * Does this path contain the given edge
     * @param edge  the given edge to test
     * @return      true if the edge is found in this path
     */
    public boolean containsEdge( final E edge ) {
        if( edge == null ) { throw new IllegalArgumentException("Attempting to test null edge."); }
        if ( edgesInOrder.isEmpty() ) return false;

        // initialize contains cache if necessary
        if ( edgesAsSet == null ) edgesAsSet = new HashSet<E>(edgesInOrder);
        return edgesAsSet.contains(edge);
    }

    /**
     * Does this path contain the given vertex?
     *
     * @param v a non-null vertex
     * @return true if v occurs within this path, false otherwise
     */
    public boolean containsVertex(final T v) {
        if ( v == null ) throw new IllegalArgumentException("Vertex cannot be null");

        // TODO -- warning this is expensive.  Need to do vertex caching
        return getVertices().contains(v);
    }

    /**
     * Checks whether a given path is a suffix of this path.
     *
     * @param other the path to compare against.
     * @throws IllegalArgumentException if <code>other</code> is <code>null</code>, or the come from
     *   different graphs.
     * @return true if <code>other</code> is a suffix of this path.
     */
    public boolean isSuffix(final Path<T, E> other) {
        if ( other == null ) throw new IllegalArgumentException("path cannot be null");
        if (other.getGraph() != this.getGraph()) throw new IllegalArgumentException("the other path most belong to the same path");
        if (!lastVertex.equals(other.lastVertex))
          return false;
        final ListIterator<E> myIt = edgesInOrder.listIterator(edgesInOrder.size());
        final ListIterator<E> otherIt = other.edgesInOrder.listIterator(other.edgesInOrder.size());
        while (myIt.hasPrevious() && otherIt.hasPrevious())
            if (otherIt.previous() != myIt.previous())
                return false;
        return !otherIt.hasPrevious();
    }

    /**
     * Check that two paths have the same edges and total score
     * @param path the other path we might be the same as
     * @return true if this and path are the same
     */
    protected boolean pathsAreTheSame(Path<T,E> path) {
        return totalScore == path.totalScore && edgesInOrder.equals(path.edgesInOrder);
    }

    @Override
    public String toString() {
        final StringBuilder b = new StringBuilder("Path{score=" + totalScore + ", path=");
        boolean first = true;
        for ( final T v : getVertices() ) {
            if ( first )
                first = false;
            else
                b.append(" -> ");
            b.append(v.getSequenceString());
        }
        b.append('}');
        return b.toString();
    }

    /**
     * Get the graph of this path
     * @return a non-null graph
     */
    public BaseGraph<T, E> getGraph() {
        return graph;
    }

    /**
     * Get the edges of this path in order
     * @return a non-null list of edges
     */
    public List<E> getEdges() { return edgesInOrder; }

    /**
     * Get the list of vertices in this path in order defined by the edges of the path
     * @return a non-null, non-empty list of vertices
     */
    public List<T> getVertices() {
        if ( getEdges().isEmpty() )
            return Collections.singletonList(lastVertex);
        else {
            final LinkedList<T> vertices = new LinkedList<T>();
            boolean first = true;
            for ( final E e : getEdges() ) {
                if ( first ) {
                    vertices.add(graph.getEdgeSource(e));
                    first = false;
                }
                vertices.add(graph.getEdgeTarget(e));
            }
            return vertices;
        }
    }

    /**
     * Get the total score of this path (bigger is better)
     * @return a positive integer
     */
    public int getScore() { return totalScore; }

    /**
     * Get the final vertex of the path
     * @return a non-null vertex
     */
    public T getLastVertex() { return lastVertex; }

    /**
     * Get the first vertex in this path
     * @return a non-null vertex
     */
    public T getFirstVertex() {
        if (edgesInOrder.size() == 0) {
            return lastVertex;
        } else {
            return getGraph().getEdgeSource(edgesInOrder.get(0));
        }
    }

    /**
     * The base sequence for this path. Pull the full sequence for source nodes and then the suffix for all subsequent nodes
     * @return  non-null sequence of bases corresponding to this path
     */
    public byte[] getBases() {
        if( getEdges().isEmpty() ) { return graph.getAdditionalSequence(lastVertex); }

        byte[] bases = graph.getAdditionalSequence(graph.getEdgeSource(edgesInOrder.get(0)));
        for( final E e : edgesInOrder ) {
            bases = ArrayUtils.addAll(bases, graph.getAdditionalSequence(graph.getEdgeTarget(e)));
        }
        return bases;
    }

    /**
     * Calculate the cigar elements for this path against the reference sequence
     *
     * @param refSeq the reference sequence that all of the bases in this path should align to
     * @return a Cigar mapping this path to refSeq, or null if no reasonable alignment could be found
     */
    public  Cigar calculateCigar(final byte[] refSeq) {
        return CigarUtils.calculateCigar(refSeq,getBases());
    }

    /**
     * Tests that this and other have the same score and vertices in the same order with the same seq
     * @param other the other path to consider.  Cannot be null
     * @return true if this and path are equal, false otherwise
     */
    public boolean equalScoreAndSequence(final Path<T,E> other) {
        if ( other == null ) throw new IllegalArgumentException("other cannot be null");
        return getScore() == other.getScore() && equalSequence(other);
    }

    /**
     * Tests that this and other have the same vertices in the same order with the same seq
     * @param other the other path to consider.  Cannot be null
     * @return true if this and path are equal, false otherwise
     */
    public boolean equalSequence(final Path<T,E> other) {
        final List<T> mine = getVertices();
        final List<T> yours = other.getVertices();
        if ( mine.size() == yours.size() ) { // hehehe
            for ( int i = 0; i < mine.size(); i++ )
                if ( ! mine.get(i).seqEquals(yours.get(i)) )
                    return false;
        }
        return true;
    }

}
