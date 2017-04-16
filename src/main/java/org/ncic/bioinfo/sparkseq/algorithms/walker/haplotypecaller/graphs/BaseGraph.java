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

import org.apache.commons.lang3.ArrayUtils;
import org.apache.log4j.Logger;
import org.jgrapht.EdgeFactory;
import org.jgrapht.alg.CycleDetector;
import org.jgrapht.graph.DefaultDirectedGraph;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: rpoplin
 * Date: 2/6/13
 */
public class BaseGraph<V extends BaseVertex, E extends BaseEdge> extends DefaultDirectedGraph<V, E> {
    protected final static Logger logger = Logger.getLogger(BaseGraph.class);
    protected final int kmerSize;

    /**
     * Construct a TestGraph with kmerSize
     * @param kmerSize
     */
    public BaseGraph(final int kmerSize, final EdgeFactory<V,E> edgeFactory) {
        super(edgeFactory);

        if ( kmerSize < 1 ) throw new IllegalArgumentException("kmerSize must be >= 1 but got " + kmerSize);
        this.kmerSize = kmerSize;
    }

    /**
     * How big of a kmer did we use to create this graph?
     * @return
     */
    public int getKmerSize() {
        return kmerSize;
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a reference node (meaning that it appears on the reference path in the graph)
     */
    public boolean isReferenceNode( final V v ) {
        if( v == null ) { throw new IllegalArgumentException("Attempting to test a null vertex."); }

        for ( final BaseEdge e : edgesOf(v) ) {
            if ( e.isRef() ) { return true; }
        }

        // edge case: if the graph only has one node then it's a ref node, otherwise it's not
        return (vertexSet().size() == 1);
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a source node (in degree == 0)
     */
    public boolean isSource( final V v ) {
        if( v == null ) { throw new IllegalArgumentException("Attempting to test a null vertex."); }
        return inDegreeOf(v) == 0;
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a sink node (out degree == 0)
     */
    public boolean isSink( final V v ) {
        if( v == null ) { throw new IllegalArgumentException("Attempting to test a null vertex."); }
        return outDegreeOf(v) == 0;
    }

    /**
     * Get the set of source vertices of this graph
     * @return a non-null set
     */
    public Set<V> getSources() {
        final Set<V> set = new LinkedHashSet<V>();
        for ( final V v : vertexSet() )
            if ( isSource(v) )
                set.add(v);
        return set;
    }

    /**
     * Get the set of sink vertices of this graph
     * @return a non-null set
     */
    public Set<V> getSinks() {
        final Set<V> set = new LinkedHashSet<V>();
        for ( final V v : vertexSet() )
            if ( isSink(v) )
                set.add(v);
        return set;
    }

    /**
     * Convert this kmer graph to a simple sequence graph.
     *
     * Each kmer suffix shows up as a distinct SeqVertex, attached in the same structure as in the kmer
     * graph.  Nodes that are sources are mapped to SeqVertex nodes that contain all of their sequence
     *
     * @return a newly allocated SequenceGraph
     */
    public SeqGraph convertToSequenceGraph() {

        final SeqGraph seqGraph = new SeqGraph(kmerSize);
        final Map<V, SeqVertex> vertexMap = new HashMap<>();


        // create all of the equivalent seq graph vertices
        for ( final V dv : vertexSet() ) {
            final SeqVertex sv = new SeqVertex(dv.getAdditionalSequence(isSource(dv)));
            sv.setAdditionalInfo(dv.additionalInfo());
            vertexMap.put(dv, sv);
            seqGraph.addVertex(sv);
        }

        // walk through the nodes and connect them to their equivalent seq vertices
        for( final E e : edgeSet() ) {
            final SeqVertex seqInV = vertexMap.get(getEdgeSource(e));
            final SeqVertex seqOutV = vertexMap.get(getEdgeTarget(e));
            //logger.info("Adding edge " + seqInV + " -> " + seqOutV);
            seqGraph.addEdge(seqInV, seqOutV, new BaseEdge(e.isRef(), e.getMultiplicity()));
        }

        return seqGraph;
    }

    /**
     * Pull out the additional sequence implied by traversing this node in the graph
     * @param v the vertex from which to pull out the additional base sequence
     * @return  non-null byte array
     */
    public byte[] getAdditionalSequence( final V v ) {
        if( v == null ) { throw new IllegalArgumentException("Attempting to pull sequence from a null vertex."); }
        return v.getAdditionalSequence(isSource(v));
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a reference source
     */
    public boolean isRefSource( final V v ) {
        if( v == null ) { throw new IllegalArgumentException("Attempting to test a null vertex."); }

        // confirm that no incoming edges are reference edges
        for ( final E edgeToTest : incomingEdgesOf(v) ) {
            if ( edgeToTest.isRef() ) { return false; }
        }

        // confirm that there is an outgoing reference edge
        for ( final E edgeToTest : outgoingEdgesOf(v) ) {
            if ( edgeToTest.isRef() ) { return true; }
        }

        // edge case: if the graph only has one node then it's a ref sink, otherwise it's not
        return (vertexSet().size() == 1);
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a reference sink
     */
    public boolean isRefSink( final V v ) {
        if( v == null ) { throw new IllegalArgumentException("Attempting to test a null vertex."); }

        // confirm that no outgoing edges are reference edges
        for ( final E edgeToTest : outgoingEdgesOf(v) ) {
            if ( edgeToTest.isRef() ) { return false; }
        }

        // confirm that there is an incoming reference edge
        for ( final E edgeToTest : incomingEdgesOf(v) ) {
            if ( edgeToTest.isRef() ) { return true; }
        }

        // edge case: if the graph only has one node then it's a ref source, otherwise it's not
        return (vertexSet().size() == 1);
    }

    /**
     * @return the reference source vertex pulled from the graph, can be null if it doesn't exist in the graph
     */
    public V getReferenceSourceVertex( ) {
        for( final V v : vertexSet() ) {
            if( isRefSource(v) ) {
                return v;
            }
        }
        return null;
    }

    /**
     * @return the reference sink vertex pulled from the graph, can be null if it doesn't exist in the graph
     */
    public V getReferenceSinkVertex( ) {
        for( final V v : vertexSet() ) {
            if( isRefSink(v) ) {
                return v;
            }
        }
        return null;
    }

    /**
     * Traverse the graph and get the next reference vertex if it exists
     * @param v the current vertex, can be null
     * @return  the next reference vertex if it exists, otherwise null
     */
    public V getNextReferenceVertex( final V v ) {
        return getNextReferenceVertex(v, false, Collections.<MultiSampleEdge>emptyList());
    }

    /**
     * Traverse the graph and get the next reference vertex if it exists
     * @param v the current vertex, can be null
     * @param allowNonRefPaths if true, allow sub-paths that are non-reference if there is only a single outgoing edge
     * @param blacklistedEdges edges to ignore in the traversal down; useful to exclude the non-reference dangling paths
     * @return the next vertex (but not necessarily on the reference path if allowNonRefPaths is true) if it exists, otherwise null
     */
    public V getNextReferenceVertex( final V v, final boolean allowNonRefPaths, final Collection<MultiSampleEdge> blacklistedEdges ) {
        if( v == null ) { return null; }

        // variable must be mutable because outgoingEdgesOf is an immutable collection
        Set<E> edges = outgoingEdgesOf(v);

        for( final E edgeToTest : edges ) {
            if( edgeToTest.isRef() ) {
                return getEdgeTarget(edgeToTest);
            }
        }

        // if we got here, then we aren't on a reference path
        if ( allowNonRefPaths ) {
            edges = new HashSet<>(edges);  // edges was immutable
            edges.removeAll(blacklistedEdges);
            if ( edges.size() == 1 )
                return getEdgeTarget(edges.iterator().next());
        }

        return null;
    }

    /**
     * Traverse the graph and get the previous reference vertex if it exists
     * @param v the current vertex, can be null
     * @return  the previous reference vertex if it exists
     */
    public V getPrevReferenceVertex( final V v ) {
        if( v == null ) { return null; }
        for( final E edgeToTest : incomingEdgesOf(v) ) {
            if( isReferenceNode(getEdgeSource(edgeToTest)) ) {
                return getEdgeSource(edgeToTest);
            }
        }
        return null;
    }

    /**
     * Does a reference path exist between the two vertices?
     * @param fromVertex    from this vertex, can be null
     * @param toVertex      to this vertex, can be null
     * @return              true if a reference path exists in the graph between the two vertices
     */
    public boolean referencePathExists(final V fromVertex, final V toVertex) {
        V v = fromVertex;
        if( v == null ) {
            return false;
        }
        v = getNextReferenceVertex(v);
        if( v == null ) {
            return false;
        }
        while( !v.equals(toVertex) ) {
            v = getNextReferenceVertex(v);
            if( v == null ) {
                return false;
            }
        }
        return true;
    }

    /**
     * Walk along the reference path in the graph and pull out the corresponding bases
     * @param fromVertex    starting vertex
     * @param toVertex      ending vertex
     * @param includeStart  should the starting vertex be included in the path
     * @param includeStop   should the ending vertex be included in the path
     * @return              byte[] array holding the reference bases, this can be null if there are no nodes between the starting and ending vertex (insertions for example)
     */
    public byte[] getReferenceBytes( final V fromVertex, final V toVertex, final boolean includeStart, final boolean includeStop ) {
        if( fromVertex == null ) { throw new IllegalArgumentException("Starting vertex in requested path cannot be null."); }
        if( toVertex == null ) { throw  new IllegalArgumentException("From vertex in requested path cannot be null."); }

        byte[] bytes = null;
        V v = fromVertex;
        if( includeStart ) {
            bytes = ArrayUtils.addAll(bytes, getAdditionalSequence(v));
        }
        v = getNextReferenceVertex(v); // advance along the reference path
        while( v != null && !v.equals(toVertex) ) {
            bytes = ArrayUtils.addAll(bytes, getAdditionalSequence(v));
            v = getNextReferenceVertex(v); // advance along the reference path
        }
        if( includeStop && v != null && v.equals(toVertex)) {
            bytes = ArrayUtils.addAll(bytes, getAdditionalSequence(v));
        }
        return bytes;
    }

    /**
     * Convenience function to add multiple vertices to the graph at once
     * @param vertices one or more vertices to add
     */
    public void addVertices(final V... vertices) {
        for ( final V v : vertices )
            addVertex(v);
    }

    /**
     * Convenience function to add multiple vertices to the graph at once
     * @param vertices one or more vertices to add
     */
    public void addVertices(final Collection<V> vertices) {
        for ( final V v : vertices )
            addVertex(v);
    }

    /**
     * Convenience function to add multiple edges to the graph
     * @param start the first vertex to connect
     * @param remaining all additional vertices to connect
     */
    public void addEdges(final V start, final V... remaining) {
        V prev = start;
        for ( final V next : remaining ) {
            addEdge(prev, next);
            prev = next;
        }
    }

    /**
     * Convenience function to add multiple edges to the graph
     * @param start the first vertex to connect
     * @param remaining all additional vertices to connect
     */
    public void addEdges(final E template, final V start, final V... remaining) {
        V prev = start;
        for ( final V next : remaining ) {
            addEdge(prev, next, (E)(template.copy())); // TODO -- is there a better way to do this?
            prev = next;
        }
    }

    /**
     * Get the set of vertices connected by outgoing edges of V
     * @param v a non-null vertex
     * @return a set of vertices connected by outgoing edges from v
     */
    public Set<V> outgoingVerticesOf(final V v) {
        final Set<V> s = new LinkedHashSet<V>();
        for ( final E e : outgoingEdgesOf(v) ) {
            s.add(getEdgeTarget(e));
        }
        return s;
    }

    /**
     * Get the set of vertices connected to v by incoming edges
     * @param v a non-null vertex
     * @return a set of vertices {X} connected X -> v
     */
    public Set<V> incomingVerticesOf(final V v) {
        final Set<V> s = new LinkedHashSet<V>();
        for ( final E e : incomingEdgesOf(v) ) {
            s.add(getEdgeSource(e));
        }
        return s;
    }

    /**
     * Get the set of vertices connected to v by incoming or outgoing edges
     * @param v a non-null vertex
     * @return a set of vertices {X} connected X -> v or v -> Y
     */
    public Set<V> neighboringVerticesOf(final V v) {
        final Set<V> s = incomingVerticesOf(v);
        s.addAll(outgoingVerticesOf(v));
        return s;
    }

    /**
     * Print out the graph in the dot language for visualization
     * @param destination File to write to
     */
    public void printGraph(final File destination, final int pruneFactor) {
        PrintStream stream = null;

        try {
            stream = new PrintStream(new FileOutputStream(destination));
            printGraph(stream, true, pruneFactor);
        } catch ( FileNotFoundException e ) {
            throw new RuntimeException(e);
        } finally {
            if ( stream != null ) stream.close();
        }
    }

    public void printGraph(final PrintStream graphWriter, final boolean writeHeader, final int pruneFactor) {
        if ( writeHeader )
            graphWriter.println("digraph assemblyGraphs {");

        for( final E edge : edgeSet() ) {
            graphWriter.println("\t" + getEdgeSource(edge).toString() + " -> " + getEdgeTarget(edge).toString() + " [" + (edge.getMultiplicity() > 0 && edge.getMultiplicity() <= pruneFactor ? "style=dotted,color=grey," : "") + "label=\"" + edge.getDotLabel() + "\"];");
            if( edge.isRef() ) {
                graphWriter.println("\t" + getEdgeSource(edge).toString() + " -> " + getEdgeTarget(edge).toString() + " [color=red];");
            }
        }

        for( final V v : vertexSet() ) {
//            graphWriter.println("\t" + v.toString() + " [label=\"" + v + "\",shape=box]");
            graphWriter.println("\t" + v.toString() + " [label=\"" + new String(getAdditionalSequence(v)) + v.additionalInfo() + "\",shape=box]");
        }

        if ( writeHeader )
            graphWriter.println("}");
    }

    /**
     * Remove edges that are connected before the reference source and after the reference sink
     *
     * Also removes all vertices that are orphaned by this process
     */
    public void cleanNonRefPaths() {
        if( getReferenceSourceVertex() == null || getReferenceSinkVertex() == null ) {
            return;
        }

        // Remove non-ref edges connected before and after the reference path
        final Set<E> edgesToCheck = new HashSet<E>();
        edgesToCheck.addAll(incomingEdgesOf(getReferenceSourceVertex()));
        while( !edgesToCheck.isEmpty() ) {
            final E e = edgesToCheck.iterator().next();
            if( !e.isRef() ) {
                edgesToCheck.addAll( incomingEdgesOf(getEdgeSource(e)) );
                removeEdge(e);
            }
            edgesToCheck.remove(e);
        }

        edgesToCheck.addAll(outgoingEdgesOf(getReferenceSinkVertex()));
        while( !edgesToCheck.isEmpty() ) {
            final E e = edgesToCheck.iterator().next();
            if( !e.isRef() ) {
                edgesToCheck.addAll( outgoingEdgesOf(getEdgeTarget(e)) );
                removeEdge(e);
            }
            edgesToCheck.remove(e);
        }

        removeSingletonOrphanVertices();
    }

    /**
     * Prune all chains from this graph where any edge in the path has multiplicity < pruneFactor
     *
     * @see LowWeightChainPruner for more information
     *
     * @param pruneFactor all edges with multiplicity < this factor that aren't ref edges will be removed
     */
    public void pruneLowWeightChains( final int pruneFactor ) {
        final LowWeightChainPruner<V,E> pruner = new LowWeightChainPruner<>(pruneFactor);
        pruner.pruneLowWeightChains(this);
    }

    /**
     * Remove all vertices in the graph that have in and out degree of 0
     */
    public void removeSingletonOrphanVertices() {
        // Run through the graph and clean up singular orphaned nodes
        final List<V> verticesToRemove = new LinkedList<>();
        for( final V v : vertexSet() ) {
            if( inDegreeOf(v) == 0 && outDegreeOf(v) == 0 && !isRefSource(v) ) {
                verticesToRemove.add(v);
            }
        }
        removeAllVertices(verticesToRemove);
    }

    /**
     * Remove all vertices on the graph that cannot be accessed by following any edge,
     * regardless of its direction, from the reference source vertex
     */
    public void removeVerticesNotConnectedToRefRegardlessOfEdgeDirection() {
        final HashSet<V> toRemove = new HashSet<>(vertexSet());

        final V refV = getReferenceSourceVertex();
        if ( refV != null ) {
            for ( final V v : new BaseGraphIterator<>(this, refV, true, true) ) {
                toRemove.remove(v);
            }
        }

        removeAllVertices(toRemove);
    }

    /**
     * Remove all vertices in the graph that aren't on a path from the reference source vertex to the reference sink vertex
     *
     * More aggressive reference pruning algorithm than removeVerticesNotConnectedToRefRegardlessOfEdgeDirection,
     * as it requires vertices to not only be connected by a series of directed edges but also prunes away
     * paths that do not also meet eventually with the reference sink vertex
     */
    public void removePathsNotConnectedToRef() {
        if ( getReferenceSourceVertex() == null || getReferenceSinkVertex() == null ) {
            throw new IllegalStateException("Graph must have ref source and sink vertices");
        }

        // get the set of vertices we can reach by going forward from the ref source
        final Set<V> onPathFromRefSource = new HashSet<>(vertexSet().size());
        for ( final V v : new BaseGraphIterator<>(this, getReferenceSourceVertex(), false, true) ) {
            onPathFromRefSource.add(v);
        }

        // get the set of vertices we can reach by going backward from the ref sink
        final Set<V> onPathFromRefSink = new HashSet<>(vertexSet().size());
        for ( final V v : new BaseGraphIterator<>(this, getReferenceSinkVertex(), true, false) ) {
            onPathFromRefSink.add(v);
        }

        // we want to remove anything that's not in both the sink and source sets
        final Set<V> verticesToRemove = new HashSet<>(vertexSet());
        onPathFromRefSource.retainAll(onPathFromRefSink);
        verticesToRemove.removeAll(onPathFromRefSource);
        removeAllVertices(verticesToRemove);

        // simple sanity checks that this algorithm is working.
        if ( getSinks().size() > 1 ) {
            throw new IllegalStateException("Should have eliminated all but the reference sink, but found " + getSinks());
        }

        if ( getSources().size() > 1 ) {
            throw new IllegalStateException("Should have eliminated all but the reference source, but found " + getSources());
        }
    }

    /**
     * Semi-lenient comparison of two graphs, truing true if g1 and g2 have similar structure
     *
     * By similar this means that both graphs have the same number of vertices, where each vertex can find
     * a vertex in the other graph that's seqEqual to it.  A similar constraint applies to the edges,
     * where all edges in g1 must have a corresponding edge in g2 where both source and target vertices are
     * seqEqual
     *
     * @param g1 the first graph to compare
     * @param g2 the second graph to compare
     * @param <T> the type of the nodes in those graphs
     * @return true if g1 and g2 are equals
     */
    public static <T extends BaseVertex, E extends BaseEdge> boolean graphEquals(final BaseGraph<T,E> g1, BaseGraph<T,E> g2) {
        final Set<T> vertices1 = g1.vertexSet();
        final Set<T> vertices2 = g2.vertexSet();
        final Set<E> edges1 = g1.edgeSet();
        final Set<E> edges2 = g2.edgeSet();

        if ( vertices1.size() != vertices2.size() || edges1.size() != edges2.size() )
            return false;

        for ( final T v1 : vertices1 ) {
            boolean found = false;
            for ( final T v2 : vertices2 )
                found = found || v1.getSequenceString().equals(v2.getSequenceString());
            if ( ! found ) return false;
        }

        for( final E e1 : g1.edgeSet() ) {
            boolean found = false;
            for( E e2 : g2.edgeSet() ) {
                if( g1.seqEquals(e1, e2, g2) ) { found = true; break; }
            }
            if( !found ) { return false; }
        }
        for( final E e2 : g2.edgeSet() ) {
            boolean found = false;
            for( E e1 : g1.edgeSet() ) {
                if( g2.seqEquals(e2, e1, g1) ) { found = true; break; }
            }
            if( !found ) { return false; }
        }
        return true;
    }

    // For use when comparing edges across graphs!
    private boolean seqEquals( final E edge1, final E edge2, final BaseGraph<V,E> graph2 ) {
        return (this.getEdgeSource(edge1).seqEquals(graph2.getEdgeSource(edge2))) && (this.getEdgeTarget(edge1).seqEquals(graph2.getEdgeTarget(edge2)));
    }


    /**
     * Get the incoming edge of v.  Requires that there be only one such edge or throws an error
     * @param v our vertex
     * @return the single incoming edge to v, or null if none exists
     */
    public E incomingEdgeOf(final V v) {
        return getSingletonEdge(incomingEdgesOf(v));
    }

    /**
     * Get the outgoing edge of v.  Requires that there be only one such edge or throws an error
     * @param v our vertex
     * @return the single outgoing edge from v, or null if none exists
     */
    public E outgoingEdgeOf(final V v) {
        return getSingletonEdge(outgoingEdgesOf(v));
    }

    /**
     * Helper function that gets the a single edge from edges, null if edges is empty, or
     * throws an error is edges has more than 1 element
     * @param edges a set of edges
     * @return a edge
     */
    private E getSingletonEdge(final Collection<E> edges) {
        if ( edges.size() > 1 ) throw new IllegalArgumentException("Cannot get a single incoming edge for a vertex with multiple incoming edges " + edges);
        return edges.isEmpty() ? null : edges.iterator().next();
    }

    /**
     * Add edge between source -> target if none exists, or add e to an already existing one if present
     *
     * @param source source vertex
     * @param target vertex
     * @param e edge to add
     */
    public void addOrUpdateEdge(final V source, final V target, final E e) {
        final E prev = getEdge(source, target);
        if ( prev != null ) {
            prev.add(e);
        } else {
            addEdge(source, target, e);
        }
    }

    @Override
    public String toString() {
        return "BaseGraph{" +
                "kmerSize=" + kmerSize +
                '}';
    }

    /**
     * Get the set of vertices within distance edges of source, regardless of edge direction
     *
     * @param source the source vertex to consider
     * @param distance the distance
     * @return a set of vertices within distance of source
     */
    protected Set<V> verticesWithinDistance(final V source, final int distance) {
        if ( distance == 0 )
            return Collections.singleton(source);

        final Set<V> found = new HashSet<>();
        found.add(source);
        for ( final V v : neighboringVerticesOf(source) ) {
            found.addAll(verticesWithinDistance(v, distance - 1));
        }

        return found;
    }

    /**
     * Get a graph containing only the vertices within distance edges of target
     * @param target a vertex in graph
     * @param distance the max distance
     * @return a non-null graph
     */
    public BaseGraph<V,E> subsetToNeighbors(final V target, final int distance) {
        if ( target == null ) throw new IllegalArgumentException("Target cannot be null");
        if ( ! containsVertex(target) ) throw new IllegalArgumentException("Graph doesn't contain vertex " + target);
        if ( distance < 0 ) throw new IllegalArgumentException("Distance must be >= 0 but got " + distance);


        final Set<V> toKeep = verticesWithinDistance(target, distance);
        final Set<V> toRemove = new HashSet<>(vertexSet());
        toRemove.removeAll(toKeep);

        final BaseGraph<V,E> result = (BaseGraph<V,E>)clone();
        result.removeAllVertices(toRemove);

        return result;
    }

    /**
     * Get a subgraph of graph that contains only vertices within 10 edges of the ref source vertex
     * @return a non-null subgraph of this graph
     */
    public BaseGraph<V,E> subsetToRefSource() {
        return subsetToNeighbors(getReferenceSourceVertex(), 10);
    }

    /**
     * Checks whether the graph contains all the vertices in a collection.
     *
     * @param vertices the vertices to check.
     *
     * @throws IllegalArgumentException if {@code vertices} is {@code null}.
     *
     * @return {@code true} if all the vertices in the input collection are present in this graph.
     * Also if the input collection is empty. Otherwise it returns {@code false}.
     */
    public boolean containsAllVertices(final Collection<? extends V> vertices) {
        if (vertices == null) throw new IllegalArgumentException("the input vertices collection cannot be null");
        for (final V vertex : vertices)
            if (!containsVertex(vertex)) return false;
        return true;
    }

    /**
     * Checks for the presence of directed cycles in the graph.
     *
     * @return {@code true} if the graph has cycles, {@code false} otherwise.
     */
    public boolean hasCycles() {
        return new CycleDetector<>(this).detectCycles();
    }
}
