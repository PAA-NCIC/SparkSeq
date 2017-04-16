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

import org.jgrapht.alg.CycleDetector;
import org.ncic.bioinfo.sparkseq.algorithms.data.basic.Pair;
import org.ncic.bioinfo.sparkseq.algorithms.utils.haplotype.Haplotype;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.*;

/**
 * Efficient algorithm to obtain the list of best haplotypes given the {@link SeqGraph instace}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class KBestHaplotypeFinder extends AbstractList<KBestHaplotype> implements Iterable<KBestHaplotype>  {

    /**
     * The search graph.
     */
    private final SeqGraph graph;

    /**
     * Map of sub-haplotype finder by their source vertex.
     */
    protected Map<SeqVertex,KBestSubHaplotypeFinder> finderByVertex;

    /**
     * Possible haplotype sink vertices.
     */
    protected Set<SeqVertex> sinks;

    /**
     * Possible haplotype source vertices.
     */
    protected Set<SeqVertex> sources;

    /**
     * The top finder.
     *
     * <p>If there is only a single source vertex, its finder is the top finder. However whent there
     * is more than one possible source, we create a composite finder that alternates between individual source vertices
     * for their best haplotypes.</p>
     */
    private final KBestSubHaplotypeFinder topFinder;

    /**
     * Constructs a new best haplotypes finder.
     *
     * @param graph the seq-graph to search.
     * @param source the source vertex for all haplotypes.
     * @param sink sink vertices for all haplotypes.
     *
     * @throws IllegalArgumentException if <ul>
     *     <li>any of {@code graph}, {@code source} or {@code sink} is {@code null} or</li>
     *     <li>either {@code source} or {@code sink} is not a vertex in {@code graph}.</li>
     * </ul>
     */
    public KBestHaplotypeFinder(final SeqGraph graph, final SeqVertex source, final SeqVertex sink) {
        this(graph,Collections.singleton(source),Collections.singleton(sink));
    }

    /**
     * Constructs a new best haplotypes finder.
     *
     * @param graph the seq-graph to search.
     * @param sources source vertices for all haplotypes.
     * @param sinks sink vertices for all haplotypes.
     *
     * @throws IllegalArgumentException if <ul>
     *     <li>any of {@code graph}, {@code sources} or {@code sinks} is {@code null} or</li>
     *     <li>any of {@code sources}' or any {@code sinks}' member is not a vertex in {@code graph}.</li>
     * </ul>
     */
    public KBestHaplotypeFinder(final SeqGraph graph, final Set<SeqVertex> sources, final Set<SeqVertex> sinks) {
        if (graph == null) throw new IllegalArgumentException("graph cannot be null");
        if (sources == null) throw new IllegalArgumentException("source cannot be null");
        if (sinks == null) throw new IllegalArgumentException("sink cannot be null");
        if (!graph.containsAllVertices(sources)) throw new IllegalArgumentException("source does not belong to the graph");
        if (!graph.containsAllVertices(sinks)) throw new IllegalArgumentException("sink does not belong to the graph");

        //TODO dealing with cycles here due to a bug in some of the graph transformations that produces cycles.
        //TODO Once that is solve, the if-else below should be substituted by a throw if there is any cycles,
        //TODO just the line commented out below if you want to trade early-bug-fail for speed.
        //this.graph = graph;
        this.graph = new CycleDetector<>(graph).detectCycles() ? removeCycles(graph,sources,sinks) : graph;

        finderByVertex = new HashMap<>(this.graph.vertexSet().size());
        this.sinks = sinks;
        this.sources = sources;
        if (sinks.size() == 0 || sources.size() == 0)
            topFinder = DeadEndKBestSubHaplotypeFinder.INSTANCE;
        else if (sources.size() == 1)
            topFinder = createVertexFinder(sources.iterator().next());
        else
            topFinder = createAggregatedFinder();
    }

    /**
     * Constructs a new best haplotype finder.
     * <p>
     *     It will consider all source and sink vertex when looking for haplotypes.
     * </p>
     *
     * @param graph the seq-graph to search for the best haplotypes.
     */
    public KBestHaplotypeFinder(SeqGraph graph) {
        this(graph,graph.getSources(),graph.getSinks());
    }

    /**
     * Creates an aggregated recursive finder to try all possible source vertices.
     *
     * @return never {@code null}.
     */
    private KBestSubHaplotypeFinder createAggregatedFinder() {
        final List<KBestSubHaplotypeFinder> sourceFinders = new ArrayList<>(sources.size());
        for (final SeqVertex source : sources)
            sourceFinders.add(createVertexFinder(source));
        return new AggregatedSubHaplotypeFinder(sourceFinders);
    }

    /**
     * Removes edges that produces cycles and also dead vertices that do not lead to any sink vertex.
     *
     * @param original graph to modify.
     * @param sources considered source vertices.
     * @param sinks considered sink vertices.
     * @return never {@code null}.
     */
    private static SeqGraph removeCycles(final SeqGraph original, final Set<SeqVertex> sources, final Set<SeqVertex> sinks) {
        final Set<BaseEdge> edgesToRemove = new HashSet<>(original.edgeSet().size());
        final Set<SeqVertex> vertexToRemove = new HashSet<>(original.vertexSet().size());

        boolean foundSomePath = false;
        for (final SeqVertex source : sources)
            foundSomePath = findGuiltyVerticesAndEdgesToRemoveCycles(original, source, sinks, edgesToRemove,
                    vertexToRemove, new HashSet<SeqVertex>(original.vertexSet().size())) | foundSomePath;

        if (!foundSomePath)
            throw new IllegalStateException("could not find any path from the source vertex to the sink vertex after removing cycles: "
                    + Arrays.toString(sources.toArray()) + " => " + Arrays.toString(sinks.toArray()));

        if (edgesToRemove.isEmpty() && vertexToRemove.isEmpty())
            throw new IllegalStateException("cannot find a way to remove the cycles");

        final SeqGraph result = (SeqGraph) original.clone();
        result.removeAllEdges(edgesToRemove);
        result.removeAllVertices(vertexToRemove);
        return result;
    }

    /**
     * Recursive call that looks for edges and vertices that need to be removed to get rid of cycles.
     *
     * @param graph the original graph.
     * @param currentVertex current search vertex.
     * @param sinks considered sink vertices.
     * @param edgesToRemove collection  of edges that need to be removed in order to get rid of cycles.
     * @param verticesToRemove collection of vertices that can be removed.
     * @param parentVertices collection of vertices that preceded the {@code currentVertex}; i.e. the it can be
     *                       reached from those vertices using edges existing in {@code graph}.
     *
     * @return {@code true} to indicate that the some sink vertex is reachable by {@code currentVertex},
     *  {@code false} otherwise.
     */
    private static boolean findGuiltyVerticesAndEdgesToRemoveCycles(final SeqGraph graph,
                                                                    final SeqVertex currentVertex,
                                                                    final Set<SeqVertex> sinks,
                                                                    final Set<BaseEdge> edgesToRemove,
                                                                    final Set<SeqVertex> verticesToRemove,
                                                                    final Set<SeqVertex> parentVertices) {
        if (sinks.contains(currentVertex)) return true;

        final Set<BaseEdge> outgoingEdges = graph.outgoingEdgesOf(currentVertex);
        boolean reachesSink = false;
        parentVertices.add(currentVertex);

        for (final BaseEdge edge : outgoingEdges) {
            final SeqVertex child = graph.getEdgeTarget(edge);
            if (parentVertices.contains(child))
                edgesToRemove.add(edge);
            else {
              final boolean childReachSink = findGuiltyVerticesAndEdgesToRemoveCycles(graph, child, sinks,
                      edgesToRemove, verticesToRemove, parentVertices);
              reachesSink = reachesSink || childReachSink;
            }
        }
        parentVertices.remove(currentVertex);
        if (!reachesSink) verticesToRemove.add(currentVertex);
        return reachesSink;
    }

    @Override
    public KBestHaplotype get(final int index) {
        if (index < 0 || index >= size())
            throw new IndexOutOfBoundsException();
        return topFinder.getKBest(index);
    }

    @Override
    public Iterator<KBestHaplotype> iterator() {
        return new Iterator<KBestHaplotype>() {
            private int nextK = 0;
            private final int maxK = topFinder.getCount();


            @Override
            public boolean hasNext() {
                return nextK < maxK;
            }

            @Override
            public KBestHaplotype next() {
                if (nextK >= maxK) throw new NoSuchElementException();
                return topFinder.getKBest(nextK++);
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    @Override
    public int size() {
        return topFinder.getCount();
    }

    /**
     * Returns an iterator on the first k best haplotypes.
     * <p>
     *     It might return less than k haplotypes if the total number of possible haplotypes is smaller.
     * </p>
     *
     * @param k the maximum number of haplotypes to return.
     * @return never {@code null}, but perhaps a iterator that return no haplotype.
     */
    public Iterator<KBestHaplotype> iterator(final int k) {

        return new Iterator<KBestHaplotype>() {
            private int nextK = 0;
            private final int maxK = Math.min(size(), k);

            @Override
            public boolean hasNext() {
                return nextK < maxK;
            }

            @Override
            public KBestHaplotype next() {
                if (nextK >= maxK) throw new NoSuchElementException();
                return topFinder.getKBest(nextK++);
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    /**
     * Creates a finder from a vertex.
     *
     * @param vertex the source vertex for the finder.
     *
     * @return never {@code null}, perhaps a finder that returns no haplotypes though.
     */
    protected KBestSubHaplotypeFinder createVertexFinder(final SeqVertex vertex) {
        KBestSubHaplotypeFinder finder = finderByVertex.get(vertex);
        if (finder == null) {
            if (sinks.contains(vertex))
                finder = new EmptyPathHaplotypeFinderNode(graph,vertex);
            else {
                final Set<BaseEdge> outgoingEdges = graph.outgoingEdgesOf(vertex);
                if (outgoingEdges.isEmpty())
                    finder = DeadEndKBestSubHaplotypeFinder.INSTANCE;
                else {
                    final Map<BaseEdge,KBestSubHaplotypeFinder> undeadChildren = createChildrenFinders(outgoingEdges);
                    finder = undeadChildren.isEmpty() ? DeadEndKBestSubHaplotypeFinder.INSTANCE :
                            new RecursiveSubHaplotypeFinder(graph,vertex,undeadChildren);
                }
            }
            finderByVertex.put(vertex, finder);
        }
        return finder;
    }

    /**
     * Creates finder for target vertices of a collection of edges.
     * <p>
     *     This peculiar signature is convenient for when we want to create finders for the children of a vertex.
     * </p>
     *
     * @param baseEdges target collection of edges.
     *
     * @return never {@code null}, perhaps an empty map if there is no children with valid paths to any sink for this
     *  finder.
     */
    private Map<BaseEdge, KBestSubHaplotypeFinder> createChildrenFinders(final Set<BaseEdge> baseEdges) {
        final Map<BaseEdge,KBestSubHaplotypeFinder> result = new LinkedHashMap<>(baseEdges.size());
        for (final BaseEdge edge : baseEdges) {
            final KBestSubHaplotypeFinder targetFinder = createVertexFinder(graph.getEdgeTarget(edge));
            if (targetFinder.getCount() == 0) continue;
            result.put(edge, targetFinder);
        }
        return result;
    }

    /**
     * Print a DOT representation of search graph.
     *
     * @param out character stream printer where to print the DOT representation to.
     *
     * @throws IllegalArgumentException if {@code out} is {@code null}.
     */
    public void printDOT(final PrintWriter out) {
        if (out == null)
            throw new IllegalArgumentException("the out writer cannot be null");
        out.println("digraph {");
        out.println("    rankdir = LR");
        out.println("    node [shape=box, margin=0.01]");
        out.println("    subgraph cluster_dummy { style = invis; x [label=\"\",shape=none,margin=0] }");
        final StringBuilder referenceCluster = new StringBuilder(1000);

        referenceCluster.append("    subgraph cluster_ref {\n");
        referenceCluster.append("        node [penwidth=2]\n");
        for (final KBestSubHaplotypeFinder finder : finderByVertex.values() ) {
            final String id = finder.id();
            final String line = String.format("    %s [label=<%s>]",id,finder.label());
            if (finder.isReference())
                referenceCluster.append("    ").append(line).append('\n');
            else
                out.println(line);
        }
        referenceCluster.append("    }");
        out.println(referenceCluster.toString());

        for (final KBestSubHaplotypeFinder finder : finderByVertex.values())
            for (final Pair<? extends KBestSubHaplotypeFinder,String> subFinderLabel : finder.subFinderLabels()) {
                final KBestSubHaplotypeFinder subFinder = subFinderLabel.getFirst();

                final String edgeLabel = subFinderLabel.getSecond();
                out.println(String.format("    %s -> %s [label=%s]",finder.id(),subFinder.id(),edgeLabel));
            }
        out.println("}");
    }

    /**
     * Print a DOT representation of search graph.
     *
     * @param file file where to print the DOT representation to.
     *
     * @throws IllegalArgumentException if {@code file} is {@code null}.
     * @throws FileNotFoundException if {@code file} cannot be created or written.
     * @throws IllegalStateException if there was some trouble when writing the DOT representation.
     */
    public void printDOT(final File file) throws FileNotFoundException {
        if (file == null)
            throw new IllegalArgumentException("the output file cannot be null");
        final PrintWriter out = new PrintWriter(file);
        printDOT(out);
        if (out.checkError())
            throw new IllegalStateException("error occurred while writing k-best haplotype search graph into file '"
                    + file.getAbsolutePath() + "'");
        out.close();
    }

    /**
     * Print a DOT representation of search graph.
     *
     * @param fileName name of the file where to print the DOT representation to.
     *
     * @throws IllegalArgumentException if {@code fileName} is {@code null}.
     * @throws FileNotFoundException if no file named {@code fileName} cannot be created or written.
     * @throws IllegalStateException if there was some trouble when writing the DOT representation.
     */
    @SuppressWarnings("unused") // Available for debugging purposes.
    public void printDOTFile(final String fileName) throws FileNotFoundException {
        printDOT(new File(fileName));
    }

    /**
     * Get the score of a give sequence of bases
     *
     * @param bases the base sequence.
     *
     * @return {@link Double#NaN} if there is no score for the sequence, i.e. there is no such a haplotype accessible
     *   throw this finder.
     */
    public double score(final byte[] bases) {
        return topFinder.score(bases,0,bases.length);
    }

    /**
     * Get the score of a give sequence of bases
     *
     * @param haplotype the haplotype.
     *
     * @return {@link Double#NaN} if there is no score for the sequence, i.e. there is no such a haplotype accessible
     *   throw this finder.
     */
    public double score(final Haplotype haplotype) {
        return score(haplotype.getBases());
    }


    /**
     * Returns a unique list of haplotypes solutions.
     * <p>
     *     The result will not contain more than one haplotype with the same base sequence. The solution of the best
     *     score is returned.
     * </p>
     * <p>
     *     This makes sense when there are more than one possible path through the graph to create the same haplotype.
     * </p>
     * <p>
     *     The resulting list is sorted by the score with more likely haplotype search results first.
     * </p>
     *
     * @param maxSize maximum number of unique results to return.
     *
     * @throws IllegalArgumentException if {@code maxSize} is negative.
     *
     * @return never {@code null}, perhaps an empty list.
     */
    public List<KBestHaplotype> unique(final int maxSize) {
        if (maxSize < 0) throw new IllegalArgumentException("maxSize cannot be negative");
        final int requiredCapacity = Math.min(maxSize,size());
        final Set<Haplotype> haplotypes = new HashSet<>(requiredCapacity);
        int resultSize = 0;
        final List<KBestHaplotype> result = new ArrayList<>(requiredCapacity);
        for (final KBestHaplotype kbh : this) {
            if (haplotypes.add(kbh.haplotype())) {
                result.add(kbh);
                if (resultSize == maxSize) break;
            }
        }
        return result;
    }

    /**
     * Returns a unique list of haplotypes solutions.
     *
     * <p>
     *     The result will not contain more than one haplotype with the same base sequence. The solution of the best
     *     score is returned.
     * </p>
     * <p>
     *     This makes sense when there are more than one possible path through the graph to create the same haplotype.
     * </p>
     * <p>
     *     The resulting list is sorted by the score with more likely haplotype search results first.
     * </p>
     *
     * @return never {@code null}, perhaps an empty list.
     */
    public List<KBestHaplotype> unique() {
        return unique(size());
    }
}
