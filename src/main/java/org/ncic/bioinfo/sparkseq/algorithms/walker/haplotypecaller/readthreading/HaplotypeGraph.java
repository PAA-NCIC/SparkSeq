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

import org.apache.log4j.Logger;
import org.ncic.bioinfo.sparkseq.algorithms.utils.SequenceComplexity;
import org.ncic.bioinfo.sparkseq.algorithms.utils.Utils;
import org.ncic.bioinfo.sparkseq.algorithms.data.basic.CountSet;
import org.ncic.bioinfo.sparkseq.algorithms.data.basic.Pair;
import org.ncic.bioinfo.sparkseq.algorithms.utils.haplotype.Haplotype;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.graphs.Kmer;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.graphs.MultiSampleEdge;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.graphs.Route;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.graphs.RouteFinder;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

/**
 *
 * Threading graph subclass used to "re-thread" haplotypes instead of reads.
 *
 * Created with IntelliJ IDEA.
 * User: valentin
 * Date: 8/23/13
 * Time: 2:42 PM
 * To change this template use File | Settings | File Templates.
 *
 * Author: wbc
 */
public class HaplotypeGraph extends ReadThreadingGraph {

    /**
     * Maximum repeat unit length considered when looking for repeats that should not be considered as
     * possible read anchor places along the reference path.
     */
    protected static final int DEFAULT_MAX_REPEAT_UNIT_LENGTH = 4;

    /**
     * Minimum repeat length to consider a region a repeat that should not be considered as possibl read anchor
     * places along the reference path.
     */
    protected static final int DEFAULT_MIN_REPEAT_LENGTH_IN_UNITS = 6;

    /**
     * Reference haplotype
     */
    private Haplotype referenceHaplotype;

    /**
     * Reference haplotype bases
     */
    private byte[] referenceBases;

    /**
     * Sets of haplotypes in the graph.
     */
    private Set<Haplotype> haplotypes;

    /**
     * Route of haplotypes in the graph.
     */
    private HaplotypeRoute referenceRoute;

    /**
     * Set of vertices along the reference route.
     */
    private Set<MultiDeBruijnVertex> referenceVertices;

    /**
     * Holds haplotype routes by haplotype.
     */
    private Map<Haplotype,HaplotypeRoute> haplotypeRouteByHaplotype;

    /**
     * Holds haplotypes by contained vertices.
     */
    private Map<MultiDeBruijnVertex,Set<HaplotypeRoute>> haplotypesByVertex;

    /**
     * Reference to the logger for this class.
     */
    private static final Logger logger = Logger.getLogger(HaplotypeGraph.class);

    /**
     * What is the maximum STR unit length.
     */
    private int maxRepeatUnitLength = DEFAULT_MAX_REPEAT_UNIT_LENGTH;

    /**
     * What is the minimum length in units for a STR.
     */
    private int minRepeatLengthInUnits = DEFAULT_MIN_REPEAT_LENGTH_IN_UNITS;


    /**
     * Indicates that the haplotype data structures need update previous to querying.
     */
    private boolean needToUpdateHaplotypeStructures = true;
    private Set<MultiDeBruijnVertex> anchorableVertices;

    /**
     * Constructs a haplotype graph from a describing string.
     *
     * <p>Used for testing</p>
     * @param string the string representation of the haplotype graph.
     */
    public HaplotypeGraph(final String string) {
        super(string);
        haplotypes = new LinkedHashSet<>(10);
        referenceVertices = Collections.emptySet();
    }

    /**
     * Constructs a new haplotype graph given its kmerSize.
     *
     * @param kmerSize 1 or greater, the targeted kmerSize
     *
     * @throws IllegalArgumentException if {@code kmerSize} is 0 or negative.
     */
    public HaplotypeGraph(final int kmerSize) {
        super(kmerSize);
        haplotypes = new LinkedHashSet<>(10);
        referenceVertices = Collections.emptySet();
    }


    /**
     * Set of vertices along the reference haplotype path.
     *
     * @return never {@code} null but perhaps empty.
     */
    public Set<MultiDeBruijnVertex> getReferenceVertices() {
        updateHaplotypeStructures();
        return referenceVertices;
    }

    /**
     * Returns the haplotype route given an haplotype.
     * @param haplotype query haplotype
     * @throws NullPointerException if {@code haplotype} is {@code null}.
     * @throws IllegalArgumentException if {@code haplotype} is not a supported haplotype in the graph.
     * @return never {@code null}.
     */
    public HaplotypeRoute getHaplotypeRoute(final Haplotype haplotype) {
        updateHaplotypeStructures();
        if (!haplotypes.contains(haplotype))
            throw new IllegalArgumentException("input haplotype must be part of the haplotype graph haplotype set");
        HaplotypeRoute result = haplotypeRouteByHaplotype.get(haplotype);
        if (result == null)
            haplotypeRouteByHaplotype.put(haplotype,result = buildHaplotypeRoute(haplotype));
        return result;
    }

    /**
     * Creates an haplotype route.
     * @param haplotype the target haplotype
     * @return {@code null} if there is no such a route in the graph.
     */
    private HaplotypeRoute buildHaplotypeRoute(final Haplotype haplotype) {
        final Route<MultiDeBruijnVertex,MultiSampleEdge> route = RouteFinder.findRoute(this,haplotype.getBases());
        if (route == null)
            return null;
        else
            return new HaplotypeRoute(route);
    }

    /**
     * Bases along the reference path.
     *
     * @return {@code null} if there is no reference.
     */
    @SuppressWarnings("unused")
    public byte[] getReferenceBases() {
        updateHaplotypeStructures();
        return referenceBases;
    }

    /**
     * Returns the reference haplotype
     * @return {@code null} if there is no such a reference.
     */
    public Haplotype getReferenceHaplotype() {
        updateHaplotypeStructures();
        return referenceHaplotype;
    }



    /**
     * Construct a haplotype graph given the haplotype list and the elected kmerSize.
     *
     * @param haplotypes whose path to add to the graph.
     * @param kmerSize   the kmerSize use to compose the graph.
     */
    public HaplotypeGraph(final int kmerSize, final List<Haplotype> haplotypes) {
        super(kmerSize);
        referenceHaplotype = findReferenceHaplotypeOrFail(haplotypes);
        this.haplotypes = new LinkedHashSet<>(haplotypes);
        addSequence("anonymous", referenceHaplotype.getBases(), true);
        for (final Haplotype h : haplotypes) {
            if (h.isReference())
                continue;
            if (h.length() < kmerSize) {
                Utils.warnUser(logger, "haplotype shorter than kmerSize " + h.length() + " < " + kmerSize + " will be dropped");
            } else
                addSequence("anonymous", h.getBases(), false);

        }
        buildGraphIfNecessary();
    }

    /**
     * Returns the reference haplotype within the input collection.
     *
     * @param haplotypes the query haplotype set.
     * @throws IllegalArgumentException if there is no reference haplotype.
     * @throws NullPointerException if {@code haplotypes} is {@code null} or contains some {@code null} value.
     * @return never {@code} null, a haplotype that is reference.
     */
    private Haplotype findReferenceHaplotypeOrFail(final List<Haplotype> haplotypes) {
        for (final Haplotype h : haplotypes)
            if (h.isReference())
                return h;
        throw new IllegalArgumentException("no reference haplotype present");
    }

    /**
     * Constructs a new haplotype graph given a template read-threading graph and set of haplotypes
     *
     * @param template the template read-threading graph.
     * @param haplotypes the haplotype set to consider
     */
    public HaplotypeGraph(final ReadThreadingGraph template, final List<Haplotype> haplotypes) {
        this(template.getKmerSize());
        referenceHaplotype = findReferenceHaplotypeOrFail(haplotypes);
        this.haplotypes = new HashSet<>(haplotypes);
        template.buildGraphIfNecessary();
        uniqueKmers = new HashMap<>();
        nonUniqueKmers = new HashSet<>();
        // Copy vertices over.
        addVertices(template.vertexSet());
        // Copy edges over.
        for (final MultiSampleEdge edge : template.edgeSet()) {
              final MultiSampleEdge newEdge = addEdge(template.getEdgeSource(edge), template.getEdgeTarget(edge));
              newEdge.setIsRef(newEdge.isRef());
              newEdge.setMultiplicity(edge.getMultiplicity());
        }
        // Copy kmer lookup tables:
        uniqueKmers.putAll(template.uniqueKmers);
        nonUniqueKmers.addAll(template.nonUniqueKmers);
        alreadyBuilt = true;
    }

    /**
     * Update the haplotype data structures based in current edges and vertices.
     */
    private void updateHaplotypeStructures() {
        if (!needToUpdateHaplotypeStructures)
            return;
        needToUpdateHaplotypeStructures = false;
        haplotypeRouteByHaplotype = new LinkedHashMap<>(haplotypes.size());
        final Iterator<Haplotype> haplotypeIterator = haplotypes.iterator();
        final Set<Haplotype> nonFoundHaplotypes = new HashSet<>(haplotypes.size());
        while (haplotypeIterator.hasNext()) {
            final Haplotype haplotype = haplotypeIterator.next();
            final HaplotypeRoute haplotypeRoute = buildHaplotypeRoute(haplotype);
            if (haplotypeRoute == null) {
                haplotypeIterator.remove();
                nonFoundHaplotypes.add(haplotype);
                if (haplotype.isReference()) {
                    referenceHaplotype = null;
                    referenceRoute = null;
                    referenceVertices = Collections.emptySet();
                    referenceBases = null;
                }
            } else {
                if (haplotype.isReference()) {
                    referenceHaplotype = haplotype;
                    referenceRoute = haplotypeRoute;
                    referenceVertices = haplotypeRoute.vertexSet();
                    referenceBases = haplotypeRoute.getBases();
                }
                haplotypeRouteByHaplotype.put(haplotype, haplotypeRoute);
            }
        }
        haplotypesByVertex = buildHaplotypesByVertex();
        anchorableVertices = calculateAnchorableVertexSet();
        logger.debug("some haplotypes do not have a path across the haplotype graph " + nonFoundHaplotypes.size());
    }

    /**
     * Builds a map for each vertex to all the haplotype routes that pass thru it.
     */
    private Map<MultiDeBruijnVertex, Set<HaplotypeRoute>> buildHaplotypesByVertex() {
        final Map<MultiDeBruijnVertex, Set<HaplotypeRoute>> result = new HashMap<>(referenceVertices.size());
        final Set<HaplotypeRoute> allHaplotypeRoutes = new LinkedHashSet<>(haplotypeRouteByHaplotype.values());
        for (final HaplotypeRoute haplotypeRoute : allHaplotypeRoutes) {
            final Set<HaplotypeRoute> singleton = Collections.singleton(haplotypeRoute);
            for (final MultiDeBruijnVertex vertex : haplotypeRoute.vertexSet())
                if (!result.containsKey(vertex))
                    result.put(vertex, singleton);
                else {
                    final Set<HaplotypeRoute> currentHrs = result.get(vertex);
                    if (currentHrs.size() == haplotypes.size() - 1)
                        result.put(vertex, allHaplotypeRoutes);
                    else if (currentHrs.size() == 1) {
                        final Set<HaplotypeRoute> newHrs = new LinkedHashSet<>(allHaplotypeRoutes.size());
                        newHrs.addAll(currentHrs);
                        newHrs.add(haplotypeRoute);
                        result.put(vertex, newHrs);
                    } else
                        currentHrs.add(haplotypeRoute);
                }
        }
        return result;
    }


    /**
     * Debug convenient method to print a graph into a file in the .dot format.
     * @param fileName name of the output file.
     * @throws NullPointerException if {@code fileName} is {@code null}.
     */
    public void printGraph(final String fileName) {
        super.printGraph(new File(fileName), 10000);
    }




    @Override
    public void printGraph(final PrintStream graphWriter, final boolean writeHeader, final int pruneFactor) {
        if ( writeHeader )
            graphWriter.println("digraph assemblyGraphs {");


        for( final MultiSampleEdge edge : edgeSet() ) {
            graphWriter.println("\t" + getEdgeSource(edge).toString() + " -> " + getEdgeTarget(edge).toString() + " [" + (edge.getMultiplicity() > 0 && edge.getMultiplicity() <= pruneFactor ? "style=dotted,color=grey," : "") + "label=\"" + edge.getDotLabel() + "\"];");
            if( edge.isRef() ) {
                graphWriter.println("\t" + getEdgeSource(edge).toString() + " -> " + getEdgeTarget(edge).toString() + " [color=red];");
            }
        }

        for( final MultiDeBruijnVertex v : vertexSet() )
            graphWriter.println("\t" + v.toString() + " [label=\"" + v.getId() + ":" +  new String(getAdditionalSequence(v)) + v.additionalInfo() + "\",shape=box]");

        if ( writeHeader )
            graphWriter.println("}");
    }

    @Override
    protected int findStart(final SequenceForKmers seqForKmers) {
        return 0;
    }

    /**
     * Checks whether the graph has some sources or sink vertices that are not reference vertices.
     *
     * @return {@code true} iff so.
     */
    public boolean hasNonReferenceEnds() {
        for (final MultiDeBruijnVertex end : getSources())
            if (!isReferenceNode(end)) return true;
        for (final MultiDeBruijnVertex end : getSinks())
            if (!isReferenceNode(end)) return true;
        return false;
    }

    /**
     * Merges vertices that share exactly the same set of outgoing vertices.
     * <p/>
     * This is done in reversed topological order and since the graph is a DAG it ensure to return a graph
     * that such merge is any longer possible. I.e. there is no need to run this method more than once.
     * <p/>
     * Notice that we will a record of distinct unique kmers that map to the same vertex that map now to the same
     * merged vertex. Thus if vertices {@code X and Y} are merged then {@code findKmer(X.sequence) == findKmer(Y.sequence)}.
     * <p/>
     * Examples:
     * <ul>
     * <li>
     *    {@code AAA -> AAC, CAA -> AAC} would become {@code NAA -> AAC}.
     * </li><li>
     *    {@code AAA -> AAC, AAA -> AAG, CAA -> AAC, CAA -> AAG} would become {@code NAA -> AAG, NAA -> AAG}
     * </li><li>
     *    {@code AAA -> AAC, AAA -> AAG, CAA -> AAC} would not change as {@code AAA} and {@code CAA}
     *    do not share {@code AAG} as outgoing vertex.
     * </li>
     * <li>
     *    {@code AAA -> AAC, AAC -> ACA, CAA -> AAC, GAC -> ACA } would become {@code NAA -> NAC, NAC -> ACA}.
     * </li>
     * </ul>
     */
    public void mergeCommonChains() {
        final int vertexCount = vertexSet().size();
        final Set<MultiDeBruijnVertex> refVertices = new HashSet<>(vertexCount);
        final Map<MultiDeBruijnVertex,Integer> indexByVertex = new HashMap<>(vertexCount);
        final int[] pendingChildren = new int[vertexCount];
        final Deque<MultiDeBruijnVertex> readyVertices = new LinkedList<>();
        final Set<MultiDeBruijnVertex> merged = new HashSet<>(1 + vertexCount / 10 );

        // Initialize traversal data structures.
        mergeCommonChainsInitialize(refVertices, indexByVertex, pendingChildren, readyVertices);

        // Traversal in inverted topological order where children nodes are processed before their parents.
        while (!readyVertices.isEmpty()) {
            final MultiDeBruijnVertex currentVertex = readyVertices.remove();
            if (merged.contains(currentVertex)) continue;

            final Set<MultiDeBruijnVertex> mergeSet = new HashSet<>(2);
            MultiDeBruijnVertex refVertex = mergeCommonChainsComposeMergeSet(refVertices, currentVertex, mergeSet);
            mergeVertices(refVertex,mergeSet,indexByVertex,pendingChildren,readyVertices);
            merged.addAll(mergeSet);
        }
        needToUpdateHaplotypeStructures = true;
    }

    /**
     * Given a seed vertex, determines the mergin set of nodes that will be collapsed into one.
     *
     * @param refVertices reference path vertices
     * @param currentVertex current vertex.
     * @param mergeSet where to store the final merging set.
     * @return the reference node if present that needs to be preserved as such. It might be {@code null}
     */
    private MultiDeBruijnVertex mergeCommonChainsComposeMergeSet(final Set<MultiDeBruijnVertex> refVertices,
                                                                 final MultiDeBruijnVertex currentVertex,
                                                                 final Set<MultiDeBruijnVertex> mergeSet) {
        final boolean currentIsSource = isSource(currentVertex);
        final Set<MultiDeBruijnVertex> children = outgoingVerticesOf(currentVertex);
        if (children.size() == 0)
            mergeSet.add(currentVertex);
        else
            for (final MultiDeBruijnVertex child : children)
                mergeSet.addAll(incomingVerticesOf(child));

        MultiDeBruijnVertex refVertex = refVertices.contains(currentVertex) ? currentVertex : null;
        final Iterator<MultiDeBruijnVertex> candidatesIt = mergeSet.iterator();
        while (candidatesIt.hasNext()) {
            final MultiDeBruijnVertex candidate = candidatesIt.next();
            if (candidate == currentVertex) continue;
            if (isSource(candidate) != currentIsSource) {
                candidatesIt.remove();
                continue;
            }
            if (currentIsSource && !candidate.getSequenceString().equals(currentVertex.getSequenceString())) {
                candidatesIt.remove();
                continue;
            }
            if (!currentIsSource && candidate.getSuffix() != currentVertex.getSuffix()) {
                candidatesIt.remove();
                continue;
            }
            final Set<MultiDeBruijnVertex> candidateChildren = outgoingVerticesOf(candidate);
            if (candidateChildren.size() != children.size())
                candidatesIt.remove();
            else {
                boolean removed = false;
                for (final MultiDeBruijnVertex candidateChild : candidateChildren)
                    if (!children.contains(candidateChild)) {
                        candidatesIt.remove();
                        removed = true;
                        break;
                    }
                if (refVertex == null && !removed && refVertices.contains(candidate)) refVertex = candidate;
            }
        }
        return refVertex;
    }

    /**
     * Initialize data-structures for {@link #mergeCommonChains}
     *
     * @param refVertices will contain reference path vertices.
     * @param indexByVertex map vertex -> index in {@code pendingChildren}.
     * @param pendingChildren number of children of a node that have not yet been processed.
     * @param readyVertices vertices that are ready to be processed (all children have been processed).
     */
    private void mergeCommonChainsInitialize(final Set<MultiDeBruijnVertex> refVertices,
                                             final Map<MultiDeBruijnVertex, Integer> indexByVertex,
                                             final int[] pendingChildren,
                                             final Deque<MultiDeBruijnVertex> readyVertices) {
        int nextIndex = 0;
        for (final MultiDeBruijnVertex v : vertexSet()) {
            indexByVertex.put(v,nextIndex++);
            if (isReferenceNode(v)) refVertices.add(v);
        }

        for (final Map.Entry<MultiDeBruijnVertex,Integer> entry : indexByVertex.entrySet())
            if ((pendingChildren[entry.getValue()] = outDegreeOf(entry.getKey())) == 0)
                readyVertices.add(entry.getKey());
    }

    // Perform the actual merge.
    private void mergeVertices(final MultiDeBruijnVertex refVertex, final Collection<MultiDeBruijnVertex> vertices, final Map<MultiDeBruijnVertex,Integer> indexByVertex, final int[] pendingChildrenCounts, final Deque<MultiDeBruijnVertex> ready) {
        if (vertices.size() == 0)
            throw new IllegalArgumentException();
        final MultiDeBruijnVertex vertexToKeep =  refVertex == null ? vertices.iterator().next() : refVertex;
        final byte[] sequence = vertexToKeep.getSequence();
        final Set<Kmer> uniqueKmersToUpdate = new HashSet<>(vertices.size());
        final Set<MultiDeBruijnVertex> parentVertices = new HashSet<>(inDegreeOf(vertexToKeep) * 2);
        parentVertices.addAll(incomingVerticesOf(vertexToKeep));
        for (final MultiDeBruijnVertex p : parentVertices)
            if (--pendingChildrenCounts[indexByVertex.get(p)] == 0)
                ready.add(p);

        final Kmer mergedKmer = new Kmer(sequence);
        if (uniqueKmers.containsKey(mergedKmer)) {
            uniqueKmersToUpdate.add(new Kmer(mergedKmer.bases().clone()));
            uniqueKmers.remove(mergedKmer);
        }
        boolean foundMergedVertex = false;
        for (final MultiDeBruijnVertex v : vertices)
            if (v == vertexToKeep)
                foundMergedVertex = true;
            else {
                final byte[] seq = v.getSequence();
                final Kmer kmer = new Kmer(seq);
                if (uniqueKmers.containsKey(kmer)) {
                    uniqueKmersToUpdate.add(kmer);
                    uniqueKmers.remove(kmer);
                }
                if (sequence.length != seq.length) throw new IllegalArgumentException("mismatched sizes " + sequence.length + " != "
                        + seq.length + " " + new String(sequence) + " " + new String(seq));
                for (int i = sequence.length - 1; i >= 0; i--) {

                    if (sequence[i] != seq[i]) sequence[i] = 'N';
                }
                for (final MultiDeBruijnVertex p : incomingVerticesOf(v)) {
                    if (--pendingChildrenCounts[indexByVertex.get(p)] == 0)
                        ready.add(p);
                    if (!parentVertices.contains(p)) {
                        parentVertices.add(p);
                        final MultiSampleEdge e = getEdge(p,v);
                        addEdge(p,vertexToKeep,new MultiSampleEdge(e.isRef(),e.getMultiplicity(),1));
                    } else {
                        getEdge(p,vertexToKeep).incMultiplicity(getEdge(p,v).getMultiplicity());
                    }
                }
                removeVertex(v);
            }
        if (!foundMergedVertex)
            throw new IllegalArgumentException("merged vertex must be contained in the input set");
        for (final Kmer kmer : uniqueKmersToUpdate)
            uniqueKmers.put(kmer,vertexToKeep);
    }

    public Map<Kmer,MultiDeBruijnVertex> uniqueKmerMap() {
        return Collections.unmodifiableMap(uniqueKmers);
    }

    @Override
    public boolean equals(Object other) {
        return (other instanceof HaplotypeGraph) && equals((HaplotypeGraph)other);
    }


    /**
     * Simple debug representation of the haplotype graph.
     * @return never {@code null}
     */
    @Override
    public String toString() {
        return getClass().getSimpleName() + "[ks=" + kmerSize + "](vs=" + vertexSet().size() + "," + edgeSet().size() + "){...}";
    }

    /**
     * Returns set of valid haplotypes.
     * @return never {@code null} but perhaps empty.
     */
    public Set<Haplotype> getHaplotypes() {
        updateHaplotypeStructures();
        return haplotypes;
    }

    /**
     * Returns a map between valid haplotypes and corresponding routes in the graph.
     * @return never {@code null} but perhaps empty.
     */
    public Map<Haplotype,HaplotypeRoute> getHaplotypeRouteMap() {
        updateHaplotypeStructures();
        return haplotypeRouteByHaplotype;
    }

    /**
     * Returns set of haplotype routes that enclose a vertex.
     * @param vertex the query vertex.
     * @return never {@code null} but perhaps empty set.
     */
    public Set<HaplotypeRoute> getEnclosingHaplotypeRoutes(final MultiDeBruijnVertex vertex) {
        updateHaplotypeStructures();
        if (haplotypesByVertex == null)
            return Collections.emptySet();
        final Set<HaplotypeRoute> result = haplotypesByVertex.get(vertex);
        if (result == null)
            return Collections.emptySet();
        else
            return result;
    }

    /**
     * Returns the reference route
     *
     * @return {@code null} if there is no valid reference haplotype.
     */
    public HaplotypeRoute getReferenceRoute() {
        updateHaplotypeStructures();
        return referenceRoute;
    }

    /***********************************************
     * deep equals implementation, used in testing. *
     ***********************************************/

    /**
     * Compare two haplotype threading graphs and it determines whether they have the same structure.
     * <p/>
     * This method goes a long way to figure out the equality and no equality of both graphs. However there
     * are "pathological" case in where it might fail to see a difference. This is due to the fact that there
     * is no guarantee of the uniqueness of sequences at source vertex.
     * <p/>
     * If there are more than one source vertex with the same sequence it try to match source vertices between both
     * graphs matching all possible paths emanating from every pair of sources.
     *
     * <p>Note: in practice this is only used in for testing purposes</p.>
     *
     * @param other  the other graph to compare against.
     * @return never {@code null}.
     */
    public boolean equals(HaplotypeGraph other) {
        updateHaplotypeStructures();
        if (other == null) return false;
        if (other == this) return true;

        if (!equals$ReferencePaths(this, other)) return false;
        final Map<String,MultiDeBruijnVertex> thisSourcesBySequence = equalsBuildSourceBySequenceMap(this);
        final Map<String,MultiDeBruijnVertex> otherSourcesBySequence = equalsBuildSourceBySequenceMap(other);
        if (thisSourcesBySequence.size() != otherSourcesBySequence.size()) return false;
        final List<MultiDeBruijnVertex> unmatchedLeft = new LinkedList<>();
        final List<MultiDeBruijnVertex> unmatchedRight = new LinkedList<>();
        final List<Pair<MultiDeBruijnVertex,MultiDeBruijnVertex>> sourcePairs = equals$matchVertexBySequenceMaps(thisSourcesBySequence,otherSourcesBySequence,unmatchedLeft,unmatchedRight);
        if (unmatchedLeft.size() > 0 || unmatchedRight.size() > 0) return false;


        final Deque<Pair<MultiDeBruijnVertex,MultiDeBruijnVertex>> pending = new LinkedList<>(sourcePairs);
        final Set<MultiDeBruijnVertex> visited = new HashSet<>(vertexSet().size());
        while (!pending.isEmpty()) {
            final Pair<MultiDeBruijnVertex,MultiDeBruijnVertex> pair = pending.removeFirst();
            final MultiDeBruijnVertex leftVertex = pair.getFirst();
            final MultiDeBruijnVertex rightVertex = pair.getSecond();
            final List<Pair<MultiDeBruijnVertex,MultiDeBruijnVertex>> childrenPairs = equals$matchVertexBySequenceMaps(equalsBuildChildrenBySuffixMap(this, leftVertex),
                    equalsBuildChildrenBySuffixMap(other, rightVertex), unmatchedLeft, unmatchedRight);
            if (unmatchedLeft.size() > 0 || unmatchedRight.size() > 0) return false;
            for (final Pair<MultiDeBruijnVertex,MultiDeBruijnVertex> childPair : childrenPairs) {
                final MultiDeBruijnVertex leftChild = childPair.getFirst();
                final MultiDeBruijnVertex rightChild = childPair.getSecond();
                final boolean leftVisited = visited.add(leftChild);
                final boolean rightVisited = visited.add(rightChild);
                if (leftVisited != rightVisited) return false; // visited before in different matchings.
                if (leftVisited) continue;
                pending.add(childPair);
                visited.add(childPair.getFirst());
                visited.add(childPair.getSecond());
            }
        }
        return true;
    }

    // Note: in practice only use by equals(HaplotypeGraph) for testing purposes.
    private boolean equals$ReferencePaths(final HaplotypeGraph g1, final HaplotypeGraph g2) {
        MultiDeBruijnVertex refVertex1 = g1.getReferenceSourceVertex();
        MultiDeBruijnVertex refVertex2 = g2.getReferenceSourceVertex();
        if (refVertex1 == null && refVertex2 == null)
            return true;
        if (refVertex1 == null || refVertex2 == null)
            return false;

        if (!refVertex1.getSequenceString().equals(refVertex2.getSequenceString()))
            return false;

        while (refVertex1 != null && refVertex2 != null) {
            if (refVertex1.getSuffix() != refVertex2.getSuffix()) return false;
            refVertex1 = g1.getNextReferenceVertex(refVertex1);
            refVertex2 = g2.getNextReferenceVertex(refVertex2);

        }
        return refVertex1 == refVertex2;
    }

    // Note: in practice only use by equals(HaplotypeGraph) for testing purposes.
    private static Map<String, MultiDeBruijnVertex> equalsBuildChildrenBySuffixMap(final HaplotypeGraph graph,
                                                                                   final MultiDeBruijnVertex vertex) {
        final Map<String,MultiDeBruijnVertex> result = new HashMap<>();
        for (final MultiDeBruijnVertex child : graph.outgoingVerticesOf(vertex))
            result.put(new String(new byte[]{child.getSuffix()}), child);
        return result;
    }

    // Note: in practice only use by equals(HaplotypeGraph) for testing purposes.
    private static List<Pair<MultiDeBruijnVertex, MultiDeBruijnVertex>> equals$matchVertexBySequenceMaps(
            final Map<String,MultiDeBruijnVertex> left, final Map<String, MultiDeBruijnVertex> right,
            final Collection<MultiDeBruijnVertex> unmatchedLeft, final Collection<MultiDeBruijnVertex> unmatchedRight) {
        final List<Pair<MultiDeBruijnVertex,MultiDeBruijnVertex>> result = new LinkedList<>();
        for (final Map.Entry<String,MultiDeBruijnVertex> leftEntry : left.entrySet())
            if (right.containsKey(leftEntry.getKey()))
                result.add(new Pair<>(leftEntry.getValue(),right.get(leftEntry.getKey())));
            else
                unmatchedLeft.add(leftEntry.getValue());
        for (final Map.Entry<String,MultiDeBruijnVertex> rightEntry : right.entrySet())
            if (!left.containsKey(rightEntry.getKey()))
                unmatchedRight.add(rightEntry.getValue());
        return result;
    }

    // Note: in practice only use by equals(HaplotypeGraph) for testing purposes.
    private static Map<String, MultiDeBruijnVertex> equalsBuildSourceBySequenceMap(final HaplotypeGraph other) {

        final Set<MultiDeBruijnVertex> sources = other.getSources();
        final Map<String,MultiDeBruijnVertex> result = new HashMap<>(sources.size());
        final Map<String,List<MultiDeBruijnVertex>> collisions = new HashMap<>(sources.size());
        for (final MultiDeBruijnVertex v : sources) {
            final String sequence = v.getSequenceString();
            if (result.containsKey(sequence)) { // we need to handle collision due to lack of uniqueness.
                final List<MultiDeBruijnVertex> collisionList;
                if (collisions.containsKey(sequence))
                    collisionList = collisions.get(sequence);
                else
                    collisions.put(sequence,collisionList = new LinkedList<>());
                collisionList.add(v);
            } else {
                result.put(sequence,v);
            }
        }
        if (collisions.size() == 0)
            return result;
        for (final String s : collisions.keySet()) {
            result.remove(s);
            final List<MultiDeBruijnVertex> vertices = collisions.remove(s);
            int number = 0;
            final List<Pair<MultiDeBruijnVertex,String>> extendedSequences = new LinkedList<>();
            for (final MultiDeBruijnVertex vertice : vertices)
                extendedSequences.add(new Pair<>(vertice, equalsCollisionResolverExtendedSequence(other, vertice)));
            Collections.sort(extendedSequences,new Comparator<Pair<MultiDeBruijnVertex,String>>(){
                public int compare(final Pair<MultiDeBruijnVertex,String> p1, final Pair<MultiDeBruijnVertex,String> p2) {
                    return p1.getSecond().compareTo(p2.getSecond());
                }
            });
            for (final Pair<MultiDeBruijnVertex,String> p : extendedSequences)
                result.put(p.getSecond() + '-' + (number++),p.getFirst());
        }
        return result;

    }

    // Note: in practice only use by equals(HaplotypeGraph) for testing purposes.
    private static String equalsCollisionResolverExtendedSequence(final HaplotypeGraph graph, final MultiDeBruijnVertex source) {
        final StringBuilder buffer = new StringBuilder(1000);
        final Set<MultiDeBruijnVertex> visited = new HashSet<>(graph.vertexSet().size());
        final Stack<MultiDeBruijnVertex> pending = new Stack<>();
        final Stack<Integer> position = new Stack<>();
        position.ensureCapacity(graph.vertexSet().size());
        pending.ensureCapacity(graph.vertexSet().size());
        pending.add(source);
        position.add(0);
        int lastPos = -1;
        while (!pending.isEmpty()) {
            final MultiDeBruijnVertex next = pending.pop();
            if (visited.contains(next)) continue;
            visited.add(next);
            final int pos = position.pop();
            final CharSequence sequence;
            if (graph.isSource(next)) {
                if (next == source) {
                    sequence = new String(next.getSequence());
                } else {
                    sequence = new StringBuffer(next.getSequence().length).append(new String(next.getSequence())).reverse().append('$');
                }
            } else {
                sequence = new String(new byte[] { next.getSuffix()});
            }

            if (pos != lastPos + 1) {
                buffer.append('[').append(Math.abs(pos)).append(']');
            }
            buffer.append(sequence);
            lastPos = pos + sequence.length() - 1;

            final List<MultiDeBruijnVertex> parents = new LinkedList<>(graph.incomingVerticesOf(next));
            Collections.sort(parents,new Comparator<MultiDeBruijnVertex>() {
                @Override
                public int compare(final MultiDeBruijnVertex o1, final MultiDeBruijnVertex o2) {
                    return Byte.compare(o1.getSuffix(),o2.getSuffix());
                }
            });
            for (final MultiDeBruijnVertex parent : parents) {
                pending.push(parent);
                position.push(lastPos + 1);
            }

            final List<MultiDeBruijnVertex> children = new LinkedList<>(graph.incomingVerticesOf(next));
            Collections.sort(children,new Comparator<MultiDeBruijnVertex>() {
                @Override
                public int compare(final MultiDeBruijnVertex o1, final MultiDeBruijnVertex o2) {
                    return Byte.compare(o1.getSuffix(),o2.getSuffix());
                }
            });
            for (final MultiDeBruijnVertex child : graph.outgoingVerticesOf(next)) {
                pending.push(child);
                position.push(lastPos + 1);
            }
        }

        return buffer.toString();
    }


    /**
     * Calculates the subset of reference path vertices that are amenable to be anchoring vertices.
     * <p/>
     * <p>
     * For a vertex to be anchorable:
     * <ul>
     * <li>Should not include bases from a repeat</li>,
     * <li>There should not be in a middle of a event block</li>
     * </ul>
     * </p>
     *
     * @return never {@code null}.
     */
    private Set<MultiDeBruijnVertex> calculateAnchorableVertexSet() {
        updateHaplotypeStructures();
        if (referenceBases == null)
            return Collections.emptySet();

        // We first check what bases in the reference path bases are part of a repeat.
        final boolean[] nonAnchorableDueToRepeats = SequenceComplexity.findBasesInShortUnitRepeats(
                referenceBases, maxRepeatUnitLength, minRepeatLengthInUnits);

        final Set<MultiDeBruijnVertex> result = new HashSet<>(100);
        final Map<MultiDeBruijnVertex, CountSet> expectedRejoins = new HashMap<>();


        MultiDeBruijnVertex currentVertex = getReferenceRoute().getFirstVertex();
        final int sourceSequenceLength = currentVertex.getSequence().length;

        // Determine whether the reference source vertex in anchorable discarding repeats:
        boolean sourceIsAnchorable = true;
        for (int i = 0; i < sourceSequenceLength; i++)
            if (nonAnchorableDueToRepeats[i]) {
                sourceIsAnchorable = false;
                break;
            }

        // Update the nonAnchorableDueToRepeats array accordingly.
        int index = currentVertex.getSequence().length - 1;
        nonAnchorableDueToRepeats[index] = !sourceIsAnchorable;


        // We keep record on all alternative path lengths:
        final CountSet pathLengths = new CountSet(haplotypes.size());
        pathLengths.setTo(0);

        // Now we go through the reference path and determine which vertices are not part of event block.
        // We keep track of open divergent paths in expectedRejoins. Thus only those vertices traversed
        // when exptectedRejoins size 0 can be anchorable:
        while (currentVertex != null) {
            int inDegree = inDegreeOf(currentVertex);
            if (inDegree > 1)
                expectedRejoins.remove(currentVertex);
            if (expectedRejoins.size() == 0 && !nonAnchorableDueToRepeats[index]) {
                currentVertex.setAdditionalInfo(currentVertex.additionalInfo() + "*");
                result.add(currentVertex);
            }
            final Set<MultiSampleEdge> nextEdges = outgoingEdgesOf(currentVertex);
            MultiDeBruijnVertex nextReferenceVertex = null;
            for (final MultiSampleEdge e : nextEdges) {
                final MultiDeBruijnVertex nextVertex = getEdgeTarget(e);
                if (e.isRef() && referenceVertices.contains(nextVertex))
                    nextReferenceVertex = nextVertex;
                else
                    calculateRejoins(nextVertex, expectedRejoins, referenceVertices, pathLengths, false, false);
            }
            currentVertex = nextReferenceVertex;
            index++;
        }
        return result;
    }



    /**
     * Returns those vertices that can be used as anchors along the refererence route.
     * @return never {@code null} but perhaps empty if there is no such a vertex.
     */
    public Set<MultiDeBruijnVertex> getAnchorableVertices() {
        updateHaplotypeStructures();
        return anchorableVertices;
    }

    /**
     * Finds non-reference wondering paths that will rejoin the reference path from a particular node.
     * <p/>
     * <p>
     * It only considers those paths that rejoin within the anchor points of a read.
     * </p>
     * <p/>
     * <p>
     * Rather than reporting explicitly the path vertice sequence, this method report the length of the paths
     * found. These are dumped into {@code expectedRejoins} where the keys are refernce path vertex where paths rejoin
     * and the value is the set of path lengths.
     * </p>
     * <p/>
     * <p> The path lengths are calculated as the length from the startVertex plus the prefix sizes {@code prefixSizes} </p>
     * <p/>
     * <p> You can also ask the method to exhaustively find all paths ({@code exhaustive == true}) or just consider
     * intermediate nodes once ({@code exhustive == false}). If the latter only the shortest paths are considered.</p>
     * <p/>
     * <p> Finally you also can check on paths backwards ({@code backwards == true}) or forwards ({@code backwards == false})</p>
     *
     * @param startVertex               the origin node for those paths.
     * @param expectedRejoins           map where to place the found paths in a form of the rejoining non-reference vertex (key) and
     *                                  set of path lengths (value).
     * @param referenceWithinBoundaries reference vertices found between read anchors. The key are the vertices, the values are
     *                                  the kmer's offset in the read.
     * @param prefixSizes               prefix path sizes to be added to the rejoin path sizes.
     * @param exhaustive                whether all paths should be considered or we only care about find out the rejoining vertices.
     * @param backwards                 whether we want to find backward paths (inverse edge traversal).
     *
     * Note: it is marked as deprecated as this method signature may change in the future. It is public just because
     *                                  is currently shared by several other classes, however it would not be surprising if
     *                                  it gets refactored out at some point. So use with care.
     */
    @Deprecated
    public void calculateRejoins(final MultiDeBruijnVertex startVertex, final Map<MultiDeBruijnVertex, CountSet> expectedRejoins,
                                  final Set<MultiDeBruijnVertex> referenceWithinBoundaries, final CountSet prefixSizes,
                                  final boolean exhaustive, final boolean backwards) {
        Queue<MultiDeBruijnVertex> queue = new LinkedList<>();
        Queue<CountSet> depths = new LinkedList<>();
        queue.add(startVertex);
        depths.add(prefixSizes);

        final Set<MultiDeBruijnVertex> visited = new HashSet<>();
        if (!exhaustive) visited.add(startVertex);
        while (!queue.isEmpty()) {
            final CountSet depth = depths.remove();
            final MultiDeBruijnVertex v = queue.remove();
            if (referenceVertices.contains(v)) {
                if (referenceWithinBoundaries.contains(v)) {
                    final CountSet previous = expectedRejoins.get(v);
                    if (previous == null)
                        expectedRejoins.put(v, depth.clone());
                    else
                        previous.addAll(depth);
                }
            } else {
                final CountSet depthPlusOne = depth.clone();
                depthPlusOne.incAll(1);
                final Set<MultiSampleEdge> nextEdges = backwards ? incomingEdgesOf(v) : outgoingEdgesOf(v);
                for (final MultiSampleEdge e : nextEdges) {
                    final MultiDeBruijnVertex w = backwards ? getEdgeSource(e) : getEdgeTarget(e);
                    if (visited.contains(w)) // avoid repetitive work.
                        continue;
                    if (!exhaustive) visited.add(w);
                    queue.add(w);
                    depths.add(depthPlusOne);
                }
            }
        }
    }

}


