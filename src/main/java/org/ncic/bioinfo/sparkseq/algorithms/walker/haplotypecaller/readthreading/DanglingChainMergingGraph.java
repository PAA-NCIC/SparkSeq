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

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.jgrapht.EdgeFactory;
import org.ncic.bioinfo.sparkseq.algorithms.utils.AlignmentUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.smithwaterman.SWPairwiseAlignment;
import org.ncic.bioinfo.sparkseq.algorithms.utils.smithwaterman.SWParameterSet;
import org.ncic.bioinfo.sparkseq.algorithms.utils.smithwaterman.SmithWaterman;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.graphs.BaseGraph;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.graphs.GraphUtils;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.graphs.MultiSampleEdge;

import java.util.*;

public abstract class DanglingChainMergingGraph extends BaseGraph<MultiDeBruijnVertex, MultiSampleEdge> {

    private static final int MAX_CIGAR_COMPLEXITY = 3;
    private int maxMismatchesInDanglingHead = -1;

    protected boolean alreadyBuilt;

    /**
     * Create a new ReadThreadingAssembler using kmerSize for matching
     * @param kmerSize must be >= 1
     */
    protected DanglingChainMergingGraph(final int kmerSize, final EdgeFactory<MultiDeBruijnVertex, MultiSampleEdge> edgeFactory) {
        super(kmerSize, edgeFactory);
    }

    protected void setMaxMismatchesInDanglingHead(final int maxMismatchesInDanglingHead) {
        this.maxMismatchesInDanglingHead = maxMismatchesInDanglingHead;
    }

    /**
     * Edge factory that encapsulates the numPruningSamples assembly parameter
     */
    protected static class MyEdgeFactory implements EdgeFactory<MultiDeBruijnVertex, MultiSampleEdge> {
        final int numPruningSamples;

        public MyEdgeFactory(int numPruningSamples) {
            this.numPruningSamples = numPruningSamples;
        }

        @Override
        public MultiSampleEdge createEdge(final MultiDeBruijnVertex sourceVertex, final MultiDeBruijnVertex targetVertex) {
            return new MultiSampleEdge(false, 1, numPruningSamples);
        }

        public MultiSampleEdge createEdge(final boolean isRef, final int multiplicity) {
            return new MultiSampleEdge(isRef, multiplicity, numPruningSamples);
        }

    }

    /**
     * Class to keep track of the important dangling chain merging data
     */
    protected static final class DanglingChainMergeHelper {
        final List<MultiDeBruijnVertex> danglingPath, referencePath;
        final byte[] danglingPathString, referencePathString;
        final Cigar cigar;

        public DanglingChainMergeHelper(final List<MultiDeBruijnVertex> danglingPath,
                                        final List<MultiDeBruijnVertex> referencePath,
                                        final byte[] danglingPathString,
                                        final byte[] referencePathString,
                                        final Cigar cigar) {
            this.danglingPath = danglingPath;
            this.referencePath = referencePath;
            this.danglingPathString = danglingPathString;
            this.referencePathString = referencePathString;
            this.cigar = cigar;
        }
    }

    /**
     * Try to recover dangling tails
     *
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
     */
    public void recoverDanglingTails(final int pruneFactor, final int minDanglingBranchLength) {
        if ( ! alreadyBuilt )  throw new IllegalStateException("recoverDanglingTails requires the graph be already built");

        int attempted = 0;
        int nRecovered = 0;
        for ( final MultiDeBruijnVertex v : vertexSet() ) {
            if ( outDegreeOf(v) == 0 && ! isRefSink(v) ) {
                attempted++;
                nRecovered += recoverDanglingTail(v, pruneFactor, minDanglingBranchLength);
            }
        }

        logger.debug("Recovered " + nRecovered + " of " + attempted + " dangling tails");
    }

    /**
     * Try to recover dangling heads
     *
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
     */
    public void recoverDanglingHeads(final int pruneFactor, final int minDanglingBranchLength) {
        if ( ! alreadyBuilt )  throw new IllegalStateException("recoverDanglingHeads requires the graph be already built");

        // we need to build a list of dangling heads because that process can modify the graph (and otherwise generate
        // a ConcurrentModificationException if we do it while iterating over the vertexes)
        final List<MultiDeBruijnVertex> danglingHeads = new ArrayList<>();

        int attempted = 0;
        int nRecovered = 0;
        for ( final MultiDeBruijnVertex v : vertexSet() ) {
            if ( inDegreeOf(v) == 0 && ! isRefSource(v) )
                danglingHeads.add(v);
        }

        // now we can try to recover the dangling heads
        for ( final MultiDeBruijnVertex v : danglingHeads ) {
            attempted++;
            nRecovered += recoverDanglingHead(v, pruneFactor, minDanglingBranchLength);
        }

        logger.debug("Recovered " + nRecovered + " of " + attempted + " dangling heads");
    }

    /**
     * Attempt to attach vertex with out-degree == 0 to the graph
     *
     * @param vertex the vertex to recover
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
     * @return 1 if we successfully recovered the vertex and 0 otherwise
     */
    protected int recoverDanglingTail(final MultiDeBruijnVertex vertex, final int pruneFactor, final int minDanglingBranchLength) {
        if ( outDegreeOf(vertex) != 0 ) throw new IllegalStateException("Attempting to recover a dangling tail for " + vertex + " but it has out-degree > 0");

        // generate the CIGAR string from Smith-Waterman between the dangling tail and reference paths
        final DanglingChainMergeHelper danglingTailMergeResult = generateCigarAgainstDownwardsReferencePath(vertex, pruneFactor, minDanglingBranchLength);

        // if the CIGAR is too complex (or couldn't be computed) then we do not allow the merge into the reference path
        if ( danglingTailMergeResult == null || ! cigarIsOkayToMerge(danglingTailMergeResult.cigar, false, true) )
            return 0;

        // merge
        return mergeDanglingTail(danglingTailMergeResult);
    }

    /**
     * Attempt to attach vertex with in-degree == 0, or a vertex on its path, to the graph
     *
     * @param vertex the vertex to recover
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
     * @return 1 if we successfully recovered a vertex and 0 otherwise
     */
    protected int recoverDanglingHead(final MultiDeBruijnVertex vertex, final int pruneFactor, final int minDanglingBranchLength) {
        if ( inDegreeOf(vertex) != 0 ) throw new IllegalStateException("Attempting to recover a dangling head for " + vertex + " but it has in-degree > 0");

        // generate the CIGAR string from Smith-Waterman between the dangling tail and reference paths
        final DanglingChainMergeHelper danglingHeadMergeResult = generateCigarAgainstUpwardsReferencePath(vertex, pruneFactor, minDanglingBranchLength);

        // if the CIGAR is too complex (or couldn't be computed) then we do not allow the merge into the reference path
        if ( danglingHeadMergeResult == null || ! cigarIsOkayToMerge(danglingHeadMergeResult.cigar, true, false) )
            return 0;

        // merge
        return mergeDanglingHead(danglingHeadMergeResult);
    }

    /**
     * Determine whether the provided cigar is okay to merge into the reference path
     *
     * @param cigar    the cigar to analyze
     * @param requireFirstElementM if true, require that the first cigar element be an M operator in order for it to be okay
     * @param requireLastElementM  if true, require that the last cigar element be an M operator in order for it to be okay
     * @return true if it's okay to merge, false otherwise
     */
    protected boolean cigarIsOkayToMerge(final Cigar cigar, final boolean requireFirstElementM, final boolean requireLastElementM) {

        final List<CigarElement> elements = cigar.getCigarElements();
        final int numElements = elements.size();

        // don't allow more than a couple of different ops
        if ( numElements == 0 || numElements > MAX_CIGAR_COMPLEXITY )
            return false;

        // the first element must be an M
        if ( requireFirstElementM && elements.get(0).getOperator() != CigarOperator.M )
            return false;

        // the last element must be an M
        if ( requireLastElementM && elements.get(numElements - 1).getOperator() != CigarOperator.M )
            return false;

        // note that there are checks for too many mismatches in the dangling branch later in the process

        return true;
    }

    /**
     * Actually merge the dangling tail if possible
     *
     * @param danglingTailMergeResult   the result from generating a Cigar for the dangling tail against the reference
     * @return 1 if merge was successful, 0 otherwise
     */
    protected int mergeDanglingTail(final DanglingChainMergeHelper danglingTailMergeResult) {

        final List<CigarElement> elements = danglingTailMergeResult.cigar.getCigarElements();
        final CigarElement lastElement = elements.get(elements.size() - 1);
        if ( lastElement.getOperator() != CigarOperator.M )
            throw new IllegalArgumentException("The last Cigar element must be an M");

        final int lastRefIndex = danglingTailMergeResult.cigar.getReferenceLength() - 1;
        final int matchingSuffix = Math.min(GraphUtils.longestSuffixMatch(danglingTailMergeResult.referencePathString, danglingTailMergeResult.danglingPathString, lastRefIndex), lastElement.getLength());
        if ( matchingSuffix == 0 )
            return 0;

        final int altIndexToMerge = Math.max(danglingTailMergeResult.cigar.getReadLength() - matchingSuffix - 1, 0);

        // there is an important edge condition that we need to handle here: Smith-Waterman correctly calculates that there is a
        // deletion, that deletion is left-aligned such that the LCA node is part of that deletion, and the rest of the dangling
        // tail is a perfect match to the suffix of the reference path.  In this case we need to push the reference index to merge
        // down one position so that we don't incorrectly cut a base off of the deletion.
        final boolean firstElementIsDeletion = elements.get(0).getOperator() == CigarOperator.D;
        final boolean mustHandleLeadingDeletionCase =  firstElementIsDeletion && (elements.get(0).getLength() + matchingSuffix == lastRefIndex + 1);
        final int refIndexToMerge = lastRefIndex - matchingSuffix + 1 + (mustHandleLeadingDeletionCase ? 1 : 0);

        // another edge condition occurs here: if Smith-Waterman places the whole tail into an insertion then it will try to
        // merge back to the LCA, which results in a cycle in the graph.  So we do not want to merge in such a case.
        if ( refIndexToMerge == 0 )
            return 0;

        // it's safe to merge now
        addEdge(danglingTailMergeResult.danglingPath.get(altIndexToMerge), danglingTailMergeResult.referencePath.get(refIndexToMerge), ((MyEdgeFactory)getEdgeFactory()).createEdge(false, 1));

        return 1;
    }

    /**
     * Actually merge the dangling head if possible
     *
     * @param danglingHeadMergeResult   the result from generating a Cigar for the dangling head against the reference
     * @return 1 if merge was successful, 0 otherwise
     */
    protected int mergeDanglingHead(final DanglingChainMergeHelper danglingHeadMergeResult) {

        final List<CigarElement> elements = danglingHeadMergeResult.cigar.getCigarElements();
        final CigarElement firstElement = elements.get(0);
        if ( firstElement.getOperator() != CigarOperator.M )
            throw new IllegalArgumentException("The first Cigar element must be an M");

        final int indexesToMerge = bestPrefixMatch(danglingHeadMergeResult.referencePathString, danglingHeadMergeResult.danglingPathString, firstElement.getLength());
        if ( indexesToMerge <= 0 )
            return 0;

        // we can't push back the reference path
        if ( indexesToMerge >= danglingHeadMergeResult.referencePath.size() - 1 )
            return 0;

        // but we can manipulate the dangling path if we need to
        if ( indexesToMerge >= danglingHeadMergeResult.danglingPath.size() &&
                ! extendDanglingPathAgainstReference(danglingHeadMergeResult, indexesToMerge - danglingHeadMergeResult.danglingPath.size() + 2) )
            return 0;

        addEdge(danglingHeadMergeResult.referencePath.get(indexesToMerge+1), danglingHeadMergeResult.danglingPath.get(indexesToMerge), ((MyEdgeFactory)getEdgeFactory()).createEdge(false, 1));

        return 1;
    }

    /**
     * Generates the CIGAR string from the Smith-Waterman alignment of the dangling path (where the
     * provided vertex is the sink) and the reference path.
     *
     * @param vertex   the sink of the dangling chain
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @return a SmithWaterman object which can be null if no proper alignment could be generated
     */
    protected DanglingChainMergeHelper generateCigarAgainstDownwardsReferencePath(final MultiDeBruijnVertex vertex, final int pruneFactor, final int minDanglingBranchLength) {
        final int minTailPathLength = Math.max(1, minDanglingBranchLength); // while heads can be 0, tails absolutely cannot

        // find the lowest common ancestor path between this vertex and the diverging master path if available
        final List<MultiDeBruijnVertex> altPath = findPathUpwardsToLowestCommonAncestor(vertex, pruneFactor);
        if ( altPath == null || isRefSource(altPath.get(0)) || altPath.size() < minTailPathLength + 1 ) // add 1 to include the LCA
            return null;

        // now get the reference path from the LCA
        final List<MultiDeBruijnVertex> refPath = getReferencePath(altPath.get(0), TraversalDirection.downwards, Arrays.asList(incomingEdgeOf(altPath.get(1))));

        // create the Smith-Waterman strings to use
        final byte[] refBases = getBasesForPath(refPath, false);
        final byte[] altBases = getBasesForPath(altPath, false);

        // run Smith-Waterman to determine the best alignment (and remove trailing deletions since they aren't interesting)
        final SmithWaterman alignment = new SWPairwiseAlignment(refBases, altBases, SWParameterSet.STANDARD_NGS, SWPairwiseAlignment.OVERHANG_STRATEGY.LEADING_INDEL);
        return new DanglingChainMergeHelper(altPath, refPath, altBases, refBases, AlignmentUtils.removeTrailingDeletions(alignment.getCigar()));
    }

    /**
     * Generates the CIGAR string from the Smith-Waterman alignment of the dangling path (where the
     * provided vertex is the source) and the reference path.
     *
     * @param vertex   the source of the dangling head
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @return a SmithWaterman object which can be null if no proper alignment could be generated
     */
    protected DanglingChainMergeHelper generateCigarAgainstUpwardsReferencePath(final MultiDeBruijnVertex vertex, final int pruneFactor, final int minDanglingBranchLength) {

        // find the highest common descendant path between vertex and the reference source if available
        final List<MultiDeBruijnVertex> altPath = findPathDownwardsToHighestCommonDescendantOfReference(vertex, pruneFactor);
        if ( altPath == null || isRefSink(altPath.get(0)) || altPath.size() < minDanglingBranchLength + 1 ) // add 1 to include the LCA
            return null;

        // now get the reference path from the LCA
        final List<MultiDeBruijnVertex> refPath = getReferencePath(altPath.get(0), TraversalDirection.upwards, Collections.<MultiSampleEdge>emptyList());

        // create the Smith-Waterman strings to use
        final byte[] refBases = getBasesForPath(refPath, true);
        final byte[] altBases = getBasesForPath(altPath, true);

        // run Smith-Waterman to determine the best alignment (and remove trailing deletions since they aren't interesting)
        final SmithWaterman alignment = new SWPairwiseAlignment(refBases, altBases, SWParameterSet.STANDARD_NGS, SWPairwiseAlignment.OVERHANG_STRATEGY.LEADING_INDEL);
        return new DanglingChainMergeHelper(altPath, refPath, altBases, refBases, AlignmentUtils.removeTrailingDeletions(alignment.getCigar()));
    }

    /**
     * Finds the path upwards in the graph from this vertex to the first diverging node, including that (lowest common ancestor) vertex.
     * Note that nodes are excluded if their pruning weight is less than the pruning factor.
     *
     * @param vertex   the original vertex
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @return the path if it can be determined or null if this vertex either doesn't merge onto another path or
     *  has an ancestor with multiple incoming edges before hitting the reference path
     */
    protected List<MultiDeBruijnVertex> findPathUpwardsToLowestCommonAncestor(final MultiDeBruijnVertex vertex, final int pruneFactor) {
        final LinkedList<MultiDeBruijnVertex> path = new LinkedList<>();

        MultiDeBruijnVertex v = vertex;
        while ( inDegreeOf(v) == 1 && outDegreeOf(v) < 2 ) {
            final MultiSampleEdge edge = incomingEdgeOf(v);
            // if it has too low a weight, don't use it (or previous vertexes) for the path
            if ( edge.getPruningMultiplicity() < pruneFactor )
                path.clear();
            // otherwise it is safe to use
            else
                path.addFirst(v);
            v = getEdgeSource(edge);
        }
        path.addFirst(v);

        return outDegreeOf(v) > 1 ? path : null;
    }

    /**
     * Finds the path downwards in the graph from this vertex to the reference sequence, including the highest common descendant vertex.
     * However note that the path is reversed so that this vertex ends up at the end of the path.
     * Also note that nodes are excluded if their pruning weight is less than the pruning factor.
     *
     * @param vertex   the original vertex
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @return the path if it can be determined or null if this vertex either doesn't merge onto the reference path or
     *  has a descendant with multiple outgoing edges before hitting the reference path
     */
    protected List<MultiDeBruijnVertex> findPathDownwardsToHighestCommonDescendantOfReference(final MultiDeBruijnVertex vertex, final int pruneFactor) {
        final LinkedList<MultiDeBruijnVertex> path = new LinkedList<>();

        MultiDeBruijnVertex v = vertex;
        while ( ! isReferenceNode(v) && outDegreeOf(v) == 1 ) {
            final MultiSampleEdge edge = outgoingEdgeOf(v);
            // if it has too low a weight, don't use it (or previous vertexes) for the path
            if ( edge.getPruningMultiplicity() < pruneFactor )
                path.clear();
                // otherwise it is safe to use
            else
                path.addFirst(v);
            v = getEdgeTarget(edge);
        }
        path.addFirst(v);

        return isReferenceNode(v) ? path : null;
    }

    private enum TraversalDirection {
        downwards,
        upwards
    }

    /**
     * Finds the path in the graph from this vertex to the reference sink, including this vertex
     *
     * @param start   the reference vertex to start from
     * @param direction describes which direction to move in the graph (i.e. down to the reference sink or up to the source)
     * @param blacklistedEdges edges to ignore in the traversal down; useful to exclude the non-reference dangling paths
     * @return the path (non-null, non-empty)
     */
    protected List<MultiDeBruijnVertex> getReferencePath(final MultiDeBruijnVertex start,
                                                         final TraversalDirection direction,
                                                         final Collection<MultiSampleEdge> blacklistedEdges) {

        final List<MultiDeBruijnVertex> path = new ArrayList<>();

        MultiDeBruijnVertex v = start;
        while ( v != null ) {
            path.add(v);
            v = (direction == TraversalDirection.downwards ? getNextReferenceVertex(v, true, blacklistedEdges) : getPrevReferenceVertex(v));
        }

        return path;
    }

    /**
     * The base sequence for the given path.
     *
     * @param path the list of vertexes that make up the path
     * @param expandSource if true and if we encounter a source node, then expand (and reverse) the character sequence for that node
     * @return  non-null sequence of bases corresponding to the given path
     */
    public byte[] getBasesForPath(final List<MultiDeBruijnVertex> path, final boolean expandSource) {
        if ( path == null ) throw new IllegalArgumentException("Path cannot be null");

        final StringBuilder sb = new StringBuilder();
        for ( final MultiDeBruijnVertex v : path ) {
            if ( expandSource && isSource(v) ) {
                final String seq = v.getSequenceString();
                sb.append(new StringBuilder(seq).reverse().toString());
            } else {
                sb.append((char)v.getSuffix());
            }
        }

        return sb.toString().getBytes();
    }

    /**
     * Finds the index of the best extent of the prefix match between the provided paths, for dangling head merging.
     * Assumes that path1.length >= maxIndex and path2.length >= maxIndex.
     *
     * @param path1  the first path
     * @param path2  the second path
     * @param maxIndex the maximum index to traverse (not inclusive)
     * @return the index of the ideal prefix match or -1 if it cannot find one, must be less than maxIndex
     */
    protected int bestPrefixMatch(final byte[] path1, final byte[] path2, final int maxIndex) {
        final int maxMismatches = getMaxMismatches(maxIndex);
        int mismatches = 0;
        int index = 0;
        int lastGoodIndex = -1;
        while ( index < maxIndex ) {
            if ( path1[index] != path2[index] ) {
                if ( ++mismatches > maxMismatches )
                    return -1;
                lastGoodIndex = index;
            }
            index++;
        }
        // if we got here then we hit the max index
        return lastGoodIndex;
    }

    /**
     * Determine the maximum number of mismatches permitted on the branch.
     * Unless it's preset (e.g. by unit tests) it should be the length of the branch divided by the kmer size.
     *
     * @param lengthOfDanglingBranch  the length of the branch itself
     * @return positive integer
     */
    private int getMaxMismatches(final int lengthOfDanglingBranch) {
        return maxMismatchesInDanglingHead > 0 ? maxMismatchesInDanglingHead : Math.max(1, (lengthOfDanglingBranch / kmerSize));
    }

    protected boolean extendDanglingPathAgainstReference(final DanglingChainMergeHelper danglingHeadMergeResult, final int numNodesToExtend) {

        final int indexOfLastDanglingNode = danglingHeadMergeResult.danglingPath.size() - 1;
        final int indexOfRefNodeToUse = indexOfLastDanglingNode + numNodesToExtend;
        if ( indexOfRefNodeToUse >= danglingHeadMergeResult.referencePath.size() )
            return false;

        final MultiDeBruijnVertex danglingSource = danglingHeadMergeResult.danglingPath.remove(indexOfLastDanglingNode);
        final StringBuilder sb = new StringBuilder();
        final byte[] refSourceSequence = danglingHeadMergeResult.referencePath.get(indexOfRefNodeToUse).getSequence();
        for ( int i = 0; i < numNodesToExtend; i++ )
            sb.append((char)refSourceSequence[i]);
        sb.append(danglingSource.getSequenceString());
        final byte[] sequenceToExtend = sb.toString().getBytes();

        // clean up the source and edge
        final MultiSampleEdge sourceEdge = outgoingEdgeOf(danglingSource);
        MultiDeBruijnVertex prevV = getEdgeTarget(sourceEdge);
        removeEdge(danglingSource, prevV);

        // extend the path
        for ( int i = numNodesToExtend; i > 0; i-- ) {
            final MultiDeBruijnVertex newV = new MultiDeBruijnVertex(Arrays.copyOfRange(sequenceToExtend, i, i+kmerSize));
            addVertex(newV);
            final MultiSampleEdge newE = addEdge(newV, prevV);
            newE.setMultiplicity(sourceEdge.getMultiplicity());
            danglingHeadMergeResult.danglingPath.add(newV);
            prevV = newV;
        }

        return true;
    }
}