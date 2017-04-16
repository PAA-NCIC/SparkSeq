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

import org.ncic.bioinfo.sparkseq.algorithms.utils.Utils;
import org.ncic.bioinfo.sparkseq.algorithms.data.basic.Pair;

import java.util.*;

/**
* General recursive sub-haplotype finder.
* <p>
*   Provides the k-best sub-haplotypes from a vertex provided map between outgoing edges and its target finders
* </p>
* <p>
*  This is done efficiently by keeping an priority-queue on best sub-haplotype solutions and pulling them on demand
*  as needed.
* </p>
* <p>
*  Solutions are cached for repeated retrieval so that we save compute at vertices that share sub-haplotypes
*     (share descendant vertices). This aspect is controlled by {@link KBestSubHaplotypeFinder} that instantiate
*     a unique {@link KBestSubHaplotypeFinder} for each vertex in the graph that belongs to a valid path
*     between the source and sink node.
* </p>
*
* @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
*/
class RecursiveSubHaplotypeFinder extends AggregatedSubHaplotypeFinder<RecursiveSubHaplotypeFinder.EdgeSubHaplotypeFinder> {


    private final SeqVertex vertex;
    private final boolean isReference;

    /**
     * Creates a recursive sub-haplotype finder give the target graph, first vertex and all possible outgoing edges
     *  with the corresponding sub-sub-haplotype finders.
     *
     * <p>For efficiency shake, it will not verify the content of {@code children} map; i.e. that indeed all keys
     * are outgoing edges from {@code vertex} on {@code graph} and that the value sub-haplotype resolver have as
     * the first vertex the adjacent vertex through that key edge.</p>
     *
     * @param vertex first vertex for all sub-haplotype solutions provided by this finder
     * @param children map from outgoing edge to the corresponding sub-sub-haplotype finder.
     */
    public RecursiveSubHaplotypeFinder(final SeqGraph graph, final SeqVertex vertex,
                                       final Map<BaseEdge, KBestSubHaplotypeFinder> children) {
        super(createChildFinderCollection(vertex, children));
        this.vertex = vertex;
        this.isReference = graph.isReferenceNode(vertex);
    }

    /**
     * Wraps the descendant vertices finders in order to take advantage of the {@link AggregatedSubHaplotypeFinder}
     * common code.
     * <p>
     * Automatically calibrates the edge score (cost) so that it takes into account the total across all edges.
     * </p>
     *
     * @param vertex the parent vertex.
     * @param finders the child vertices indexed by the connecting edge.
     * @return never {@code null} but potentially an empty collection if there is child returning some sub-haplotype
     *    solution.
     */
    private static Collection<EdgeSubHaplotypeFinder> createChildFinderCollection(final SeqVertex vertex,
                                                             final Map<BaseEdge,KBestSubHaplotypeFinder> finders) {
        if (finders == null) throw new IllegalArgumentException("the edge to child map cannot be null");
        final ArrayList<EdgeSubHaplotypeFinder> result = new ArrayList<>(finders.size());
        for (final Map.Entry<BaseEdge,KBestSubHaplotypeFinder> e : finders.entrySet()) {
            final EdgeSubHaplotypeFinder subFinder = new EdgeSubHaplotypeFinder(vertex,e.getKey(), e.getValue());
            if (subFinder.getCount() == 0) continue;
            result.add(subFinder);
        }
        if (result.size() == 0)
            return Collections.emptySet();
        else if (result.size() == 1) // no calibration needed, by default edgeScore is 0.
            return Collections.singleton(result.get(0));
        else {
            double totalEdgeMultiplicityAcrossEdges = 0;
            for (final EdgeSubHaplotypeFinder finder : result)
                totalEdgeMultiplicityAcrossEdges += Math.max(0.5,finder.edge.getMultiplicity());
            final double log10TotalEdgeMultiplicityAcrossEdges = Math.log10(totalEdgeMultiplicityAcrossEdges);
            for (final EdgeSubHaplotypeFinder finder : result)
                finder.calibrateEdgeScore(log10TotalEdgeMultiplicityAcrossEdges);
            return result;
        }
    }

    @Override
    public boolean isReference() {
        return isReference;
    }

    @Override
    public String label() {
        return vertex.getSequenceString();
    }

    @Override
    public Set<Pair<? extends KBestSubHaplotypeFinder, String>> subFinderLabels() {
        final Set<Pair<? extends KBestSubHaplotypeFinder,String>> result = new LinkedHashSet<>(subFinders.size());
        for (final EdgeSubHaplotypeFinder subFinder : subFinders)
            result.add(new Pair<>(subFinder,simplifyZeros(String.format("%.4f", subFinder.edgeScore))));
        return result;
    }

    /**
     * Removes zeros decimal positions from edge-labels.
     *
     * @param edgeLabel the original label to reformat.
     * @return never {@code null}, the reformatted label.
     */
    private String simplifyZeros(final String edgeLabel) {
        if (edgeLabel.equals("0.000") || edgeLabel.equals("-0.000") )
            return "0.";
        int i = edgeLabel.length() - 1;
        while (edgeLabel.charAt(i) == '0')
            i--;
        return (i == edgeLabel.length() - 1) ? edgeLabel : edgeLabel.substring(0,i);
    }

    protected static class EdgeSubHaplotypeFinder implements KBestSubHaplotypeFinder {

        private final KBestSubHaplotypeFinder childFinder;

        private final SeqVertex vertex;

        private final BaseEdge edge;

        private double edgeScore = 0;

        private EdgeSubHaplotypeFinder(final SeqVertex vertex, final BaseEdge edge, final KBestSubHaplotypeFinder childFinder) {
            this.childFinder = childFinder;
            this.edge = edge;
            this.vertex = vertex;
            this.edgeScore = 0;
        }

        private void calibrateEdgeScore(final double log10TotalMultiplicityAcrossOutgoingEdges) {
            edgeScore = Math.log10(Math.max(edge.getMultiplicity(),0.5)) - log10TotalMultiplicityAcrossOutgoingEdges;
        }

        @Override
        public String id() {
            return childFinder.id();
        }

        @Override
        public String label() {
            return childFinder.label();
        }

        @Override
        public Set<Pair<? extends KBestSubHaplotypeFinder, String>> subFinderLabels() {
            return childFinder.subFinderLabels();
        }

        @Override
        public int getCount() {
            return childFinder.getCount();
        }

        @Override
        public KBestHaplotype getKBest(int k) {
            return new ChildKBestSubHaplotype(vertex,edge,childFinder.getKBest(k),edgeScore);
        }

        @Override
        public boolean isReference() {
            return childFinder.isReference();
        }

        @Override
        public double score(final byte[] bases, final int offset, final int length) {
            if (length == 0)
                return 0;
            final byte[] vertexSequence = vertex.getSequence();
            if (length < vertexSequence.length) // query is not long enough to have any score.
                return Double.NaN;
            else if (!Utils.equalRange(vertexSequence,0,bases,offset,vertexSequence.length))
                return Double.NaN;
            else
                return edgeScore + childFinder.score(bases,offset + vertexSequence.length,length - vertexSequence.length);
        }
    }

    @Override
    public String id() {
        return "v" + vertex.getId();
    }

    /**
     * Custom extension of the {@link KBestHaplotype} used for solutions generated by this class.
     *
     * <p>
     *     These by delegating on the encapsulated solution from outgoing edge's finder by adding
     *     the edge score and prefixing this outer finder
     *     source vertex.
     * </p>
     */
    private static class ChildKBestSubHaplotype extends KBestHaplotype {

        private final double score;
        private final KBestHaplotype child;
        private final SeqVertex vertex;
        private final boolean isReference;


        public ChildKBestSubHaplotype(final SeqVertex vertex, final BaseEdge edge, final KBestHaplotype child, final double edgeScore) {
            this.score = edgeScore + child.score();
            this.vertex = vertex;
            this.child = child;
            this.isReference = edge.isRef() && child.isReference();
        }

        @Override
        public SeqGraph graph() {
            return child.graph();
        }

        @Override
        public double score() {
            return score;
        }

        @Override
        public int rank() {
            return child.rank();
        }

        @Override
        protected SeqVertex head() {
            return vertex;
        }

        @Override
        protected KBestHaplotype tail() {
            return child;
        }

        @Override
        public boolean isReference() {
            return isReference;
        }
    }
}
