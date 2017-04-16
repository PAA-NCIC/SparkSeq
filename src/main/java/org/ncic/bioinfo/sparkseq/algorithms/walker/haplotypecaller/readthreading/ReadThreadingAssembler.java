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
import org.ncic.bioinfo.sparkseq.algorithms.utils.MathUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.haplotype.Haplotype;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.GATKSAMRecord;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.AssemblyResult;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.graphs.SeqGraph;

import java.io.File;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

/**
 * Author: wbc
 */
public class ReadThreadingAssembler extends LocalAssemblyEngine {
    private final static Logger logger = Logger.getLogger(ReadThreadingAssembler.class);

    private final static int DEFAULT_NUM_PATHS_PER_GRAPH = 128;
    private final static int GGA_MODE_ARTIFICIAL_COUNTS = 1000;
    private final static int KMER_SIZE_ITERATION_INCREASE = 10;
    private final static int MAX_KMER_ITERATIONS_TO_ATTEMPT = 6;

    /**
     * The min and max kmer sizes to try when building the graph.
     */
    private final List<Integer> kmerSizes;

    private final boolean dontIncreaseKmerSizesForCycles;
    private final boolean allowNonUniqueKmersInRef;
    private final int numPruningSamples;
    protected boolean removePathsNotConnectedToRef = true;
    private boolean justReturnRawGraph = false;

    /**
     * for testing only
     */
    public ReadThreadingAssembler() {
        this(DEFAULT_NUM_PATHS_PER_GRAPH, Arrays.asList(25));
    }

    public ReadThreadingAssembler(final int maxAllowedPathsForReadThreadingAssembler, final List<Integer> kmerSizes, final boolean dontIncreaseKmerSizesForCycles, final boolean allowNonUniqueKmersInRef, final int numPruningSamples) {
        super(maxAllowedPathsForReadThreadingAssembler);
        this.kmerSizes = kmerSizes;
        this.dontIncreaseKmerSizesForCycles = dontIncreaseKmerSizesForCycles;
        this.allowNonUniqueKmersInRef = allowNonUniqueKmersInRef;
        this.numPruningSamples = numPruningSamples;
    }

    protected ReadThreadingAssembler(final int maxAllowedPathsForReadThreadingAssembler, final List<Integer> kmerSizes) {
        this(maxAllowedPathsForReadThreadingAssembler, kmerSizes, true, true, 1);
    }

    /**
     * for testing purposes
     */
    protected void setJustReturnRawGraph(boolean justReturnRawGraph) {
        this.justReturnRawGraph = justReturnRawGraph;
    }

    private void addResult(final List<AssemblyResult> results, final AssemblyResult maybeNullResult) {
        if (maybeNullResult != null)
            results.add(maybeNullResult);
    }

    @Override
    public List<AssemblyResult> assemble(final List<GATKSAMRecord> reads, final Haplotype refHaplotype, final List<Haplotype> givenHaplotypes) {
        final List<AssemblyResult> results = new LinkedList<>();

        // first, try using the requested kmer sizes
        for (final int kmerSize : kmerSizes) {
            addResult(results, createGraph(reads, refHaplotype, kmerSize, givenHaplotypes, dontIncreaseKmerSizesForCycles, allowNonUniqueKmersInRef));
        }

        // if none of those worked, iterate over larger sizes if allowed to do so
        if (results.isEmpty() && !dontIncreaseKmerSizesForCycles) {
            int kmerSize = MathUtils.arrayMaxInt(kmerSizes) + KMER_SIZE_ITERATION_INCREASE;
            int numIterations = 1;
            while (results.isEmpty() && numIterations <= MAX_KMER_ITERATIONS_TO_ATTEMPT) {
                // on the last attempt we will allow low complexity graphs
                final boolean lastAttempt = numIterations == MAX_KMER_ITERATIONS_TO_ATTEMPT;
                addResult(results, createGraph(reads, refHaplotype, kmerSize, givenHaplotypes, lastAttempt, lastAttempt));
                kmerSize += KMER_SIZE_ITERATION_INCREASE;
                numIterations++;
            }
        }

        return results;
    }

    /**
     * Creates the sequence graph for the given kmerSize
     *
     * @param reads                    reads to use
     * @param refHaplotype             reference haplotype
     * @param kmerSize                 kmer size
     * @param activeAlleleHaplotypes   the GGA haplotypes to inject into the graph
     * @param allowLowComplexityGraphs if true, do not check for low-complexity graphs
     * @param allowNonUniqueKmersInRef if true, do not fail if the reference has non-unique kmers
     * @return sequence graph or null if one could not be created (e.g. because it contains cycles or too many paths or is low complexity)
     */
    protected AssemblyResult createGraph(final List<GATKSAMRecord> reads,
                                         final Haplotype refHaplotype,
                                         final int kmerSize,
                                         final List<Haplotype> activeAlleleHaplotypes,
                                         final boolean allowLowComplexityGraphs,
                                         final boolean allowNonUniqueKmersInRef) {
        if (refHaplotype.length() < kmerSize) {
            // happens in cases where the assembled region is just too small
            return new AssemblyResult(AssemblyResult.Status.FAILED, null);
        }

        if (!allowNonUniqueKmersInRef && !ReadThreadingGraph.determineNonUniqueKmers(new SequenceForKmers("ref", refHaplotype.getBases(), 0, refHaplotype.getBases().length, 1, true), kmerSize).isEmpty()) {
            if (debug)
                logger.info("Not using kmer size of " + kmerSize + " in read threading assembler because reference contains non-unique kmers");
            return null;
        }

        final ReadThreadingGraph rtgraph = new ReadThreadingGraph(kmerSize, debugGraphTransformations, minBaseQualityToUseInAssembly, numPruningSamples);

        rtgraph.setThreadingStartOnlyAtExistingVertex(!recoverDanglingBranches);

        // add the reference sequence to the graph
        rtgraph.addSequence("ref", refHaplotype.getBases(), true);

        // add the artificial GGA haplotypes to the graph
        int hapCount = 0;
        for (final Haplotype h : activeAlleleHaplotypes) {
            rtgraph.addSequence("activeAllele" + hapCount++, h.getBases(), GGA_MODE_ARTIFICIAL_COUNTS, false);
        }

        // Next pull kmers out of every read and throw them on the graph
        for (final GATKSAMRecord read : reads) {
            rtgraph.addRead(read);
        }

        // actually build the read threading graph
        rtgraph.buildGraphIfNecessary();

        // sanity check: make sure there are no cycles in the graph
        if (rtgraph.hasCycles()) {
            if (debug)
                logger.info("Not using kmer size of " + kmerSize + " in read threading assembler because it contains a cycle");
            return null;
        }

        // sanity check: make sure the graph had enough complexity with the given kmer
        if (!allowLowComplexityGraphs && rtgraph.isLowComplexity()) {
            if (debug)
                logger.info("Not using kmer size of " + kmerSize + " in read threading assembler because it does not produce a graph with enough complexity");
            return null;
        }

        printDebugGraphTransform(rtgraph, new File("" + refHaplotype.getGenomeLocation() + "-sequenceGraph." + kmerSize + ".0.0.raw_readthreading_graph.dot"));

        // go through and prune all of the chains where all edges have <= pruneFactor.  This must occur
        // before recoverDanglingTails in the graph, so that we don't spend a ton of time recovering
        // tails that we'll ultimately just trim away anyway, as the dangling tail edges have weight of 1
        rtgraph.pruneLowWeightChains(pruneFactor);

        // look at all chains in the graph that terminate in a non-ref node (dangling sources and sinks) and see if
        // we can recover them by merging some N bases from the chain back into the reference
        if (recoverDanglingBranches) {
            rtgraph.recoverDanglingTails(pruneFactor, minDanglingBranchLength);
            rtgraph.recoverDanglingHeads(pruneFactor, minDanglingBranchLength);
        }

        // remove all heading and trailing paths
        if (removePathsNotConnectedToRef) rtgraph.removePathsNotConnectedToRef();

        printDebugGraphTransform(rtgraph, new File("" + refHaplotype.getGenomeLocation() + "-sequenceGraph." + kmerSize + ".0.1.cleaned_readthreading_graph.dot"));

        final SeqGraph initialSeqGraph = rtgraph.convertToSequenceGraph();
        if (debugGraphTransformations)
            initialSeqGraph.printGraph(new File("" + refHaplotype.getGenomeLocation() + "-sequenceGraph." + kmerSize + ".0.1.initial_seqgraph.dot"), 10000);

        // if the unit tests don't want us to cleanup the graph, just return the raw sequence graph
        if (justReturnRawGraph)
            return new AssemblyResult(AssemblyResult.Status.ASSEMBLED_SOME_VARIATION, initialSeqGraph);

        if (debug)
            logger.info("Using kmer size of " + rtgraph.getKmerSize() + " in read threading assembler");
        printDebugGraphTransform(initialSeqGraph, new File("" + refHaplotype.getGenomeLocation() + "-sequenceGraph." + kmerSize + ".0.2.initial_seqgraph.dot"));
        initialSeqGraph.cleanNonRefPaths(); // TODO -- I don't this is possible by construction

        final AssemblyResult cleaned = cleanupSeqGraph(initialSeqGraph);
        final AssemblyResult.Status status = cleaned.getStatus();
        final AssemblyResult result = new AssemblyResult(status, cleaned.getGraph());
        result.setThreadingGraph(rtgraph);
        return result;
    }

    @Override
    public String toString() {
        return "ReadThreadingAssembler{" +
                "kmerSizes=" + kmerSizes +
                '}';
    }
}