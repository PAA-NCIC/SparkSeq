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
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.log4j.Logger;
import org.ncic.bioinfo.sparkseq.algorithms.utils.CigarUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLoc;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenotypingGivenAllelesUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.haplotype.Haplotype;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.GATKSAMRecord;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.ActiveRegion;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.AssemblyResult;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.AssemblyResultSet;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.ReadErrorCorrector;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.graphs.BaseEdge;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.graphs.BaseGraph;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.graphs.BaseVertex;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.graphs.KBestHaplotype;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.graphs.KBestHaplotypeFinder;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.graphs.SeqGraph;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.graphs.SeqVertex;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Author: wbc
 */
public abstract class LocalAssemblyEngine {
    private final static Logger logger = Logger.getLogger(LocalAssemblyEngine.class);

    /**
     * If false, we will only write out a region around the reference source
     */
    private final static boolean PRINT_FULL_GRAPH_FOR_DEBUGGING = true;
    public static final byte DEFAULT_MIN_BASE_QUALITY_TO_USE = (byte) 10;
    private static final int MIN_HAPLOTYPE_REFERENCE_LENGTH = 30;

    protected final int numBestHaplotypesPerGraph;

    protected boolean debug = false;
    protected boolean allowCyclesInKmerGraphToGeneratePaths = false;
    protected boolean debugGraphTransformations = false;
    protected boolean recoverDanglingBranches = true;
    protected int minDanglingBranchLength = 0;

    protected byte minBaseQualityToUseInAssembly = DEFAULT_MIN_BASE_QUALITY_TO_USE;
    protected int pruneFactor = 2;
    protected boolean errorCorrectKmers = false;

    private PrintStream graphWriter = null;

    /**
     * Create a new LocalAssemblyEngine with all default parameters, ready for use
     *
     * @param numBestHaplotypesPerGraph the number of haplotypes to generate for each assembled graph
     */
    protected LocalAssemblyEngine(final int numBestHaplotypesPerGraph) {
        if (numBestHaplotypesPerGraph < 1)
            throw new IllegalArgumentException("numBestHaplotypesPerGraph should be >= 1 but got " + numBestHaplotypesPerGraph);
        this.numBestHaplotypesPerGraph = numBestHaplotypesPerGraph;
    }

    /**
     * Main subclass function: given reads and a reference haplotype give us graphs to use for constructing
     * non-reference haplotypes.
     *
     * @param reads        the reads we're going to assemble
     * @param refHaplotype the reference haplotype
     * @return a non-null list of reads
     */
    protected abstract List<AssemblyResult> assemble(List<GATKSAMRecord> reads, Haplotype refHaplotype, List<Haplotype> givenHaplotypes);

    /**
     * Main entry point into the assembly engine. Build a set of deBruijn graphs out of the provided reference sequence and list of reads
     *
     * @param activeRegion             ActiveRegion object holding the reads which are to be used during assembly
     * @param refHaplotype             reference haplotype object
     * @param fullReferenceWithPadding byte array holding the reference sequence with padding
     * @param refLoc                   GenomeLoc object corresponding to the reference sequence with padding
     * @param givenAlleles             the alleles to inject into the haplotypes during GGA mode
     * @param readErrorCorrector       a ReadErrorCorrector object, if read are to be corrected before assembly. Can be null if no error corrector is to be used.
     * @return the resulting assembly-result-set
     */
    public AssemblyResultSet runLocalAssembly(final ActiveRegion activeRegion,
                                              final Haplotype refHaplotype,
                                              final byte[] fullReferenceWithPadding,
                                              final GenomeLoc refLoc,
                                              final List<VariantContext> givenAlleles,
                                              final ReadErrorCorrector readErrorCorrector) {
        if (activeRegion == null) {
            throw new IllegalArgumentException("Assembly engine cannot be used with a null ActiveRegion.");
        }
        if (activeRegion.getExtendedLoc() == null) {
            throw new IllegalArgumentException("Active region must have an extended location.");
        }
        if (refHaplotype == null) {
            throw new IllegalArgumentException("Reference haplotype cannot be null.");
        }
        if (fullReferenceWithPadding.length != refLoc.size()) {
            throw new IllegalArgumentException("Reference bases and reference loc must be the same size.");
        }
        if (pruneFactor < 0) {
            throw new IllegalArgumentException("Pruning factor cannot be negative");
        }

        // create the list of artificial haplotypes that should be added to the graph for GGA mode
        final List<Haplotype> givenHaplotypes = GenotypingGivenAllelesUtils.composeGivenHaplotypes(refHaplotype, givenAlleles, activeRegion.getExtendedLoc());

        // error-correct reads before clipping low-quality tails: some low quality bases might be good and we want to recover them
        final List<GATKSAMRecord> correctedReads;
        if (readErrorCorrector != null) {
            // now correct all reads in active region after filtering/downsampling
            // Note that original reads in active region are NOT modified by default, since they will be used later for GL computation,
            // and we only want the read-error corrected reads for graph building.
            readErrorCorrector.addReadsToKmers(activeRegion.getReads());
            correctedReads = new ArrayList<>(readErrorCorrector.correctReads(activeRegion.getReads()));
        } else {
            correctedReads = activeRegion.getReads();
        }

        final List<SeqGraph> nonRefGraphs = new LinkedList<>();
        final AssemblyResultSet resultSet = new AssemblyResultSet();
        resultSet.setRegionForGenotyping(activeRegion);
        resultSet.setFullReferenceWithPadding(fullReferenceWithPadding);
        resultSet.setPaddedReferenceLoc(refLoc);
        final GenomeLoc activeRegionExtendedLocation = activeRegion.getExtendedLoc();
        refHaplotype.setGenomeLocation(activeRegionExtendedLocation);
        resultSet.add(refHaplotype);
        final Map<SeqGraph, AssemblyResult> assemblyResultByGraph = new HashMap<>();
        // create the graphs by calling our subclass assemble method
        for (final AssemblyResult result : assemble(correctedReads, refHaplotype, givenHaplotypes)) {
            if (result.getStatus() == AssemblyResult.Status.ASSEMBLED_SOME_VARIATION) {
                // do some QC on the graph
                sanityCheckGraph(result.getGraph(), refHaplotype);
                // add it to graphs with meaningful non-reference features
                assemblyResultByGraph.put(result.getGraph(), result);
                nonRefGraphs.add(result.getGraph());
            }

        }

        findBestPaths(nonRefGraphs, refHaplotype, refLoc, activeRegionExtendedLocation, assemblyResultByGraph, resultSet);

        // print the graphs if the appropriate debug option has been turned on
        if (graphWriter != null) {
            printGraphs(nonRefGraphs);
        }

        return resultSet;
    }

    protected List<Haplotype> findBestPaths(final List<SeqGraph> graphs, final Haplotype refHaplotype, final GenomeLoc refLoc, final GenomeLoc activeRegionWindow,
                                            final Map<SeqGraph, AssemblyResult> assemblyResultByGraph, final AssemblyResultSet assemblyResultSet) {
        // add the reference haplotype separately from all the others to ensure that it is present in the list of haplotypes
        final Set<Haplotype> returnHaplotypes = new LinkedHashSet<>();

        final int activeRegionStart = refHaplotype.getAlignmentStartHapwrtRef();
        final ArrayList<KBestHaplotypeFinder> finders = new ArrayList<>(graphs.size());
        int failedCigars = 0;

        for (final SeqGraph graph : graphs) {
            final SeqVertex source = graph.getReferenceSourceVertex();
            final SeqVertex sink = graph.getReferenceSinkVertex();
            if (source == null || sink == null)
                throw new IllegalArgumentException("Both source and sink cannot be null but got " + source + " and sink " + sink + " for graph " + graph);
            final KBestHaplotypeFinder haplotypeFinder = new KBestHaplotypeFinder(graph, source, sink);
            finders.add(haplotypeFinder);
            final Iterator<KBestHaplotype> bestHaplotypes = haplotypeFinder.iterator(numBestHaplotypesPerGraph);

            while (bestHaplotypes.hasNext()) {
                final KBestHaplotype kBestHaplotype = bestHaplotypes.next();
                final Haplotype h = kBestHaplotype.haplotype();
                if (!returnHaplotypes.contains(h)) {
                    final Cigar cigar = CigarUtils.calculateCigar(refHaplotype.getBases(), h.getBases());

                    if (cigar == null) {
                        failedCigars++; // couldn't produce a meaningful alignment of haplotype to reference, fail quietly
                        continue;
                    } else if (cigar.isEmpty()) {
                        throw new IllegalStateException("Smith-Waterman alignment failure. Cigar = " + cigar + " with reference length " + cigar.getReferenceLength() +
                                " but expecting reference length of " + refHaplotype.getCigar().getReferenceLength());
                    } else if (pathIsTooDivergentFromReference(cigar) || cigar.getReferenceLength() < MIN_HAPLOTYPE_REFERENCE_LENGTH) {
                        // N cigar elements means that a bubble was too divergent from the reference so skip over this path
                        continue;
                    } else if (cigar.getReferenceLength() != refHaplotype.getCigar().getReferenceLength()) { // SW failure
                        throw new IllegalStateException("Smith-Waterman alignment failure. Cigar = " + cigar + " with reference length "
                                + cigar.getReferenceLength() + " but expecting reference length of " + refHaplotype.getCigar().getReferenceLength()
                                + " ref = " + refHaplotype + " path " + new String(h.getBases()));
                    }

                    h.setCigar(cigar);
                    h.setAlignmentStartHapwrtRef(activeRegionStart);
                    h.setGenomeLocation(activeRegionWindow);
                    returnHaplotypes.add(h);
                    assemblyResultSet.add(h, assemblyResultByGraph.get(graph));

                    if (debug)
                        logger.info("Adding haplotype " + h.getCigar() + " from graph with kmer " + graph.getKmerSize());
                }
            }
        }

        // Make sure that the ref haplotype is amongst the return haplotypes and calculate its score as
        // the first returned by any finder.
        if (!returnHaplotypes.contains(refHaplotype)) {
            double refScore = Double.NaN;
            for (final KBestHaplotypeFinder finder : finders) {
                final double candidate = finder.score(refHaplotype);
                if (Double.isNaN(candidate)) continue;
                refScore = candidate;
                break;
            }
            refHaplotype.setScore(refScore);
            returnHaplotypes.add(refHaplotype);
        }

        if (failedCigars != 0)
            logger.debug(String.format("failed to align some haplotypes (%d) back to the reference (loc=%s); these will be ignored.", failedCigars, refLoc.toString()));

        if (debug) {
            if (returnHaplotypes.size() > 1) {
                logger.info("Found " + returnHaplotypes.size() + " candidate haplotypes of " + returnHaplotypes.size() + " possible combinations to evaluate every read against.");
            } else {
                logger.info("Found only the reference haplotype in the assembly graph.");
            }
            for (final Haplotype h : returnHaplotypes) {
                logger.info(h.toString());
                logger.info("> Cigar = " + h.getCigar() + " : " + h.getCigar().getReferenceLength() + " score " + h.getScore() + " ref " + h.isReference());
            }
        }

        return new ArrayList<>(returnHaplotypes);

    }

    /**
     * We use CigarOperator.N as the signal that an incomplete or too divergent bubble was found during bubble traversal
     *
     * @param c the cigar to test
     * @return true if we should skip over this path
     */
    private boolean pathIsTooDivergentFromReference(final Cigar c) {
        for (final CigarElement ce : c.getCigarElements()) {
            if (ce.getOperator().equals(CigarOperator.N)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Print graph to file if debugGraphTransformations is enabled
     *
     * @param graph the graph to print
     * @param file  the destination file
     */
    protected void printDebugGraphTransform(final BaseGraph graph, final File file) {
        if (debugGraphTransformations) {
            if (PRINT_FULL_GRAPH_FOR_DEBUGGING)
                graph.printGraph(file, pruneFactor);
            else
                graph.subsetToRefSource().printGraph(file, pruneFactor);
        }
    }

    protected AssemblyResult cleanupSeqGraph(final SeqGraph seqGraph) {
        printDebugGraphTransform(seqGraph, new File("sequenceGraph.1.dot"));

        // the very first thing we need to do is zip up the graph, or pruneGraph will be too aggressive
        seqGraph.zipLinearChains();
        printDebugGraphTransform(seqGraph, new File("sequenceGraph.2.zipped.dot"));

        // now go through and prune the graph, removing vertices no longer connected to the reference chain
        seqGraph.removeSingletonOrphanVertices();
        seqGraph.removeVerticesNotConnectedToRefRegardlessOfEdgeDirection();

        printDebugGraphTransform(seqGraph, new File("sequenceGraph.3.pruned.dot"));
        seqGraph.simplifyGraph();
        printDebugGraphTransform(seqGraph, new File("sequenceGraph.4.merged.dot"));

        // The graph has degenerated in some way, so the reference source and/or sink cannot be id'd.  Can
        // happen in cases where for example the reference somehow manages to acquire a cycle, or
        // where the entire assembly collapses back into the reference sequence.
        if (seqGraph.getReferenceSourceVertex() == null || seqGraph.getReferenceSinkVertex() == null)
            return new AssemblyResult(AssemblyResult.Status.JUST_ASSEMBLED_REFERENCE, seqGraph);

        seqGraph.removePathsNotConnectedToRef();
        seqGraph.simplifyGraph();
        if (seqGraph.vertexSet().size() == 1) {
            // we've perfectly assembled into a single reference haplotype, add a empty seq vertex to stop
            // the code from blowing up.
            // TODO -- ref properties should really be on the vertices, not the graph itself
            final SeqVertex complete = seqGraph.vertexSet().iterator().next();
            final SeqVertex dummy = new SeqVertex("");
            seqGraph.addVertex(dummy);
            seqGraph.addEdge(complete, dummy, new BaseEdge(true, 0));
        }
        printDebugGraphTransform(seqGraph, new File("sequenceGraph.5.final.dot"));
        return new AssemblyResult(AssemblyResult.Status.ASSEMBLED_SOME_VARIATION, seqGraph);
    }

    /**
     * Perform general QC on the graph to make sure something hasn't gone wrong during assembly
     *
     * @param graph        the graph to check
     * @param refHaplotype the reference haplotype
     */
    private <T extends BaseVertex, E extends BaseEdge> void sanityCheckGraph(final BaseGraph<T, E> graph, final Haplotype refHaplotype) {
        sanityCheckReferenceGraph(graph, refHaplotype);
    }

    /**
     * Make sure the reference sequence is properly represented in the provided graph
     *
     * @param graph        the graph to check
     * @param refHaplotype the reference haplotype
     */
    private <T extends BaseVertex, E extends BaseEdge> void sanityCheckReferenceGraph(final BaseGraph<T, E> graph, final Haplotype refHaplotype) {
        if (graph.getReferenceSourceVertex() == null) {
            throw new IllegalStateException("All reference graphs must have a reference source vertex.");
        }
        if (graph.getReferenceSinkVertex() == null) {
            throw new IllegalStateException("All reference graphs must have a reference sink vertex.");
        }
        if (!Arrays.equals(graph.getReferenceBytes(graph.getReferenceSourceVertex(), graph.getReferenceSinkVertex(), true, true), refHaplotype.getBases())) {
            throw new IllegalStateException("Mismatch between the reference haplotype and the reference assembly graph path. for graph " + graph +
                    " graph = " + new String(graph.getReferenceBytes(graph.getReferenceSourceVertex(), graph.getReferenceSinkVertex(), true, true)) +
                    " haplotype = " + new String(refHaplotype.getBases())
            );
        }
    }

    /**
     * Print the generated graphs to the graphWriter
     *
     * @param graphs a non-null list of graphs to print out
     */
    private void printGraphs(final List<SeqGraph> graphs) {
        final int writeFirstGraphWithSizeSmallerThan = 50;

        graphWriter.println("digraph assemblyGraphs {");
        for (final SeqGraph graph : graphs) {
            if (debugGraphTransformations && graph.getKmerSize() >= writeFirstGraphWithSizeSmallerThan) {
                logger.info("Skipping writing of graph with kmersize " + graph.getKmerSize());
                continue;
            }

            graph.printGraph(graphWriter, false, pruneFactor);

            if (debugGraphTransformations)
                break;
        }

        graphWriter.println("}");
    }

    // -----------------------------------------------------------------------------------------------
    //
    // getter / setter routines for generic assembler properties
    //
    // -----------------------------------------------------------------------------------------------

    public int getPruneFactor() {
        return pruneFactor;
    }

    public void setPruneFactor(int pruneFactor) {
        this.pruneFactor = pruneFactor;
    }

    public boolean shouldErrorCorrectKmers() {
        return errorCorrectKmers;
    }

    public void setErrorCorrectKmers(boolean errorCorrectKmers) {
        this.errorCorrectKmers = errorCorrectKmers;
    }

    public void setGraphWriter(PrintStream graphWriter) {
        this.graphWriter = graphWriter;
    }

    public byte getMinBaseQualityToUseInAssembly() {
        return minBaseQualityToUseInAssembly;
    }

    public void setMinBaseQualityToUseInAssembly(byte minBaseQualityToUseInAssembly) {
        this.minBaseQualityToUseInAssembly = minBaseQualityToUseInAssembly;
    }

    public boolean isDebug() {
        return debug;
    }

    public void setDebug(boolean debug) {
        this.debug = debug;
    }

    public boolean isAllowCyclesInKmerGraphToGeneratePaths() {
        return allowCyclesInKmerGraphToGeneratePaths;
    }

    public void setAllowCyclesInKmerGraphToGeneratePaths(boolean allowCyclesInKmerGraphToGeneratePaths) {
        this.allowCyclesInKmerGraphToGeneratePaths = allowCyclesInKmerGraphToGeneratePaths;
    }

    public boolean isDebugGraphTransformations() {
        return debugGraphTransformations;
    }

    public void setDebugGraphTransformations(boolean debugGraphTransformations) {
        this.debugGraphTransformations = debugGraphTransformations;
    }

    public boolean isRecoverDanglingBranches() {
        return recoverDanglingBranches;
    }

    public void setRecoverDanglingBranches(final boolean recoverDanglingBranches) {
        this.recoverDanglingBranches = recoverDanglingBranches;
    }

    public void setMinDanglingBranchLength(final int minDanglingBranchLength) {
        this.minDanglingBranchLength = minDanglingBranchLength;
    }
}
