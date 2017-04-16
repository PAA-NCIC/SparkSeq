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
package org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.apache.commons.lang.ArrayUtils;
import org.apache.log4j.Logger;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RODNames;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.Tags;
import org.ncic.bioinfo.sparkseq.algorithms.engine.ActiveRegionWalker;
import org.ncic.bioinfo.sparkseq.algorithms.engine.Walker;
import org.ncic.bioinfo.sparkseq.algorithms.utils.AlignmentUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.FragmentUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLoc;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLocParser;
import org.ncic.bioinfo.sparkseq.algorithms.utils.QualityUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.ReadUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.clip.ReadClipper;
import org.ncic.bioinfo.sparkseq.algorithms.utils.downsampling.DownsamplingUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.fragments.FragmentCollection;
import org.ncic.bioinfo.sparkseq.algorithms.utils.haplotype.Haplotype;
import org.ncic.bioinfo.sparkseq.algorithms.utils.pairhmm.PairHMM;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefMetaDataTracker;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.GATKSAMRecord;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.BadMateFilter;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RodBinding;
import org.ncic.bioinfo.sparkseq.algorithms.utils.vcfWriter.GVCFWriter;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.afcalculate.FixedAFCalculatorProvider;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.VariantAnnotatorEngine;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.interfaces.AnnotatorCompatible;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.argcollection.DbsnpArgumentCollection;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.argcollection.HaplotypeCallerArgumentCollection;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.genotyper.HaplotypeCallerGenotypingEngine;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.genotyper.LDMerger;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.genotyper.MergeVariantsAcrossHaplotypes;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.model.ReferenceConfidenceMode;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.model.ReferenceConfidenceModel;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.readlikelihood.PairHMMLikelihoodCalculationEngine;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.readlikelihood.ReadLikelihoodCalculationEngine;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.readlikelihood.ReadLikelihoods;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.readthreading.LocalAssemblyEngine;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.ncic.bioinfo.sparkseq.exceptions.PipelineException;
import org.ncic.bioinfo.sparkseq.exceptions.UserException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

/**
 * Author: wbc
 */
public class HaplotypeCaller implements AnnotatorCompatible {

    final protected static Logger logger = Logger.getLogger(HaplotypeCaller.class);

    // -----------------------------------------------------------------------------------------------
    // general haplotype caller arguments
    // -----------------------------------------------------------------------------------------------

    protected ReadLikelihoodCalculationEngine.Implementation likelihoodEngineImplementation = ReadLikelihoodCalculationEngine.Implementation.PairHMM;

    protected HeterogeneousKmerSizeResolution heterogeneousKmerSizeResolution = HeterogeneousKmerSizeResolution.COMBO_MIN;

    /**
     * If set, certain "early exit" optimizations in HaplotypeCaller, which aim to save compute and time by skipping
     * calculations if an ActiveRegion is determined to contain no variants, will be disabled. This is most likely to be useful if
     * you're using the -bamout argument to examine the placement of reads following reassembly and are interested in seeing the mapping of
     * reads in regions with no variations. Setting the -forceActive and -dontTrimActiveRegions flags may also be necessary.
     */
    private boolean disableOptimizations = false;

    /**
     * rsIDs from this file are used to populate the ID column of the output. Also, the DB INFO flag will be set when appropriate.
     * dbSNP is not used in any way for the calculations themselves.
     */
    public static DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();
    private double log10GlobalReadMismappingRate;

    /**
     * Active region trimmer reference.
     */
    protected ActiveRegionTrimmer trimmer = new ActiveRegionTrimmer();

    public RodBinding<VariantContext> getDbsnpRodBinding() {
        return dbsnp.dbsnp;
    }

    /**
     * If a call overlaps with a record from the provided comp track, the INFO field will be annotated
     * as such in the output with the track name (e.g. -comp:FOO will have 'FOO' in the INFO field). Records that are
     * filtered in the comp track will be ignored. Note that 'dbSNP' has been special-cased (see the --dbsnp argument).
     */
    public List<RodBinding<VariantContext>> comps = Collections.emptyList();

    public List<RodBinding<VariantContext>> getCompRodBindings() {
        return comps;
    }

    // The following are not used by the Unified Genotyper
    public RodBinding<VariantContext> getSnpEffRodBinding() {
        return null;
    }

    public List<RodBinding<VariantContext>> getResourceRodBindings() {
        return Collections.emptyList();
    }

    public boolean alwaysAppendDbsnpId() {
        return false;
    }

    protected List<String> infoFieldAnnotations = Arrays.asList(new String[]{"BaseQualityRankSumTest", "ChromosomeCounts", "ClippingRankSumTest", "Coverage",
            "FisherStrand", "HaplotypeScore", "InbreedingCoeff", "MappingQualityRankSumTest", "MappingQualityZero", "QualByDepth", "RMSMappingQuality",
            "ReadPosRankSumTest", "SpanningDeletions", "StrandOddsRatio", "TandemRepeatAnnotator"});

    protected List<String> genotypeAnnotations = new ArrayList<>(Arrays.asList(new String[]{"DepthPerAlleleBySample", "DepthPerSampleHC", "StrandBiasBySample"}));

    /**
     * Which annotations to exclude from output in the VCF file.  Note that this argument has higher priority than the -A or -G arguments,
     * so these annotations will be excluded even if they are explicitly included with the other options.
     */
    protected List<String> annotationsToExclude = new ArrayList<>(Arrays.asList(new String[]{"SpanningDeletions", "TandemRepeatAnnotator"}));

    /**
     * You can use this argument to specify that HC should process a single sample out of a multisample BAM file. This
     * is especially useful if your samples are all in the same file but you need to run them individually through HC
     * in -ERC GVC mode (which is the recommended usage). Note that the name is case-sensitive.
     */
    protected String sampleNameToUse = null;

    // -----------------------------------------------------------------------------------------------
    // arguments to control internal behavior of the read threading assembler
    // -----------------------------------------------------------------------------------------------

    /**
     * Multiple kmer sizes can be specified, using e.g. `-kmerSize 10 -kmerSize 25`.
     */
    protected List<Integer> kmerSizes = Arrays.asList(10, 25);

    /**
     * When graph cycles are detected, the normal behavior is to increase kmer sizes iteratively until the cycles are
     * resolved. Disabling this behavior may cause the program to give up on assembling the ActiveRegion.
     */
    protected boolean dontIncreaseKmerSizesForCycles = false;

    /**
     * By default, the program does not allow processing of reference sections that contain non-unique kmers. Disabling
     * this check may cause problems in the assembly graph.
     */
    protected boolean allowNonUniqueKmersInRef = false;

    /**
     * If fewer samples than the specified number pass the minPruning threshold for a given path, that path will be eliminated from the graph.
     */
    protected int numPruningSamples = 1;

    /**
     * As of version 3.3, this argument is no longer needed because dangling end recovery is now the default behavior. See GATK 3.3 release notes for more details.
     */
    @Deprecated
    protected boolean DEPRECATED_RecoverDanglingHeads = false;

    /**
     * By default, the read threading assembler will attempt to recover dangling heads and tails. See the `minDanglingBranchLength` argument documentation for more details.
     */
    protected boolean doNotRecoverDanglingBranches = false;

    /**
     * When constructing the assembly graph we are often left with "dangling" branches.  The assembly engine attempts to rescue these branches
     * by merging them back into the main graph.  This argument describes the minimum length of a dangling branch needed for the engine to
     * try to rescue it.  A smaller number here will lead to higher sensitivity to real variation but also to a higher number of false positives.
     */
    protected int minDanglingBranchLength = 4;

    /**
     * This argument is specifically intended for 1000G consensus analysis mode. Setting this flag will inject all
     * provided alleles to the assembly graph but will not forcibly genotype all of them.
     */
    protected boolean consensusMode = false;

    // -----------------------------------------------------------------------------------------------
    // general advanced arguments to control haplotype caller behavior
    // -----------------------------------------------------------------------------------------------


    /**
     * When HC is run in reference confidence mode with banding compression enabled (-ERC GVCF), homozygous-reference
     * sites are compressed into bands of similar genotype quality (GQ) that are emitted as a single VCF record. See
     * the FAQ documentation for more details about the GVCF format.
     * <p>
     * This argument allows you to set the GQ boundaries. HC expects a list of multiple GQ threshold values. To pass
     * multiple values, you provide them one by one with the argument, as in `-GQB 10 -GQB 20 -GQB 30` and so on. Note
     * that GQ values are capped at 99 in the GATK.
     */
    protected List<Integer> GVCFGQBands = new ArrayList<Integer>(70) {{
        for (int i = 1; i <= 60; ++i) add(i);
        add(70);
        add(80);
        add(90);
        add(99);
    }};

    /**
     * This parameter determines the maximum size of an indel considered as potentially segregating in the
     * reference model.  It is used to eliminate reads from being indel informative at a site, and determines
     * by that mechanism the certainty in the reference base.  Conceptually, setting this parameter to
     * X means that each informative read is consistent with any indel of size < X being present at a specific
     * position in the genome, given its alignment to the reference.
     */
    protected int indelSizeToEliminateInRefModel = 10;

    // -----------------------------------------------------------------------------------------------
    // general advanced arguments to control haplotype caller behavior
    // -----------------------------------------------------------------------------------------------

    /**
     * Bases with a quality below this threshold will not be used for calling.
     */
    public byte MIN_BASE_QUALTY_SCORE = 10;

    /**
     * Paths with fewer supporting kmers than the specified threshold will be pruned from the graph.
     * <p>
     * Be aware that this argument can dramatically affect the results of variant calling and should only be used with great caution.
     * Using a prune factor of 1 (or below) will prevent any pruning from the graph, which is generally not ideal; it can make the
     * calling much slower and even less accurate (because it can prevent effective merging of "tails" in the graph).  Higher values
     * tend to make the calling much faster, but also lowers the sensitivity of the results (because it ultimately requires higher
     * depth to produce calls).
     */
    protected int MIN_PRUNE_FACTOR = 2;

    protected int gcpHMM = 10;

    /**
     * If this flag is provided, the haplotype caller will include unmapped reads (that have chromosomal coordinates) in the assembly and calling
     * when these reads occur in the region being analyzed.  Typically, for paired end analyses, one pair of the
     * read can map, but if its pair is too divergent then it may be unmapped and placed next to its mate, taking
     * the mates contig and alignment start.  If this flag is provided the haplotype caller will see such reads,
     * and may make use of them in assembly and calling, where possible.
     */
    protected boolean includeUnmappedReads = false;

    protected boolean USE_ALLELES_TRIGGER = false;

    /**
     * The phredScaledGlobalReadMismappingRate reflects the average global mismapping rate of all reads, regardless of their
     * mapping quality.  This term effects the probability that a read originated from the reference haplotype, regardless of
     * its edit distance from the reference, in that the read could have originated from the reference haplotype but
     * from another location in the genome.  Suppose a read has many mismatches from the reference, say like 5, but
     * has a very high mapping quality of 60.  Without this parameter, the read would contribute 5 * Q30 evidence
     * in favor of its 5 mismatch haplotype compared to reference, potentially enough to make a call off that single
     * read for all of these events.  With this parameter set to Q30, though, the maximum evidence against any haplotype
     * that this (and any) read could contribute is Q30.
     * <p>
     * Set this term to any negative number to turn off the global mapping rate.
     */
    protected int phredScaledGlobalReadMismappingRate = 45;

    /**
     * The assembly graph can be quite complex, and could imply a very large number of possible haplotypes.  Each haplotype
     * considered requires N PairHMM evaluations if there are N reads across all samples.  In order to control the
     * run of the haplotype caller we only take maxNumHaplotypesInPopulation paths from the graph, in order of their
     * weights, no matter how many paths are possible to generate from the graph.  Putting this number too low
     * will result in dropping true variation because paths that include the real variant are not even considered.
     * You can consider increasing this number when calling organisms with high heterozygosity.
     */
    protected int maxNumHaplotypesInPopulation = 128;

    protected boolean mergeVariantsViaLD = false;

    /**
     * As of GATK 3.3, HaplotypeCaller outputs physical information (see release notes and documentation for details). This argument disables that behavior.
     */
    protected boolean doNotRunPhysicalPhasing = false;

    public static final String HAPLOTYPE_CALLER_PHASING_ID_KEY = "PID";
    public static final String HAPLOTYPE_CALLER_PHASING_GT_KEY = "PGT";

    // -----------------------------------------------------------------------------------------------
    // arguments for debugging / developing the haplotype caller
    // -----------------------------------------------------------------------------------------------
    /**
     * The PairHMM implementation to use for genotype likelihood calculations. The various implementations balance a tradeoff of accuracy and runtime.
     */
    public PairHMM.HMM_IMPLEMENTATION pairHMM = PairHMM.HMM_IMPLEMENTATION.VECTOR_LOGLESS_CACHING;

    protected String keepRG = null;

    /**
     * This argument is intended for benchmarking and scalability testing.
     */
    protected boolean justDetermineActiveRegions = false;

    /**
     * This argument is intended for benchmarking and scalability testing.
     */
    protected boolean dontGenotype = false;

    /**
     * Enabling this argument may cause fundamental problems with the assembly graph itself.
     */
    protected boolean errorCorrectKmers = false;

    protected boolean debugGraphTransformations = false;

    protected boolean dontUseSoftClippedBases = false;

    protected boolean captureAssemblyFailureBAM = false;

    protected boolean allowCyclesInKmerGraphToGeneratePaths = false;

    protected boolean noFpga = false;

    // Parameters to control read error correction
    /**
     * Enabling this argument may cause fundamental problems with the assembly graph itself.
     */
    protected boolean errorCorrectReads = false;

    /**
     * Enabling this argument may cause fundamental problems with the assembly graph itself.
     */
    protected int kmerLengthForReadErrorCorrection = 25;

    protected int minObservationsForKmerToBeSolid = 20;

    /**
     * When calculating the likelihood of variants, we can try to correct for PCR errors that cause indel artifacts.
     * The correction is based on the reference context, and acts specifically around repetitive sequences that tend
     * to cause PCR errors). The variant likelihoods are penalized in increasing scale as the context around a
     * putative indel is more repetitive (e.g. long homopolymer). The correction can be disabling by specifying
     * '-pcrModel NONE'; in that case the default base insertion/deletion qualities will be used (or taken from the
     * read if generated through the BaseRecalibrator). <b>VERY IMPORTANT: when using PCR-free sequencing data we
     * definitely recommend setting this argument to NONE</b>.
     */
    public PairHMMLikelihoodCalculationEngine.PCR_ERROR_MODEL pcrErrorModel = PairHMMLikelihoodCalculationEngine.PCR_ERROR_MODEL.CONSERVATIVE;


    // reference base padding size
    public static final int REFERENCE_PADDING = 500;

    /**
     * When downsampling, level the coverage of the reads in each sample to no more than maxReadsInRegionPerSample reads,
     * not reducing coverage at any read start to less than minReadsPerAlignmentStart
     */
    protected int maxReadsInRegionPerSample = 1000;

    protected int minReadsPerAlignmentStart = 5;

    private byte MIN_TAIL_QUALITY;
    private static final byte MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION = 6;

    // the minimum length of a read we'd consider using for genotyping
    private final static int MIN_READ_LENGTH = 10;

    private HaplotypeCallerArgumentCollection SCAC = new HaplotypeCallerArgumentCollection();
    private SampleList samplesList;

    // the assembly engine
    private LocalAssemblyEngine assemblyEngine = null;

    ReferenceConfidenceModel referenceConfidenceModel = null;

    // the likelihoods engine
    private ReadLikelihoodCalculationEngine likelihoodCalculationEngine = null;

    // the genotyping engine
    private HaplotypeCallerGenotypingEngine genotypingEngine = null;

    private GenomeLocParser genomeLocParser = null;

    public HaplotypeCaller(GenomeLocParser genomeLocParser,
                           List<String> samples,
                           boolean useGVCF) {
        this.genomeLocParser = genomeLocParser;
        setGVCFArg(useGVCF);

        this.samplesList = new IndexedSampleList(samples);
        initialize();
    }

    protected void initialize() {
        if (emitReferenceConfidence()) {

            if (SCAC.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES)
                throw new UserException.BadArgumentValue("ERC/gt_mode", "you cannot request reference confidence output and GENOTYPE_GIVEN_ALLELES at the same time");

            SCAC.genotypeArgs.STANDARD_CONFIDENCE_FOR_EMITTING = -0.0;
            SCAC.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING = -0.0;


            // also, we don't need to output several of the annotations
            annotationsToExclude.add("ChromosomeCounts");
            annotationsToExclude.add("FisherStrand");
            annotationsToExclude.add("StrandOddsRatio");
            annotationsToExclude.add("QualByDepth");

            SCAC.annotateAllSitesWithPLs = true;
        }

        genotypingEngine = new HaplotypeCallerGenotypingEngine(SCAC, samplesList, genomeLocParser, FixedAFCalculatorProvider.createThreadSafeProvider(SCAC, logger), !doNotRunPhysicalPhasing);

        final VariantAnnotatorEngine annotationEngine = new VariantAnnotatorEngine(infoFieldAnnotations, genotypeAnnotations, annotationsToExclude, this, genomeLocParser);

        referenceConfidenceModel = new ReferenceConfidenceModel(genomeLocParser, samplesList, indelSizeToEliminateInRefModel);

        // create and setup the assembler
        assemblyEngine = new ReadThreadingAssembler(maxNumHaplotypesInPopulation, kmerSizes, dontIncreaseKmerSizesForCycles, allowNonUniqueKmersInRef, numPruningSamples);

        assemblyEngine.setErrorCorrectKmers(errorCorrectKmers);
        assemblyEngine.setPruneFactor(MIN_PRUNE_FACTOR);
        assemblyEngine.setDebug(SCAC.DEBUG);
        assemblyEngine.setDebugGraphTransformations(debugGraphTransformations);
        assemblyEngine.setAllowCyclesInKmerGraphToGeneratePaths(allowCyclesInKmerGraphToGeneratePaths);
        assemblyEngine.setRecoverDanglingBranches(!doNotRecoverDanglingBranches);
        assemblyEngine.setMinDanglingBranchLength(minDanglingBranchLength);
        assemblyEngine.setMinBaseQualityToUseInAssembly(MIN_BASE_QUALTY_SCORE);

        MIN_TAIL_QUALITY = (byte) (MIN_BASE_QUALTY_SCORE - 1);
        // setup the likelihood calculation engine
        if (phredScaledGlobalReadMismappingRate < 0) phredScaledGlobalReadMismappingRate = -1;

        // configure the global mismapping rate
        if (phredScaledGlobalReadMismappingRate < 0) {
            log10GlobalReadMismappingRate = -Double.MAX_VALUE;
        } else {
            log10GlobalReadMismappingRate = QualityUtils.qualToErrorProbLog10(phredScaledGlobalReadMismappingRate);
        }

        //static member function - set number of threads
        PairHMM.setNumberOfThreads(1);
        // create our likelihood calculation engine
        likelihoodCalculationEngine = createLikelihoodCalculationEngine();

        final MergeVariantsAcrossHaplotypes variantMerger = mergeVariantsViaLD ? new LDMerger(SCAC.DEBUG, 10, 1) : new MergeVariantsAcrossHaplotypes();

        genotypingEngine.setCrossHaplotypeEventMerger(variantMerger);

        genotypingEngine.setAnnotationEngine(annotationEngine);

        trimmer.initialize(genomeLocParser, SCAC.DEBUG,
                SCAC.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES, emitReferenceConfidence());
    }

    /**
     * Instantiates the appropriate likelihood calculation engine.
     *
     * @return never {@code null}.
     */
    private ReadLikelihoodCalculationEngine createLikelihoodCalculationEngine() {
        switch (likelihoodEngineImplementation) {
            case PairHMM:
                return new PairHMMLikelihoodCalculationEngine((byte) gcpHMM, pairHMM, log10GlobalReadMismappingRate, noFpga, pcrErrorModel);
            default:
                //Note: we do not include in the error message list as it is of no grand public interest.
                throw new UserException("Unsupported likelihood calculation engine '" + likelihoodCalculationEngine +
                        "'. Please use one of the following instead: 'PairHMM' or 'GraphBased'.");
        }
    }

    private void setGVCFArg(boolean useGVCF) {
        if (useGVCF) {
            SCAC.emitReferenceConfidence = ReferenceConfidenceMode.GVCF;
        }
    }

    public List<VariantContext> map(ActiveRegionMapData activeRegionMapData) {
        List<VariantContext> resultVCFRecords = new ArrayList<>();
        VariantContextWriter vcfWriter = null;

        if (SCAC.emitReferenceConfidence == ReferenceConfidenceMode.GVCF) {
            vcfWriter = new GVCFWriter(GVCFGQBands, SCAC.genotypeArgs.samplePloidy);
        }

        List<VariantContext> result = null;
        ActiveRegion originalActiveRegion = activeRegionMapData.activeRegion;
        RefMetaDataTracker metaDataTracker = activeRegionMapData.tracker;
        byte[] fullReferenceWithPadding = activeRegionMapData.fullReferenceWithPadding;
        byte[] refBases = activeRegionMapData.refBases;

        if (!originalActiveRegion.isActive()) {
            // Not active so nothing to do!
            result = referenceModelForNoVariation(originalActiveRegion, true, refBases);
            addAllIntoWriter(result, vcfWriter, resultVCFRecords);
            return resultVCFRecords;
        } else if (originalActiveRegion.size() == 0) {
            result = referenceModelForNoVariation(originalActiveRegion, true, refBases);
            addAllIntoWriter(result, vcfWriter, resultVCFRecords);
            return resultVCFRecords;
        }

        // run the local assembler, getting back a collection of information on how we should proceed
        final List<VariantContext> givenAlleles = new ArrayList<>();
        final AssemblyResultSet untrimmedAssemblyResult = assembleReads(originalActiveRegion, givenAlleles, fullReferenceWithPadding, refBases);

        final TreeSet<VariantContext> allVariationEvents = untrimmedAssemblyResult.getVariationEvents();
        // TODO - line bellow might be unnecessary : it might be that assemblyResult will always have those alleles anyway
        // TODO - so check and remove if that is the case:
        allVariationEvents.addAll(givenAlleles);

        final ActiveRegionTrimmer.Result trimmingResult = trimmer.trim(originalActiveRegion, allVariationEvents);

        if (!trimmingResult.isVariationPresent() && !disableOptimizations) {
            result = referenceModelForNoVariation(originalActiveRegion, false, refBases);
            addAllIntoWriter(result, vcfWriter, resultVCFRecords);
            return resultVCFRecords;
        }

        final AssemblyResultSet assemblyResult =
                trimmingResult.needsTrimming() ? untrimmedAssemblyResult.trimTo(trimmingResult.getCallableRegion()) : untrimmedAssemblyResult;

        final ActiveRegion regionForGenotyping = assemblyResult.getRegionForGenotyping();

        // filter out reads from genotyping which fail mapping quality based criteria
        //TODO - why don't do this before any assembly is done? Why not just once at the beginning of this method
        //TODO - on the originalActiveRegion?
        //TODO - if you move this up you might have to consider to change referenceModelForNoVariation
        //TODO - that does also filter reads.
        final Collection<GATKSAMRecord> filteredReads = filterNonPassingReads(regionForGenotyping);
        final Map<String, List<GATKSAMRecord>> perSampleFilteredReadList = splitReadsBySample(filteredReads);

        // abort early if something is out of the acceptable range
        // TODO is this ever true at this point??? perhaps GGA. Need to check.
        if (!assemblyResult.isVariationPresent() && !disableOptimizations) {
            result = referenceModelForNoVariation(originalActiveRegion, false, refBases);

            addAllIntoWriter(result, vcfWriter, resultVCFRecords);
            return resultVCFRecords;
        }

        // TODO is this ever true at this point??? perhaps GGA. Need to check.
        if (regionForGenotyping.size() == 0 && !disableOptimizations) {
            // no reads remain after filtering so nothing else to do!
            result = referenceModelForNoVariation(originalActiveRegion, false, refBases);

            addAllIntoWriter(result, vcfWriter, resultVCFRecords);
            return resultVCFRecords;
        }

        // evaluate each sample's reads against all haplotypes
        //logger.info("Computing read likelihoods with " + assemblyResult.regionForGenotyping.size() + " reads");
        final List<Haplotype> haplotypes = assemblyResult.getHaplotypeList();
        final Map<String, List<GATKSAMRecord>> reads = splitReadsBySample(regionForGenotyping.getReads());

        // Calculate the likelihoods: CPU intensive part.
        final ReadLikelihoods<Haplotype> readLikelihoods =
                likelihoodCalculationEngine.computeReadLikelihoods(assemblyResult, samplesList, reads);

        // Realign reads to their best haplotype.
        final Map<GATKSAMRecord, GATKSAMRecord> readRealignments = realignReadsToTheirBestHaplotype(readLikelihoods, assemblyResult.getPaddedReferenceLoc());
        readLikelihoods.changeReads(readRealignments);

        // Note: we used to subset down at this point to only the "best" haplotypes in all samples for genotyping, but there
        //  was a bad interaction between that selection and the marginalization that happens over each event when computing
        //  GLs.  In particular, for samples that are heterozygous non-reference (B/C) the marginalization for B treats the
        //  haplotype containing C as reference (and vice versa).  Now this is fine if all possible haplotypes are included
        //  in the genotyping, but we lose information if we select down to a few haplotypes.  [EB]

        final HaplotypeCallerGenotypingEngine.CalledHaplotypes calledHaplotypes = genotypingEngine.assignGenotypeLikelihoods(
                haplotypes,
                readLikelihoods,
                perSampleFilteredReadList,
                assemblyResult.getFullReferenceWithPadding(),
                assemblyResult.getPaddedReferenceLoc(),
                regionForGenotyping.getLocation(),
                genomeLocParser,
                metaDataTracker,
                (consensusMode ? Collections.<VariantContext>emptyList() : givenAlleles),
                emitReferenceConfidence());

        if (emitReferenceConfidence()) {
            if (!containsCalls(calledHaplotypes)) {
                // no called all of the potential haplotypes
                result = referenceModelForNoVariation(originalActiveRegion, false, refBases);
            } else {
                result = new ArrayList<>();
                // output left-flanking non-variant section:
                if (trimmingResult.hasLeftFlankingRegion()) {
                    ActiveRegion leftRegion = trimmingResult.nonVariantLeftFlankRegion();
                    byte[] trimedRef = getTrimedRefBases(originalActiveRegion, leftRegion, refBases);
                    result.addAll(referenceModelForNoVariation(leftRegion, false, trimedRef));
                }
                // output variant containing region.
                result.addAll(referenceConfidenceModel.calculateRefConfidence(assemblyResult.getReferenceHaplotype(),
                        calledHaplotypes.getCalledHaplotypes(), assemblyResult.getPaddedReferenceLoc(), regionForGenotyping,
                        readLikelihoods, genotypingEngine.getPloidyModel(), genotypingEngine.getGenotypingModel(), calledHaplotypes.getCalls()));
                // output right-flanking non-variant section:
                if (trimmingResult.hasRightFlankingRegion()) {
                    ActiveRegion rightRegion = trimmingResult.nonVariantRightFlankRegion();
                    byte[] trimedRef = getTrimedRefBases(originalActiveRegion, rightRegion, refBases);
                    result.addAll(referenceModelForNoVariation(rightRegion, false, trimedRef));
                }
            }
        } else {
            result = calledHaplotypes.getCalls();
        }

        addAllIntoWriter(result, vcfWriter, resultVCFRecords);

        if (SCAC.emitReferenceConfidence == ReferenceConfidenceMode.GVCF) {
            ((GVCFWriter) vcfWriter).close(false);
            return ((GVCFWriter) vcfWriter).getResultGVCFWriter();
        } else {
            return resultVCFRecords;
        }
    }

    private byte[] getTrimedRefBases(ActiveRegion originalRegion, ActiveRegion trimedRegion, byte[] originalRefBases) {
        GenomeLoc originalLoc = originalRegion.getExtendedLoc();
        GenomeLoc trimedLoc = trimedRegion.getExtendedLoc();
        if(originalLoc.getStart() > trimedLoc.getStart() || originalLoc.getStop() < trimedLoc.getStop()) {
            throw new PipelineException("Error !! implement error when get trimed ref base");
        }
        int startIdx = trimedLoc.getStart() - originalLoc.getStart();
        int endIdx = originalRefBases.length + trimedLoc.getStop() - originalLoc.getStop();
        return ArrayUtils.subarray(originalRefBases, startIdx, endIdx);
    }

    /**
     * Returns a map with the original read as a key and the realigned read as the value.
     * <p>
     * Missing keys or equivalent key and value pairs mean that the read was not realigned.
     * </p>
     *
     * @return never {@code null}
     */
    private Map<GATKSAMRecord, GATKSAMRecord> realignReadsToTheirBestHaplotype(final ReadLikelihoods<Haplotype> originalReadLikelihoods, final GenomeLoc paddedReferenceLoc) {

        final Collection<ReadLikelihoods<Haplotype>.BestAllele> bestAlleles = originalReadLikelihoods.bestAlleles();
        final Map<GATKSAMRecord, GATKSAMRecord> result = new HashMap<>(bestAlleles.size());

        for (final ReadLikelihoods<Haplotype>.BestAllele bestAllele : bestAlleles) {
            final GATKSAMRecord originalRead = bestAllele.read;
            final Haplotype bestHaplotype = bestAllele.allele;
            final boolean isInformative = bestAllele.isInformative();
            final GATKSAMRecord realignedRead = AlignmentUtils.createReadAlignedToRef(originalRead, bestHaplotype, paddedReferenceLoc.getStart(), isInformative);
            result.put(originalRead, realignedRead);
        }
        return result;
    }

    private boolean containsCalls(final HaplotypeCallerGenotypingEngine.CalledHaplotypes calledHaplotypes) {
        final List<VariantContext> calls = calledHaplotypes.getCalls();
        if (calls.isEmpty()) return false;
        for (final VariantContext call : calls)
            for (final Genotype genotype : call.getGenotypes())
                if (genotype.isCalled())
                    return true;
        return false;
    }

    /**
     * Create an ref model result (ref model or no calls depending on mode) for an active region without any variation
     * (not is active, or assembled to just ref)
     *
     * @param region             the region to return a no-variation result
     * @param needsToBeFinalized should the region be finalized before computing the ref model (should be false if already done)
     * @return a list of variant contexts (can be empty) to emit for this ref region
     */
    private List<VariantContext> referenceModelForNoVariation(final ActiveRegion region, final boolean needsToBeFinalized, byte[] refBases) {
        if (emitReferenceConfidence()) {
            //TODO - why the activeRegion cannot manage its own one-time finalization and filtering?
            //TODO - perhaps we can remove the last parameter of this method and the three lines bellow?
            if (needsToBeFinalized)
                finalizeActiveRegion(region);
            filterNonPassingReads(region);

            final GenomeLoc paddedLoc = region.getExtendedLoc();
            final Haplotype refHaplotype = createReferenceHaplotype(region, paddedLoc, refBases);
            final List<Haplotype> haplotypes = Collections.singletonList(refHaplotype);
            return referenceConfidenceModel.calculateRefConfidence(refHaplotype, haplotypes,
                    paddedLoc, region, createDummyStratifiedReadMap(refHaplotype, samplesList, region),
                    genotypingEngine.getPloidyModel(), genotypingEngine.getGenotypingModel(), Collections.<VariantContext>emptyList());
        } else
            return Collections.emptyList();
    }

    /**
     * High-level function that runs the assembler on the active region reads,
     * returning a data structure with the resulting information needed
     * for further HC steps
     *
     * @param activeRegion the region we should assemble
     * @param giveAlleles  additional alleles we might need to genotype (can be empty)
     * @return the AssemblyResult describing how to proceed with genotyping
     */
    protected AssemblyResultSet assembleReads(final ActiveRegion activeRegion, final List<VariantContext> giveAlleles, byte[] fullReferenceWithPadding, byte[] refBases) {
        // Create the reference haplotype which is the bases from the reference that make up the active region
        finalizeActiveRegion(activeRegion); // handle overlapping fragments, clip adapter and low qual tails

        final GenomeLoc paddedReferenceLoc = getPaddedLoc(activeRegion);
        final Haplotype referenceHaplotype = createReferenceHaplotype(activeRegion, paddedReferenceLoc, refBases);

        // Create ReadErrorCorrector object if requested - will be used within assembly engine.
        ReadErrorCorrector readErrorCorrector = null;

        try {
            final AssemblyResultSet assemblyResultSet = assemblyEngine.runLocalAssembly(activeRegion, referenceHaplotype, fullReferenceWithPadding, paddedReferenceLoc, giveAlleles, readErrorCorrector);
            assemblyResultSet.debugDump(logger);
            return assemblyResultSet;

        } catch (final Exception e) {
            throw e;
        }
    }

    /**
     * Helper function to create the reference haplotype out of the active region and a padded loc
     *
     * @param activeRegion       the active region from which to generate the reference haplotype
     * @param paddedReferenceLoc the GenomeLoc which includes padding and shows how big the reference haplotype should be
     * @return a non-null haplotype
     */
    private Haplotype createReferenceHaplotype(final ActiveRegion activeRegion, final GenomeLoc paddedReferenceLoc, byte[] refBases) {
        return ReferenceConfidenceModel.createReferenceHaplotype(activeRegion, refBases, paddedReferenceLoc);
    }

    /**
     * Create a context that maps each read to the reference haplotype with log10 L of 0
     *
     * @param refHaplotype a non-null reference haplotype
     * @param samples      a list of all samples
     * @param region       the active region containing reads
     * @return a map from sample -> PerReadAlleleLikelihoodMap that maps each read to ref
     */
    public static ReadLikelihoods<Haplotype> createDummyStratifiedReadMap(final Haplotype refHaplotype,
                                                                          final SampleList samples,
                                                                          final ActiveRegion region) {
        return new ReadLikelihoods<>(samples, new IndexedAlleleList<>(refHaplotype),
                splitReadsBySample(samples, region.getReads()), false);
    }


    private void finalizeActiveRegion(final ActiveRegion activeRegion) {
        if (activeRegion.isFinalized()) return;

        // Loop through the reads hard clipping the adaptor and low quality tails
        final List<GATKSAMRecord> readsToUse = new ArrayList<>(activeRegion.getReads().size());
        for (final GATKSAMRecord myRead : activeRegion.getReads()) {
            GATKSAMRecord clippedRead;
            if (errorCorrectReads)
                clippedRead = ReadClipper.hardClipLowQualEnds(myRead, MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION);
            else  // default case: clip low qual ends of reads
                clippedRead = ReadClipper.hardClipLowQualEnds(myRead, MIN_TAIL_QUALITY);

            if (dontUseSoftClippedBases || !ReadUtils.hasWellDefinedFragmentSize(clippedRead)) {
                // remove soft clips if we cannot reliably clip off adapter sequence or if the user doesn't want to use soft clips at all
                clippedRead = ReadClipper.hardClipSoftClippedBases(clippedRead);
            } else {
                // revert soft clips so that we see the alignment start and end assuming the soft clips are all matches
                // TODO -- WARNING -- still possibility that unclipping the soft clips will introduce bases that aren't
                // TODO -- truly in the extended region, as the unclipped bases might actually include a deletion
                // TODO -- w.r.t. the reference.  What really needs to happen is that kmers that occur before the
                // TODO -- reference haplotype start must be removed
                clippedRead = ReadClipper.revertSoftClippedBases(clippedRead);
            }

            clippedRead = (clippedRead.getReadUnmappedFlag() ? clippedRead : ReadClipper.hardClipAdaptorSequence(clippedRead));
            if (!clippedRead.isEmpty() && clippedRead.getCigar().getReadLength() > 0) {
                clippedRead = ReadClipper.hardClipToRegion(clippedRead, activeRegion.getExtendedLoc().getStart(), activeRegion.getExtendedLoc().getStop());
                if (activeRegion.readOverlapsRegion(clippedRead) && clippedRead.getReadLength() > 0) {
                    //logger.info("Keeping read " + clippedRead + " start " + clippedRead.getAlignmentStart() + " end " + clippedRead.getAlignmentEnd());
                    readsToUse.add(clippedRead);
                }
            }
        }

        // TODO -- Performance optimization: we partition the reads by sample 4 times right now; let's unify that code.

        final List<GATKSAMRecord> downsampledReads = DownsamplingUtils.levelCoverageByPosition(ReadUtils.sortReadsByCoordinate(readsToUse), maxReadsInRegionPerSample, minReadsPerAlignmentStart);

        // handle overlapping read pairs from the same fragment
        cleanOverlappingReadPairs(downsampledReads);

        activeRegion.clearReads();
        activeRegion.addAll(downsampledReads);
        activeRegion.setFinalized(true);
    }

    private Set<GATKSAMRecord> filterNonPassingReads(final ActiveRegion activeRegion) {
        final Set<GATKSAMRecord> readsToRemove = new LinkedHashSet<>();
        for (final GATKSAMRecord rec : activeRegion.getReads()) {
            if (rec.getReadLength() < MIN_READ_LENGTH || rec.getMappingQuality() < 20 || BadMateFilter.hasBadMate(rec) || (keepRG != null && !rec.getReadGroup().getId().equals(keepRG))) {
                readsToRemove.add(rec);
            }
        }
        activeRegion.removeAll(readsToRemove);
        return readsToRemove;
    }

    private GenomeLoc getPaddedLoc(final ActiveRegion activeRegion) {
        final int padLeft = Math.max(activeRegion.getExtendedLoc().getStart() - REFERENCE_PADDING, 1);
        final int padRight = Math.min(activeRegion.getExtendedLoc().getStop() + REFERENCE_PADDING,
                genomeLocParser.getContigInfo(activeRegion.getExtendedLoc().getContig()).getSequenceLength());
        return genomeLocParser.createGenomeLoc(activeRegion.getExtendedLoc().getContig(), padLeft, padRight);
    }

    private Map<String, List<GATKSAMRecord>> splitReadsBySample(final Collection<GATKSAMRecord> reads) {
        return splitReadsBySample(samplesList, reads);
    }

    private static Map<String, List<GATKSAMRecord>> splitReadsBySample(final SampleList samplesList, final Collection<GATKSAMRecord> reads) {
        final Map<String, List<GATKSAMRecord>> returnMap = new HashMap<>();
        final int sampleCount = samplesList.sampleCount();
        for (int i = 0; i < sampleCount; i++)
            returnMap.put(samplesList.sampleAt(i), new ArrayList<>());

        for (final GATKSAMRecord read : reads)
            returnMap.get(read.getReadGroup().getSample()).add(read);

        return returnMap;
    }

    /**
     * Are we emitting a reference confidence in some form, or not?
     *
     * @return true if HC must emit reference confidence.
     */
    private boolean emitReferenceConfidence() {
        return SCAC.emitReferenceConfidence != ReferenceConfidenceMode.NONE;
    }

    /**
     * Clean up reads/bases that overlap within read pairs
     *
     * @param reads the list of reads to consider
     */
    private void cleanOverlappingReadPairs(final List<GATKSAMRecord> reads) {
        for (final List<GATKSAMRecord> perSampleReadList : splitReadsBySample(reads).values()) {
            final FragmentCollection<GATKSAMRecord> fragmentCollection = FragmentUtils.create(perSampleReadList);
            for (final List<GATKSAMRecord> overlappingPair : fragmentCollection.getOverlappingPairs())
                FragmentUtils.adjustQualsOfOverlappingPairedFragments(overlappingPair);
        }
    }

    private void addAllIntoWriter(List<VariantContext> variantContexts,
                                  VariantContextWriter vcfWriter,
                                  List<VariantContext> resultVCFRecords) {
        if (SCAC.emitReferenceConfidence == ReferenceConfidenceMode.GVCF) {
            for (VariantContext context : variantContexts) {
                vcfWriter.add(context);
            }
        } else {
            resultVCFRecords.addAll(variantContexts);
        }
    }
}
