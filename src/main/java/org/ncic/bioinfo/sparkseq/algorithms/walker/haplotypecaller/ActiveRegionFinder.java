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
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.DuplicateReadFilter;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.FailsVendorQualityCheckFilter;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.FilterUtils;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.NotPrimaryAlignmentFilter;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.UnmappedReadFilter;
import org.ncic.bioinfo.sparkseq.algorithms.engine.LocusWalker;
import org.ncic.bioinfo.sparkseq.algorithms.utils.AlignmentContextUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GATKVariantContextUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLoc;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLocParser;
import org.ncic.bioinfo.sparkseq.algorithms.utils.MathUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.QualityUtils;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefMetaDataTracker;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.ReferenceContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.AlignmentContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.GATKSAMRecord;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.RegionSamTraverser;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.SamContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.Filter;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.HCMappingQualityFilter;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.MappingQualityUnavailableFilter;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RODContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RefMetaTrackerTraverser;
import org.ncic.bioinfo.sparkseq.algorithms.walker.baserecalibrator.BAQ;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.afcalculate.FixedAFCalculatorProvider;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.argcollection.HaplotypeCallerArgumentCollection;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.argcollection.UnifiedArgumentCollection;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.genotyper.HaplotypeCallerGenotypingEngine;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.genotyper.UnifiedGenotypingEngine;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.genotyper.VariantCallContext;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.model.GenotypeLikelihoodsCalculationModel;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.model.GenotypingModel;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.model.ReferenceConfidenceMode;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.model.ReferenceConfidenceModel;
import org.ncic.bioinfo.sparkseq.exceptions.UserException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Map;

/**
 * Author: wbc
 */
public class ActiveRegionFinder extends LocusWalker {

    private HaplotypeCallerArgumentCollection SCAC = new HaplotypeCallerArgumentCollection();

    // -----------------------------------------------------------------------------------------------
    // general advanced arguments to control haplotype caller behavior
    // -----------------------------------------------------------------------------------------------

    /**
     * Bases with a quality below this threshold will not be used for calling.
     */
    public byte MIN_BASE_QUALTY_SCORE = 10;

    /**
     * As of GATK 3.3, HaplotypeCaller outputs physical information (see release notes and documentation for details). This argument disables that behavior.
     */
    protected boolean doNotRunPhysicalPhasing = false;

    public static final String HAPLOTYPE_CALLER_PHASING_ID_KEY = "PID";
    public static final String HAPLOTYPE_CALLER_PHASING_GT_KEY = "PGT";

    // -----------------------------------------------------------------------------------------------
    // done with Haplotype caller parameters
    // -----------------------------------------------------------------------------------------------

    // the UG engines
    private UnifiedGenotypingEngine activeRegionEvaluationGenotyperEngine = null;

    // the genotyping engine
    private HaplotypeCallerGenotypingEngine genotypingEngine = null;

    private SampleList samplesList;

    ReferenceConfidenceModel referenceConfidenceModel = null;

    private final static Allele FAKE_REF_ALLELE = Allele.create("N", true); // used in isActive function to call into UG Engine. Should never appear anywhere in a VCF file
    private final static Allele FAKE_ALT_ALLELE = Allele.create("<FAKE_ALT>", false); // used in isActive function to call into UG Engine. Should never appear anywhere in a VCF file

    /////////////////////////////////////////
    //  与active region traverse相关变量    //
    /////////////////////////////////////////
    private boolean walkerHasPresetRegions = false;
    private int activeRegionExtension = -1;
    private int maxRegionSize = -1;
    private int minRegionSize = -1;

    public Integer maxProbPropagationDistance = 50;
    public Double activeProbThreshold = 0.002;
    protected int indelSizeToEliminateInRefModel = 10;

    private List<ActiveRegionMapData> resultActiveRegions = new ArrayList<>();
    private ActivityProfile activityProfile = null;

    public ActiveRegionFinder(GenomeLocParser genomeLocParser,
                              RefContentProvider refContentProvider,
                              SamContentProvider samContentProvider,
                              List<RODContentProvider> rodContentProviderList,
                              boolean useGVCF) {
        super(genomeLocParser, refContentProvider, samContentProvider, rodContentProviderList);
        setGVCFArg(useGVCF);
        initialize();
    }

    private void setGVCFArg(boolean useGVCF) {
        if (useGVCF) {
            SCAC.emitReferenceConfidence = ReferenceConfidenceMode.GVCF;
        }
    }

    @Override
    protected void initialize() {
        if (emitReferenceConfidence()) {

            if (SCAC.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES)
                throw new UserException.BadArgumentValue("ERC/gt_mode", "you cannot request reference confidence output and GENOTYPE_GIVEN_ALLELES at the same time");

            SCAC.genotypeArgs.STANDARD_CONFIDENCE_FOR_EMITTING = -0.0;
            SCAC.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING = -0.0;

            SCAC.annotateAllSitesWithPLs = true;
        }

        samplesList = getSamplesListFromContentProvider();

        final UnifiedArgumentCollection simpleUAC = SCAC.cloneTo(UnifiedArgumentCollection.class);
        simpleUAC.outputMode = OutputMode.EMIT_VARIANTS_ONLY;
        simpleUAC.genotypingOutputMode = GenotypingOutputMode.DISCOVERY;
        simpleUAC.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING = Math.min(4.0, SCAC.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING); // low values used for isActive determination only, default/user-specified values used for actual calling
        simpleUAC.genotypeArgs.STANDARD_CONFIDENCE_FOR_EMITTING = Math.min(4.0, SCAC.genotypeArgs.STANDARD_CONFIDENCE_FOR_EMITTING); // low values used for isActive determination only, default/user-specified values used for actual calling
        simpleUAC.CONTAMINATION_FRACTION = 0.0;
        simpleUAC.CONTAMINATION_FRACTION_FILE = null;
        simpleUAC.exactCallsLog = null;
        // Seems that at least with some test data we can lose genuine haploid variation if we use
        // UGs engine with ploidy == 1
        simpleUAC.genotypeArgs.samplePloidy = Math.max(2, SCAC.genotypeArgs.samplePloidy);

        activeRegionEvaluationGenotyperEngine = new UnifiedGenotypingEngine(simpleUAC,
                FixedAFCalculatorProvider.createThreadSafeProvider(simpleUAC, logger),
                samplesList, genomeLocParser, BAQ.CalculationMode.OFF);
        activeRegionEvaluationGenotyperEngine.setLogger(logger);

        genotypingEngine = new HaplotypeCallerGenotypingEngine(SCAC, samplesList, genomeLocParser,
                FixedAFCalculatorProvider.createThreadSafeProvider(SCAC, logger), !doNotRunPhysicalPhasing);

        referenceConfidenceModel = new ReferenceConfidenceModel(genomeLocParser, samplesList, indelSizeToEliminateInRefModel);

        // 关于active region遍历识别的initialize
        walkerHasPresetRegions = false;
        activeRegionExtension = 100;
        maxRegionSize = 300;
        minRegionSize = 50;
        double bandPassSigma = 17.0;
        activityProfile = new BandPassActivityProfile(genomeLocParser, null, maxProbPropagationDistance, activeProbThreshold,
                BandPassActivityProfile.MAX_FILTER_SIZE, bandPassSigma);
    }

    private SampleList getSamplesListFromContentProvider() {
        List<SAMReadGroupRecord> readGroupInfos =
                samContentProvider.getSamFileHeader().getReadGroups();
        List<String> samples = new ArrayList<>();
        for (SAMReadGroupRecord readGroup : readGroupInfos) {
            samples.add(readGroup.getSample());
        }
        SampleList sampleList = new IndexedSampleList(samples);
        return sampleList;
    }

    private boolean emitReferenceConfidence() {
        return SCAC.emitReferenceConfidence != ReferenceConfidenceMode.NONE;
    }

    /**
     * 对应HaplotypeCaller中的isActive
     *
     * @param tracker
     * @param ref
     * @param context
     */
    @Override
    public void map(RefMetaDataTracker tracker,
                    ReferenceContext ref,
                    AlignmentContext context) {

        ActivityProfileState result = null;

        if (context == null || context.getBasePileup().isEmpty()) {
            // if we don't have any data, just abort early
            result = new ActivityProfileState(ref.getLocus(), 0.0);
        } else {

            final int ploidy = activeRegionEvaluationGenotyperEngine.getConfiguration().genotypeArgs.samplePloidy;
            final List<Allele> noCall = GATKVariantContextUtils.noCallAlleles(ploidy); // used to noCall all genotypes until the exact model is applied
            final Map<String, AlignmentContext> splitContexts = AlignmentContextUtils.splitContextBySampleName(context);
            final GenotypesContext genotypes = GenotypesContext.create(splitContexts.keySet().size());
            final MathUtils.RunningAverage averageHQSoftClips = new MathUtils.RunningAverage();
            final GenotypingModel genotypingModel = genotypingEngine.getGenotypingModel();
            for (final Map.Entry<String, AlignmentContext> sample : splitContexts.entrySet()) {
                final String sampleName = sample.getKey();
                // The ploidy here is not dictated by the sample but by the simple genotyping-engine used to determine whether regions are active or not.
                final int activeRegionDetectionHackishSamplePloidy = activeRegionEvaluationGenotyperEngine.getConfiguration().genotypeArgs.samplePloidy;
                final double[] genotypeLikelihoods = referenceConfidenceModel.calcGenotypeLikelihoodsOfRefVsAny(sampleName, activeRegionDetectionHackishSamplePloidy, genotypingModel, sample.getValue().getBasePileup(), ref.getBase(), MIN_BASE_QUALTY_SCORE, averageHQSoftClips).genotypeLikelihoods;
                genotypes.add(new GenotypeBuilder(sample.getKey()).alleles(noCall).PL(genotypeLikelihoods).make());
            }

            final List<Allele> alleles = Arrays.asList(FAKE_REF_ALLELE, FAKE_ALT_ALLELE);
            final VariantCallContext vcOut = activeRegionEvaluationGenotyperEngine.calculateGenotypes(new VariantContextBuilder("HCisActive!", context.getContig(), context.getLocation().getStart(), context.getLocation().getStop(), alleles).genotypes(genotypes).make(), GenotypeLikelihoodsCalculationModel.Model.SNP);
            final double isActiveProb = vcOut == null ? 0.0 : QualityUtils.qualToProb(vcOut.getPhredScaledQual());

            result = new ActivityProfileState(ref.getLocus(), isActiveProb, averageHQSoftClips.mean() > 6.0 ? ActivityProfileState.Type.HIGH_QUALITY_SOFT_CLIPS : ActivityProfileState.Type.NONE, averageHQSoftClips.mean());
        }

        activityProfile.add(result);
    }

    @Override
    protected List<Filter> getFilter() {
        List<Filter> filters = new ArrayList<>();
        filters.add(new HCMappingQualityFilter());
        filters.add(new MappingQualityUnavailableFilter());
        return filters;
    }

    @Override
    protected void onTraversalDone() {

        FilterUtils filterUtils = new FilterUtils();
        filterUtils.addFilter(new UnmappedReadFilter());
        filterUtils.addFilter(new NotPrimaryAlignmentFilter());
        filterUtils.addFilter(new DuplicateReadFilter());
        filterUtils.addFilter(new FailsVendorQualityCheckFilter());
        filterUtils.addFilter(new HCMappingQualityFilter());
        filterUtils.addFilter(new MappingQualityUnavailableFilter());

        RefMetaTrackerTraverser rodMetaTrackerTraverser = new RefMetaTrackerTraverser(rodContentProviderList);
        RegionSamTraverser regionSamTraverser = new RegionSamTraverser(this.samContentProvider, filterUtils);
        prepActiveRegionsForProcessing(regionSamTraverser, rodMetaTrackerTraverser);
    }

    /**
     * Take the individual isActive calls and integrate them into contiguous active regions and
     * add these blocks of work to the work queue
     * band-pass filter the list of isActive probabilities and turn into active regions
     * <p>
     * forceAllRegionsToBeActive 强制都是active的，以前的只是为了避免他们提前pop出来，现在遍历完了再走的，没关系
     *
     * @param regionSamTraverser
     * @param rodMetaTrackerTraverser
     */
    private void prepActiveRegionsForProcessing(final RegionSamTraverser regionSamTraverser,
                                                final RefMetaTrackerTraverser rodMetaTrackerTraverser) {
        // We don't have preset regions, so we get our regions from the activity profile
        final Collection<ActiveRegion> activeRegions = activityProfile
                .popReadyActiveRegions(activeRegionExtension, minRegionSize, maxRegionSize, true);

        for (ActiveRegion activeRegion : activeRegions) {
            resultActiveRegions.add(prepActiveRegionForProcessing(
                    activeRegion, regionSamTraverser, rodMetaTrackerTraverser));
        }
    }

    private ActiveRegionMapData prepActiveRegionForProcessing(final ActiveRegion activeRegion,
                                                              final RegionSamTraverser regionSamTraverser,
                                                              final RefMetaTrackerTraverser rodMetaTrackerTraverser) {
        GenomeLoc extendLoc = activeRegion.getExtendedLoc();
        List<GATKSAMRecord> overlappedRecords = regionSamTraverser.getOverlappedReads(extendLoc);
        for (GATKSAMRecord record : overlappedRecords) {
            activeRegion.add(record);
        }
        // prepare the RefMetaDataTracker information
        final GenomeLoc loc = activeRegion.getLocation();
        // get all of the RODs that cover the active region (without extension)
        final RefMetaDataTracker tracker = rodMetaTrackerTraverser.getOverlappedTracker(loc);

        return new ActiveRegionMapData(activeRegion, tracker, refContentProvider);
    }

    public List<ActiveRegionMapData> getResultActiveRegions() {
        return resultActiveRegions;
    }
}
