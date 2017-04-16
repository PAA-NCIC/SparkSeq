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
package org.ncic.bioinfo.sparkseq.algorithms.walker.mutect;

import htsjdk.samtools.SAMRecord;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.collections.ListUtils;
import org.apache.commons.lang.math.NumberUtils;
import org.apache.commons.math.MathException;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefMetaDataTracker;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.ReferenceContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.AlignmentContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.GATKSAMRecord;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.IntervalLocusSamTraverser;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.PileupElement;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.ReadBackedPileup;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.ReadBackedPileupImpl;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.SamContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.Filter;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.FilterUtils;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.GATKFeature;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RODContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RODNames;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RodBinding;
import org.ncic.bioinfo.sparkseq.algorithms.engine.IntervalLocusWalker;
import org.ncic.bioinfo.sparkseq.algorithms.utils.CGAAlignmentUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLoc;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLocParser;
import org.ncic.bioinfo.sparkseq.exceptions.GATKException;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.TreeMap;

/**
 * Author: wbc
 */
public class Mutect extends IntervalLocusWalker {

    public MuTectArgumentCollection MTAC = new MuTectArgumentCollection();

    public int MIN_QSUM_QSCORE = 13;
    public boolean USE_MAPQ0_IN_NORMAL_QSCORE = true;

    private boolean hasNormalBam = false;

    private double contaminantAlternateFraction;

    private TumorPowerCalculator tumorPowerCalculator;
    private NormalPowerCalculator normalNovelSitePowerCalculator;
    private NormalPowerCalculator normalDbSNPSitePowerCalculator;
    private TumorPowerCalculator strandArtifactPowerCalculator;

    private CallStatsGenerator callStatsGenerator;

    public List<RodBinding<VariantContext>> normalPanelRod = Collections.emptyList();

    private IntervalLocusSamTraverser normalSamTraverser = null;

    private List<VariantContext> resultVCFRecords = new ArrayList<>();
    private List<String> resultVCFOutInfos = new ArrayList<>();

    public Mutect(GenomeLocParser genomeLocParser,
                  RefContentProvider refContentProvider,
                  SamContentProvider tumorContentProvider,
                  SamContentProvider normalContentProvider,
                  List<RODContentProvider> rodContentProviderList,
                  List<GenomeLoc> intervals) {
        super(genomeLocParser, refContentProvider, tumorContentProvider, rodContentProviderList, intervals);
        if (normalContentProvider != null) {
            hasNormalBam = true;
            normalSamTraverser = getNormalLocusSamTraverser(normalContentProvider, intervals);
        } else {
            normalSamTraverser = getNormalLocusSamTraverser(
                    new SamContentProvider(ListUtils.EMPTY_LIST, tumorContentProvider.getSamFileHeader()), intervals);
        }

        initialize();
    }

    private IntervalLocusSamTraverser getNormalLocusSamTraverser(
            SamContentProvider normalContentProvider, List<GenomeLoc> intervals) {
        FilterUtils filterUtils = getLocusWalkerFilterUtils();
        List<Filter> filters = getFilter();
        if (filters != null) {
            for (Filter filter : filters) {
                filterUtils.addFilter(filter);
            }
        }

        GenomeLoc allLocus = refContentProvider.getLocus();
        IntervalLocusSamTraverser traverser = new IntervalLocusSamTraverser(normalContentProvider, allLocus, intervals, filterUtils);
        traverser.rewind();
        return traverser;
    }

    @Override
    protected void initialize() {
        callStatsGenerator = new CallStatsGenerator(MTAC.ENABLE_QSCORE_OUTPUT);

        if (!hasNormalBam) {
            MTAC.NORMAL_LOD_THRESHOLD = -1 * Float.MAX_VALUE;
            MTAC.NORMAL_DBSNP_LOD_THRESHOLD = -1 * Float.MAX_VALUE;
            MTAC.NORMAL_ARTIFACT_LOD_THRESHOLD = Float.MAX_VALUE;
            MTAC.NORMAL_SAMPLE_NAME = "none";
        }

        this.contaminantAlternateFraction = Math.max(MTAC.MINIMUM_MUTATION_CELL_FRACTION, MTAC.FRACTION_CONTAMINATION);

        // coverage related initialization
        double powerConstantEps = Math.pow(10, -1 * (MTAC.POWER_CONSTANT_QSCORE / 10));

        this.tumorPowerCalculator = new TumorPowerCalculator(powerConstantEps, MTAC.TUMOR_LOD_THRESHOLD, this.contaminantAlternateFraction);
        this.normalNovelSitePowerCalculator = new NormalPowerCalculator(powerConstantEps, MTAC.NORMAL_LOD_THRESHOLD);
        this.normalDbSNPSitePowerCalculator = new NormalPowerCalculator(powerConstantEps, MTAC.NORMAL_DBSNP_LOD_THRESHOLD);
        this.strandArtifactPowerCalculator = new TumorPowerCalculator(powerConstantEps, MTAC.STRAND_ARTIFACT_LOD_THRESHOLD, 0.0f);
    }

    @Override
    protected List<Filter> getFilter() {
        List<Filter> filters = new ArrayList<>();
        return filters;
    }

    private List<VariantContext> getVCInTrackerInLocus(String rodName, final RefMetaDataTracker tracker) {
        List<Feature> features = tracker.getValues(rodName);
        List<VariantContext> vcs = new ArrayList<>(features.size());
        for (Feature feature : features) {
            vcs.add((VariantContext) ((GATKFeature.TribbleGATKFeature) feature).getUnderlyingObject());
        }
        return vcs;
    }

    private ReadBackedPileup cleanNoneRefPileupElement(ReadBackedPileup rawPileup) {
        ArrayList<PileupElement> pileupElements = new ArrayList<PileupElement>(rawPileup.depthOfCoverage());

        for (PileupElement p : rawPileup ) {
            final byte base = p.getBase();
            if (base != 'N') {
                pileupElements.add(p);
            }
        }
        return new ReadBackedPileupImpl(rawPileup.getLocation(), pileupElements);
    }

    @Override
    protected void map(RefMetaDataTracker tracker,
                       ReferenceContext ref,
                       AlignmentContext rawContext) {

        final char upRef = Character.toUpperCase((char) ref.getBase());
        if (upRef != 'A' && upRef != 'C' && upRef != 'G' && upRef != 'T') {
            return;
        }

        ReadBackedPileup tumorPileup = cleanNoneRefPileupElement(rawContext.getBasePileup());
        ReadBackedPileup normalPileup = cleanNoneRefPileupElement(normalSamTraverser.next().getBasePileup());
        // an optimization to speed things up when there is no coverage
        if (tumorPileup.depthOfCoverage() == 0 && normalPileup.depthOfCoverage() == 0) {
            return;
        }

        TreeMap<Double, CandidateMutation> messageByTumorLod = new TreeMap<Double, CandidateMutation>();
        // get sequence context around mutation
        String sequenceContext = SequenceUtils.createSequenceContext(this.refContentProvider, ref, 3);

        try {
            final LocusReadPile tumorReadPile = new LocusReadPile(tumorPileup, upRef, MTAC.MIN_QSCORE, MIN_QSUM_QSCORE, false, MTAC.ARTIFACT_DETECTION_MODE, MTAC.ENABLE_QSCORE_OUTPUT);
            final LocusReadPile normalReadPile = new LocusReadPile(normalPileup, upRef, MTAC.MIN_QSCORE, 0, this.USE_MAPQ0_IN_NORMAL_QSCORE, true, MTAC.ENABLE_QSCORE_OUTPUT);

            Collection<VariantContext> panelOfNormalsVC = tracker.getValues(normalPanelRod, rawContext.getLocation());
            Collection<VariantContext> cosmicVC = getVCInTrackerInLocus(RODNames.COSMIC, tracker);
            Collection<VariantContext> dbsnpVC = getVCInTrackerInLocus(RODNames.DBSNP, tracker);

            // remove the effect of cosmic from dbSNP
            boolean germlineAtRisk = (!dbsnpVC.isEmpty() && cosmicVC.isEmpty());

            // compute coverage flags
            int tumorCoveredDepthThreshold = 14;
            int normalCoveredDepthThreshold = (germlineAtRisk) ? 19 : 8;
            if (!hasNormalBam) {
                normalCoveredDepthThreshold = 0;
            }

            int tumorBaseCount = tumorReadPile.finalPileupReads.size();
            int normalBaseCount = normalReadPile.finalPileupReads.size();
            boolean isTumorCovered = tumorBaseCount >= tumorCoveredDepthThreshold;
            boolean isNormalCovered = normalBaseCount >= normalCoveredDepthThreshold;
            boolean isBaseCovered = isTumorCovered && isNormalCovered;
            if (!hasNormalBam) {
                isBaseCovered = isTumorCovered;
            }

            int tumorQ20BaseCount = tumorReadPile.getFilteredBaseCount(20);
            int normalQ20BaseCount = normalReadPile.getFilteredBaseCount(20);

            // calculate power
            double tumorPower = tumorPowerCalculator.cachingPowerCalculation(tumorBaseCount, MTAC.POWER_CONSTANT_AF);

            double normalPowerNoSNPPrior = normalNovelSitePowerCalculator.cachingPowerCalculation(normalBaseCount);
            double normalPowerWithSNPPrior = normalDbSNPSitePowerCalculator.cachingPowerCalculation(normalBaseCount);

            double normalPower = (germlineAtRisk) ? normalPowerWithSNPPrior : normalPowerNoSNPPrior;

            double combinedPower = tumorPower * normalPower;
            if (!hasNormalBam) {
                combinedPower = tumorPower;
            }

            int mapQ0Reads =
                    tumorReadPile.qualityScoreFilteredPileup.getNumberOfMappingQualityZeroReads() +
                            normalReadPile.qualityScoreFilteredPileup.getNumberOfMappingQualityZeroReads();

            int totalReads =
                    tumorReadPile.qualityScoreFilteredPileup.depthOfCoverage() +
                            normalReadPile.qualityScoreFilteredPileup.depthOfCoverage();

            // Test each of the possible alternate alleles
            for (final char altAllele : new char[]{'A', 'C', 'G', 'T'}) {
                if (altAllele == upRef) {
                    continue;
                }
                if (!MTAC.FORCE_OUTPUT && tumorReadPile.qualitySums.getCounts(altAllele) == 0) {
                    continue;
                }

                CandidateMutation candidate = new CandidateMutation(rawContext.getLocation(), upRef);
                candidate.setSequenceContext(sequenceContext);
                candidate.setTumorSampleName(MTAC.TUMOR_SAMPLE_NAME);
                candidate.setNormalSampleName(MTAC.NORMAL_SAMPLE_NAME);
                candidate.setCovered(isBaseCovered);
                candidate.setPower(combinedPower);
                candidate.setTumorPower(tumorPower);
                candidate.setNormalPower(normalPower);
                candidate.setNormalPowerWithSNPPrior(normalPowerWithSNPPrior);
                candidate.setNormalPowerNoSNPPrior(normalPowerNoSNPPrior);
                candidate.setTumorQ20Count(tumorQ20BaseCount);
                candidate.setNormalQ20Count(normalQ20BaseCount);
                candidate.setInitialTumorNonRefQualitySum(tumorReadPile.qualitySums.getOtherQualities(upRef));
                candidate.setAltAllele(altAllele);
                candidate.setMapQ0Reads(mapQ0Reads);
                candidate.setTotalReads(totalReads);
                candidate.setContaminationFraction(MTAC.FRACTION_CONTAMINATION);
                candidate.setPanelOfNormalsVC(panelOfNormalsVC.isEmpty() ? null : panelOfNormalsVC.iterator().next()); // if there are multiple, we're just grabbing the first
                candidate.setCosmicSite(!cosmicVC.isEmpty());
                candidate.setDbsnpSite(!dbsnpVC.isEmpty());
                candidate.setDbsnpVC(dbsnpVC.isEmpty() ? null : dbsnpVC.iterator().next());
                candidate.setTumorF(tumorReadPile.estimateAlleleFraction(upRef, altAllele));

                if (!MTAC.FORCE_OUTPUT && candidate.getTumorF() < MTAC.TUMOR_F_PRETEST) {
                    continue;
                }

                candidate.setInitialTumorAltCounts(tumorReadPile.qualitySums.getCounts(altAllele));
                candidate.setInitialTumorRefCounts(tumorReadPile.qualitySums.getCounts(upRef));
                candidate.setInitialTumorAltQualitySum(tumorReadPile.qualitySums.getQualitySum(altAllele));
                candidate.setInitialTumorRefQualitySum(tumorReadPile.qualitySums.getQualitySum(upRef));

                double tumorLod = tumorReadPile.calculateAltVsRefLOD((byte) altAllele, candidate.getTumorF(), 0);
                candidate.setTumorLodFStar(tumorLod);

                candidate.setInitialTumorReadDepth(tumorReadPile.finalPileupReads.size());
                candidate.setTumorInsertionCount(tumorReadPile.getInsertionsCount());
                candidate.setTumorDeletionCount(tumorReadPile.getDeletionsCount());

                if (candidate.getTumorLodFStar() < MTAC.INITIAL_TUMOR_LOD_THRESHOLD) {
                    continue;
                }

                // calculate lod of contaminant
                double contaminantF = Math.min(contaminantAlternateFraction, candidate.getTumorF());
                VariableAllelicRatioGenotypeLikelihoods contaminantLikelihoods
                        = new VariableAllelicRatioGenotypeLikelihoods(upRef, contaminantF);

                List<PileupElement> peList = new ArrayList<PileupElement>(tumorReadPile.finalPileup.depthOfCoverage());
                for (PileupElement pe : tumorReadPile.finalPileup) {
                    peList.add(pe);
                }

                Collections.sort(peList, new PileupComparatorByAltRefQual((byte) altAllele));
                int readsToKeep = (int) (peList.size() * contaminantAlternateFraction);

                for (PileupElement pe : peList) {
                    byte base = pe.getBase();
                    if (pe.getBase() == altAllele) {
                        // if we've retained all we need, then turn the remainder of alts to ref
                        if (readsToKeep == 0) {
                            base = (byte) upRef;
                        } else {
                            readsToKeep--;
                        }
                    }

                    contaminantLikelihoods.add(base, pe.getQual());
                }
                double[] refHetHom = LocusReadPile.extractRefHetHom(contaminantLikelihoods, upRef, altAllele);
                double contaminantLod = refHetHom[1] - refHetHom[0];
                candidate.setContaminantLod(contaminantLod);

                final QualitySums normQs = normalReadPile.qualitySums;


                VariableAllelicRatioGenotypeLikelihoods normalGl = normalReadPile.calculateLikelihoods(normalReadPile.qualityScoreFilteredPileup); // use MAPQ0 reads
                candidate.setInitialNormalBestGenotype(normalReadPile.getBestGenotype(normalGl));
                candidate.setInitialNormalLod(LocusReadPile.getRefVsAlt(normalGl, upRef, altAllele));


                double normalF = Math.max(LocusReadPile.estimateAlleleFraction(normalReadPile.qualityScoreFilteredPileup, upRef, altAllele), MTAC.MINIMUM_NORMAL_ALLELE_FRACTION);
                candidate.setNormalF(normalF);


                candidate.setInitialNormalAltQualitySum(normQs.getQualitySum(altAllele));
                candidate.setInitialNormalRefQualitySum(normQs.getQualitySum(upRef));

                candidate.setNormalAltQualityScores(normQs.getBaseQualityScores(altAllele));
                candidate.setNormalRefQualityScores(normQs.getBaseQualityScores(upRef));

                candidate.setInitialNormalAltCounts(normQs.getCounts(altAllele));
                candidate.setInitialNormalRefCounts(normQs.getCounts(upRef));
                candidate.setInitialNormalReadDepth(normalReadPile.finalPileupReads.size());

                // TODO: parameterize filtering Mate-Rescued Reads (if someone wants to disable this)
                final LocusReadPile t2 = filterReads(ref, tumorReadPile.finalPileup, true);

                // if there are no reads remaining, abandon this theory
                if (!MTAC.FORCE_OUTPUT && t2.finalPileupReads.size() == 0) {
                    continue;
                }

                candidate.setInitialTumorAltCounts(t2.qualitySums.getCounts(altAllele));
                candidate.setInitialTumorRefCounts(t2.qualitySums.getCounts(upRef));
                candidate.setInitialTumorAltQualitySum(t2.qualitySums.getQualitySum(altAllele));
                candidate.setInitialTumorRefQualitySum(t2.qualitySums.getQualitySum(upRef));

                candidate.setTumorAltQualityScores(t2.qualitySums.getBaseQualityScores(altAllele));
                candidate.setTumorRefQualityScores(t2.qualitySums.getBaseQualityScores(upRef));

                VariableAllelicRatioGenotypeLikelihoods t2Gl = t2.calculateLikelihoods(t2.finalPileup);
                candidate.setInitialTumorLod(t2.getAltVsRef(t2Gl, upRef, altAllele));
                candidate.setInitialTumorReadDepth(t2.finalPileupReads.size());

                candidate.setTumorF(t2.estimateAlleleFraction(upRef, altAllele));
                double tumorLod2 = t2.calculateAltVsRefLOD((byte) altAllele, candidate.getTumorF(), 0);
                candidate.setTumorLodFStar(tumorLod2);

                //TODO: clean up use of forward/reverse vs positive/negative (prefer the latter since GATK uses it)
                ReadBackedPileup forwardPileup = filterReads(ref, tumorReadPile.finalPileupPositiveStrand, true).finalPileupPositiveStrand;
                double f2forward = LocusReadPile.estimateAlleleFraction(forwardPileup, upRef, altAllele);
                candidate.setTumorLodFStarForward(t2.calculateAltVsRefLOD(forwardPileup, (byte) altAllele, f2forward, 0.0));

                ReadBackedPileup reversePileup = filterReads(ref, tumorReadPile.finalPileupNegativeStrand, true).finalPileupNegativeStrand;
                double f2reverse = LocusReadPile.estimateAlleleFraction(reversePileup, upRef, altAllele);
                candidate.setTumorLodFStarReverse(t2.calculateAltVsRefLOD(reversePileup, (byte) altAllele, f2reverse, 0.0));

                // calculate strand bias power
                candidate.setPowerToDetectPositiveStrandArtifact(
                        strandArtifactPowerCalculator.cachingPowerCalculation(reversePileup.depthOfCoverage(), candidate.getTumorF())
                );
                candidate.setPowerToDetectNegativeStrandArtifact(
                        strandArtifactPowerCalculator.cachingPowerCalculation(forwardPileup.depthOfCoverage(), candidate.getTumorF())
                );

                candidate.setStrandContingencyTable(SequenceUtils.getStrandContingencyTable(forwardPileup, reversePileup, (byte) upRef, (byte) altAllele));

                ArrayList<PileupElement> mutantPileupElements = new ArrayList<PileupElement>();
                ArrayList<PileupElement> referencePileupElements = new ArrayList<PileupElement>();


                for (PileupElement p : t2.finalPileup) {
                    final SAMRecord read = p.getRead();
                    final int offset = p.getOffset();

                    if (read.getReadString().charAt(offset) == altAllele) {
                        mutantPileupElements.add(p);
                    } else if (read.getReadString().charAt(offset) == upRef) {
                        referencePileupElements.add(p);
                    } else {
                        // just drop the read...
                    }
                }

                ReadBackedPileup mutantPileup =
                        new ReadBackedPileupImpl(rawContext.getLocation(), mutantPileupElements);

                ReadBackedPileup referencePileup =
                        new ReadBackedPileupImpl(rawContext.getLocation(), referencePileupElements);

                // TODO: shouldn't this be refAllele here?
                final LocusReadPile mutantPile = new LocusReadPile(mutantPileup, altAllele, 0, 0, MTAC.ENABLE_QSCORE_OUTPUT);
                final LocusReadPile refPile = new LocusReadPile(referencePileup, altAllele, 0, 0, MTAC.ENABLE_QSCORE_OUTPUT);

                // Set the maximum observed mapping quality score for the reference and alternate alleles
                int[] rmq = referencePileup.getMappingQuals();
                candidate.setTumorRefMaxMapQ((rmq.length == 0) ? 0 : NumberUtils.max(rmq));

                int[] amq = mutantPileup.getMappingQuals();
                candidate.setTumorAltMaxMapQ((amq.length == 0) ? 0 : NumberUtils.max(amq));


                // start with just the tumor pile
                candidate.setTumorAltForwardOffsetsInRead(SequenceUtils.getForwardOffsetsInRead(mutantPileup));
                candidate.setTumorAltReverseOffsetsInRead(SequenceUtils.getReverseOffsetsInRead(mutantPileup));

                if (candidate.getTumorAltForwardOffsetsInRead().size() > 0) {
                    double[] offsets = MuTectStats.convertIntegersToDoubles(candidate.getTumorAltForwardOffsetsInRead());
                    double median = MuTectStats.getMedian(offsets);
                    candidate.setTumorForwardOffsetsInReadMedian(median);
                    candidate.setTumorForwardOffsetsInReadMad(MuTectStats.calculateMAD(offsets, median));
                }


                if (candidate.getTumorAltReverseOffsetsInRead().size() > 0) {
                    double[] offsets = MuTectStats.convertIntegersToDoubles(candidate.getTumorAltReverseOffsetsInRead());
                    double median = MuTectStats.getMedian(offsets);
                    candidate.setTumorReverseOffsetsInReadMedian(median);
                    candidate.setTumorReverseOffsetsInReadMad(MuTectStats.calculateMAD(offsets, median));
                }

                // test to see if the candidate should be rejected
                performRejection(candidate);

                messageByTumorLod.put(candidate.getInitialTumorLod(), candidate);

            }

            // if more than one site passes the tumor lod threshold for KEEP the fail the tri_allelic Site filter
            int passingCandidates = 0;
            for (CandidateMutation c : messageByTumorLod.values()) {
                if (c.getTumorLodFStar() >= MTAC.TUMOR_LOD_THRESHOLD) {
                    passingCandidates++;
                }
            }

            if (passingCandidates > 1) {
                for (CandidateMutation c : messageByTumorLod.values()) {
                    c.addRejectionReason("triallelic_site");
                }
            }

            // write out the call stats for the "best" candidate
            if (!messageByTumorLod.isEmpty()) {
                CandidateMutation m = messageByTumorLod.lastEntry().getValue();

                // only output passing calls OR rejected sites if ONLY_PASSING_CALLS is not specified
                if (!m.isRejected() || (m.isRejected() && !MTAC.ONLY_PASSING_CALLS)) {

                    //out.println(callStatsGenerator.generateCallStats(m));
                    resultVCFOutInfos.add(callStatsGenerator.generateCallStats(m));
                    resultVCFRecords.add(VCFGenerator.generateVC(m));
                }
            }
        } catch (MathException me) {
            throw new GATKException(me.getMessage());
        }
    }

    int MAX_READ_MISMATCH_QUALITY_SCORE_SUM = 100;
    private static Character MAPPED_BY_MATE = 'M';

    private LocusReadPile filterReads(final ReferenceContext ref, final ReadBackedPileup pile, boolean filterMateRescueReads) {
        ArrayList<PileupElement> newPileupElements = new ArrayList<PileupElement>();

        for (PileupElement p : pile) {
            final GATKSAMRecord read = p.getRead();

            int mismatchQualitySum =
                    CGAAlignmentUtils.mismatchesInRefWindow(p, ref, refContentProvider, false, true);

            // do we have to many mismatches overall?
            if (mismatchQualitySum > this.MAX_READ_MISMATCH_QUALITY_SCORE_SUM) {
                continue;
            }

            // is this a heavily clipped read?
            if (SequenceUtils.isReadHeavilySoftClipped(read, MTAC.HEAVILY_CLIPPED_READ_FRACTION)) {
                continue;
            }

            // was this read ONLY placed because it's mate was uniquely placed? (supplied by BWA)
            if (filterMateRescueReads && MAPPED_BY_MATE.equals(read.getAttribute("XT"))) {
                continue;
            }

            // if we're here... we passed all the read filters!
            newPileupElements.add(new PileupElement(p));

        }
        ReadBackedPileup newPileup =
                new ReadBackedPileupImpl(ref.getLocus(), newPileupElements);


        return new LocusReadPile(newPileup, (char) ref.getBase(), 0, 0, MTAC.ENABLE_QSCORE_OUTPUT);
    }

    private void performRejection(CandidateMutation candidate) {
        if (candidate.getTumorLodFStar() < MTAC.TUMOR_LOD_THRESHOLD) {
            candidate.addRejectionReason("fstar_tumor_lod");
        }

        if (MTAC.ARTIFACT_DETECTION_MODE) {
            return;
        }


        if (candidate.getTumorInsertionCount() >= MTAC.GAP_EVENTS_THRESHOLD ||
                candidate.getTumorDeletionCount() >= MTAC.GAP_EVENTS_THRESHOLD) {
            candidate.addRejectionReason("nearby_gap_events");
        }

        if (MTAC.FRACTION_CONTAMINATION + MTAC.MINIMUM_MUTATION_CELL_FRACTION > 0 && candidate.getTumorLodFStar() <= MTAC.TUMOR_LOD_THRESHOLD + Math.max(0, candidate.getContaminantLod())) {
            candidate.addRejectionReason("possible_contamination");
        }

        if (candidate.isGermlineAtRisk() && candidate.getInitialNormalLod() < MTAC.NORMAL_DBSNP_LOD_THRESHOLD) {
            candidate.addRejectionReason("germline_risk");
        }

        if (candidate.getInitialNormalLod() < MTAC.NORMAL_LOD_THRESHOLD) {
            candidate.addRejectionReason("normal_lod");
        }

        if ((candidate.getInitialNormalAltCounts() >= MTAC.MAX_ALT_ALLELES_IN_NORMAL_COUNT || candidate.getNormalF() >= MTAC.MAX_ALT_ALLELE_IN_NORMAL_FRACTION) && candidate.getInitialNormalAltQualitySum() > MTAC.MAX_ALT_ALLELES_IN_NORMAL_QSCORE_SUM) {
            candidate.addRejectionReason("alt_allele_in_normal");
        }

        if ((candidate.getTumorForwardOffsetsInReadMedian() != null && candidate.getTumorForwardOffsetsInReadMedian() <= MTAC.PIR_MEDIAN_THRESHOLD && candidate.getTumorForwardOffsetsInReadMad() != null && candidate.getTumorForwardOffsetsInReadMad() <= MTAC.PIR_MAD_THRESHOLD) ||
                candidate.getTumorReverseOffsetsInReadMedian() != null && candidate.getTumorReverseOffsetsInReadMedian() <= MTAC.PIR_MEDIAN_THRESHOLD && candidate.getTumorReverseOffsetsInReadMad() != null && candidate.getTumorReverseOffsetsInReadMad() <= MTAC.PIR_MAD_THRESHOLD) {
            candidate.addRejectionReason("clustered_read_position");

        }

        // TODO: sync naming (is it positive or forward)?
        if (
                (candidate.getPowerToDetectNegativeStrandArtifact() >= MTAC.STRAND_ARTIFACT_POWER_THRESHOLD && candidate.getTumorLodFStarForward() < MTAC.STRAND_ARTIFACT_LOD_THRESHOLD) ||
                        (candidate.getPowerToDetectPositiveStrandArtifact() >= MTAC.STRAND_ARTIFACT_POWER_THRESHOLD && candidate.getTumorLodFStarReverse() < MTAC.STRAND_ARTIFACT_LOD_THRESHOLD)
                ) {
            candidate.addRejectionReason("strand_artifact");
        }

        if (candidate.getTotalReads() > 0 && ((float) candidate.getMapQ0Reads() / (float) candidate.getTotalReads()) >= MTAC.FRACTION_MAPQ0_THRESHOLD) {
            candidate.addRejectionReason("poor_mapping_region_mapq0");
        }

        if (candidate.getTumorAltMaxMapQ() < MTAC.REQUIRED_MAXIMUM_ALT_ALLELE_MAPPING_QUALITY_SCORE) {
            candidate.addRejectionReason("poor_mapping_region_alternate_allele_mapq");
        }

        if (candidate.isSeenInPanelOfNormals()) {
            if (candidate.isCosmicSite()) {
                // if we saw it in the panel of normals, retain the call it was a COSMIC, but non-dbsnp site,
            } else {
                // otherwise, reject it
                candidate.addRejectionReason("seen_in_panel_of_normals");
            }
        }

    }


    @Override
    protected void onTraversalDone() {

    }

    public List<VariantContext> getResultVCFRecords() {
        return resultVCFRecords;
    }

    public List<String> getResultVCFOutInfos() {
        return resultVCFOutInfos;
    }

    public static class PileupComparatorByAltRefQual implements Comparator<PileupElement> {
        private byte alt;

        public PileupComparatorByAltRefQual(byte alt) {
            this.alt = alt;
        }

        public int compare(PileupElement o1, PileupElement o2) {
            return internalCompare(o1.getBase(), o1.getQual(), o2.getBase(), o2.getQual());
        }

        public int internalCompare(byte base1, byte qual1, byte base2, byte qual2) {
            // if the bases are the same, the higher quality score comes first
            if (base1 == base2) {
                if (qual1 == qual2) {
                    return 0;
                }
                return (qual1 > qual2) ? -1 : 1;

                // if the bases are not the same, then the alternate is first
            } else {
                if (base1 == alt) {
                    return -1;
                } else if (base2 == alt) {
                    return 1;
                } else {
                    return base1 < base2 ? -1 : 1;
                }

            }
        }

    }
}
