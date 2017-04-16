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
package org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.readlikelihood;

import htsjdk.samtools.SAMUtils;
import htsjdk.variant.variantcontext.Allele;
import org.apache.log4j.Logger;
import org.ncic.bioinfo.sparkseq.algorithms.utils.MathUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.QualityUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.haplotype.Haplotype;
import org.ncic.bioinfo.sparkseq.algorithms.utils.pairhmm.ArrayLoglessPairHMM;
import org.ncic.bioinfo.sparkseq.algorithms.utils.pairhmm.Log10PairHMM;
import org.ncic.bioinfo.sparkseq.algorithms.utils.pairhmm.LoglessPairHMM;
import org.ncic.bioinfo.sparkseq.algorithms.utils.pairhmm.PairHMM;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.GATKSAMRecord;
import org.ncic.bioinfo.sparkseq.algorithms.walker.baserecalibrator.covariate.RepeatCovariate;
import org.ncic.bioinfo.sparkseq.algorithms.walker.baserecalibrator.covariate.RepeatLengthCovariate;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.AlleleList;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.AssemblyResultSet;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.IndexedAlleleList;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.SampleList;
import org.ncic.bioinfo.sparkseq.exceptions.UserException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Set;

/**
 * Author: wbc
 */
public class PairHMMLikelihoodCalculationEngine implements ReadLikelihoodCalculationEngine {
    private final static Logger logger = Logger.getLogger(PairHMMLikelihoodCalculationEngine.class);

    public static final byte BASE_QUALITY_SCORE_THRESHOLD = (byte) 18; // Base quals less than this value are squashed down to min possible qual

    private final byte constantGCP;

    private final double log10globalReadMismappingRate;

    private final PairHMM.HMM_IMPLEMENTATION hmmType;
    private final boolean noFpga;

    private final ThreadLocal<PairHMM> pairHMMThreadLocal = new ThreadLocal<PairHMM>() {
        @Override
        protected PairHMM initialValue() {
            switch (hmmType) {
                case EXACT:
                    return new Log10PairHMM(true);
                case ORIGINAL:
                    return new Log10PairHMM(false);
                case LOGLESS_CACHING:
                    return new LoglessPairHMM();
                // TODO by wbc 暂时不支持这俩，因为要外部代码，这里只先做一个纯的java工程
                case VECTOR_LOGLESS_CACHING:
                    return new LoglessPairHMM();
                case DEBUG_VECTOR_LOGLESS_CACHING:
                    return new LoglessPairHMM();
                case ARRAY_LOGLESS:
                    return new ArrayLoglessPairHMM();
                default:
                    throw new UserException.BadArgumentValue("pairHMM", "Specified pairHMM implementation is unrecognized or incompatible with the HaplotypeCaller. Acceptable options are ORIGINAL, EXACT, CACHING, LOGLESS_CACHING, and ARRAY_LOGLESS.");
            }
        }
    };
//    Attempted to do as below, to avoid calling pairHMMThreadLocal.get() later on, but it resulted in a NullPointerException
//    private final PairHMM pairHMM = pairHMMThreadLocal.get();

    private final static boolean WRITE_LIKELIHOODS_TO_FILE = false;
    private final static String LIKELIHOODS_FILENAME = "likelihoods.txt";
    private final PrintStream likelihoodsStream;

    public enum PCR_ERROR_MODEL {
        /**
         * no specialized PCR error model will be applied; if base insertion/deletion qualities are present they will be used
         */
        NONE,
        /**
         * a more aggressive model will be applied that sacrifices true positives in order to remove more false positives
         */
        AGGRESSIVE,
        /**
         * a less aggressive model will be applied that tries to maintain a high true positive rate at the expense of allowing more false positives
         */
        CONSERVATIVE
    }

    private final PCR_ERROR_MODEL pcrErrorModel;

    /**
     * The expected rate of random sequencing errors for a read originating from its true haplotype.
     * <p>
     * For example, if this is 0.01, then we'd expect 1 error per 100 bp.
     */
    private final static double EXPECTED_ERROR_RATE_PER_BASE = 0.02;

    /**
     * Create a new PairHMMLikelihoodCalculationEngine using provided parameters and hmm to do its calculations
     *
     * @param constantGCP                   the gap continuation penalty to use with the PairHMM
     * @param hmmType                       the type of the HMM to use
     * @param log10globalReadMismappingRate the global mismapping probability, in log10(prob) units.  A value of
     *                                      -3 means that the chance that a read doesn't actually belong at this
     *                                      location in the genome is 1 in 1000.  The effect of this parameter is
     *                                      to cap the maximum likelihood difference between the reference haplotype
     *                                      and the best alternative haplotype by -3 log units.  So if the best
     *                                      haplotype is at -10 and this parameter has a value of -3 then even if the
     *                                      reference haplotype gets a score of -100 from the pairhmm it will be
     *                                      assigned a likelihood of -13.
     * @param noFpga                        disable FPGA acceleration
     */
    public PairHMMLikelihoodCalculationEngine(final byte constantGCP, final PairHMM.HMM_IMPLEMENTATION hmmType, final double log10globalReadMismappingRate, final boolean noFpga, final PCR_ERROR_MODEL pcrErrorModel) {
        this.hmmType = hmmType;
        this.constantGCP = constantGCP;
        this.log10globalReadMismappingRate = log10globalReadMismappingRate;
        this.noFpga = noFpga;
        this.pcrErrorModel = pcrErrorModel;

        initializePCRErrorModel();

        if (WRITE_LIKELIHOODS_TO_FILE) {
            try {
                likelihoodsStream = new PrintStream(new FileOutputStream(new File(LIKELIHOODS_FILENAME)));
            } catch (FileNotFoundException e) {
                throw new RuntimeException(e);
            }
        } else {
            likelihoodsStream = null;
        }
    }

    @Override
    public void close() {
        if (likelihoodsStream != null) likelihoodsStream.close();
        pairHMMThreadLocal.get().close();
    }

    private void capMinimumReadQualities(GATKSAMRecord read, byte[] readQuals, byte[] readInsQuals, byte[] readDelQuals) {
        for (int kkk = 0; kkk < readQuals.length; kkk++) {
            readQuals[kkk] = (byte) Math.min(0xff & readQuals[kkk], read.getMappingQuality()); // cap base quality by mapping quality, as in UG
            readQuals[kkk] = (readQuals[kkk] < BASE_QUALITY_SCORE_THRESHOLD ? QualityUtils.MIN_USABLE_Q_SCORE : readQuals[kkk]);
            readInsQuals[kkk] = (readInsQuals[kkk] < QualityUtils.MIN_USABLE_Q_SCORE ? QualityUtils.MIN_USABLE_Q_SCORE : readInsQuals[kkk]);
            readDelQuals[kkk] = (readDelQuals[kkk] < QualityUtils.MIN_USABLE_Q_SCORE ? QualityUtils.MIN_USABLE_Q_SCORE : readDelQuals[kkk]);
        }
    }

    /**
     * Pre-processing of the reads to be evaluated at the current location from the current sample.
     * We apply the PCR Error Model, and cap the minimum base, insertion, and deletion qualities of each read.
     * Modified copies of reads are packed into a new list, while original reads are retained for downstream use
     *
     * @param reads The original list of unmodified reads
     * @return processedReads. A new list of reads, in the same order, whose qualities have been altered by PCR error model and minimal quality thresholding
     */
    private List<GATKSAMRecord> modifyReadQualities(final List<GATKSAMRecord> reads) {
        final List<GATKSAMRecord> result = new ArrayList<>(reads.size());

        for (final GATKSAMRecord read : reads) {
            final byte[] readBases = read.getReadBases();

            // NOTE -- must clone anything that gets modified here so we don't screw up future uses of the read
            final byte[] readQuals = read.getBaseQualities().clone();
            final byte[] readInsQuals = read.getBaseInsertionQualities().clone();
            final byte[] readDelQuals = read.getBaseDeletionQualities().clone();

            applyPCRErrorModel(readBases, readInsQuals, readDelQuals);
            capMinimumReadQualities(read, readQuals, readInsQuals, readDelQuals);

            // Create a new copy of the read and sets its base qualities to the modified versions.
            // Pack this into a new list for return
            result.add(GATKSAMRecord.createQualityModifiedRead(read, readBases, readQuals, readInsQuals, readDelQuals));
        }
        return result;
    }

    /**
     * Initialize our pairHMM with parameters appropriate to the haplotypes and reads we're going to evaluate
     * <p>
     * After calling this routine the PairHMM will be configured to best evaluate all reads in the samples
     * against the set of haplotypes
     *
     * @param haplotypes        a non-null list of haplotypes
     * @param perSampleReadList a mapping from sample -> reads
     */
    private void initializePairHMM(final List<Haplotype> haplotypes, final Map<String, List<GATKSAMRecord>> perSampleReadList) {
        int X_METRIC_LENGTH = 0;
        for (final Map.Entry<String, List<GATKSAMRecord>> sample : perSampleReadList.entrySet()) {
            for (final GATKSAMRecord read : sample.getValue()) {
                final int readLength = read.getReadLength();
                if (readLength > X_METRIC_LENGTH) {
                    X_METRIC_LENGTH = readLength;
                }
            }
        }
        int Y_METRIC_LENGTH = 0;
        for (final Haplotype h : haplotypes) {
            final int haplotypeLength = h.getBases().length;
            if (haplotypeLength > Y_METRIC_LENGTH) {
                Y_METRIC_LENGTH = haplotypeLength;
            }
        }

        // initialize arrays to hold the probabilities of being in the match, insertion and deletion cases
        pairHMMThreadLocal.get().initialize(haplotypes, perSampleReadList, X_METRIC_LENGTH, Y_METRIC_LENGTH);
    }

    private void finalizePairHMM() {
        pairHMMThreadLocal.get().finalizeRegion();
    }

    @Override
    public ReadLikelihoods<Haplotype> computeReadLikelihoods(final AssemblyResultSet assemblyResultSet, final SampleList samples, final Map<String, List<GATKSAMRecord>> perSampleReadList) {

        final List<Haplotype> haplotypeList = assemblyResultSet.getHaplotypeList();
        final AlleleList<Haplotype> haplotypes = new IndexedAlleleList<>(haplotypeList);

        // configure the HMM
        initializePairHMM(haplotypeList, perSampleReadList);

        // Add likelihoods for each sample's reads to our result
        final ReadLikelihoods<Haplotype> result = new ReadLikelihoods<>(samples, haplotypes, perSampleReadList, false);
        final int sampleCount = result.sampleCount();
        for (int s = 0; s < sampleCount; s++) {
            final ReadLikelihoods.Matrix<Haplotype> sampleLikelihoods = result.sampleMatrix(s);
            computeReadLikelihoods(sampleLikelihoods);
        }

        result.normalizeLikelihoods(false, log10globalReadMismappingRate);
        result.filterPoorlyModeledReads(EXPECTED_ERROR_RATE_PER_BASE);
        finalizePairHMM();
        return result;
    }

    private void computeReadLikelihoods(final ReadLikelihoods.Matrix<Haplotype> likelihoods) {

        // Modify the read qualities by applying the PCR error model and capping the minimum base,insertion,deletion qualities
        final List<GATKSAMRecord> processedReads = modifyReadQualities(likelihoods.reads());

        final Map<GATKSAMRecord, byte[]> gapContinuationPenalties = buildGapContinuationPenalties(processedReads, constantGCP);
        // Run the PairHMM to calculate the log10 likelihood of each (processed) reads' arising from each haplotype
        pairHMMThreadLocal.get().computeLikelihoods(likelihoods, processedReads, gapContinuationPenalties);

        if (WRITE_LIKELIHOODS_TO_FILE)
            writeDebugLikelihoods(likelihoods);
    }

    private Map<GATKSAMRecord, byte[]> buildGapContinuationPenalties(final List<GATKSAMRecord> processedReads, final byte gcp) {
        final Map<GATKSAMRecord, byte[]> result = new HashMap<>(processedReads.size());
        for (final GATKSAMRecord read : processedReads) {
            final byte[] readGcpArray = new byte[read.getReadLength()];
            Arrays.fill(readGcpArray, gcp);
            result.put(read, readGcpArray);
        }
        return result;
    }

    private void writeDebugLikelihoods(final ReadLikelihoods.Matrix<Haplotype> likelihoods) {
        final List<GATKSAMRecord> reads = likelihoods.reads();
        final List<Haplotype> haplotypes = likelihoods.alleles();
        final int haplotypeCount = haplotypes.size();
        final int readCount = reads.size();
        for (int r = 0; r < readCount; r++)
            for (int a = 0; a < haplotypeCount; a++)
                writeDebugLikelihoods(reads.get(r), haplotypes.get(a), likelihoods.get(a, r));
        likelihoodsStream.flush();
    }

    private void writeDebugLikelihoods(final GATKSAMRecord processedRead, final Haplotype haplotype, final double log10l) {
        likelihoodsStream.printf("%s %s %s %s %s %s %f%n",
                haplotype.getBaseString(),
                new String(processedRead.getReadBases()),
                SAMUtils.phredToFastq(processedRead.getBaseQualities()),
                SAMUtils.phredToFastq(processedRead.getBaseInsertionQualities()),
                SAMUtils.phredToFastq(processedRead.getBaseDeletionQualities()),
                SAMUtils.phredToFastq(constantGCP),
                log10l);
    }

    @Deprecated
    public static double[][] computeDiploidHaplotypeLikelihoods(final String sample,
                                                                final ReadLikelihoods readLikelihoods,
                                                                final List alleleOrdering,
                                                                final boolean normalize) {
        return computeDiploidHaplotypeLikelihoods(Collections.singleton(sample), readLikelihoods, alleleOrdering, normalize);
    }

    @Deprecated
    private static double[][] computeDiploidHaplotypeLikelihoods(final Set<String> samples,
                                                                 final ReadLikelihoods readLikelihoods,
                                                                 final List alleleOrdering,
                                                                 final boolean normalize) {

        final int numHaplotypes = alleleOrdering.size();
        final int[] alleleIndices = new int[alleleOrdering.size()];
        final ListIterator alleleIterator = alleleOrdering.listIterator();
        int nextAlleleIndex = 0;
        while (alleleIterator.hasNext())
            if ((alleleIndices[nextAlleleIndex++] = readLikelihoods.alleleIndex((Allele) alleleIterator.next())) == -1)
                throw new IllegalArgumentException("allele " + alleleIterator.previous() + " not found in likelihood collection ");

        final double[][] haplotypeLikelihoodMatrix = new double[numHaplotypes][numHaplotypes];

        // compute the diploid haplotype likelihoods
        for (final String sample : samples) {
            final int sampleIndex = readLikelihoods.sampleIndex(sample);
            if (sampleIndex == -1)
                throw new IllegalArgumentException("the sample provided is not in the likelihood collection");
            final ReadLikelihoods.Matrix sampleLikelihoods = readLikelihoods.sampleMatrix(sampleIndex);
            final int sampleReadCount = readLikelihoods.sampleReadCount(sampleIndex);
            for (int iii = 0; iii < numHaplotypes; iii++) {
                final int iii_allele = alleleIndices[iii];
                for (int jjj = 0; jjj <= iii; jjj++) {
                    final int jjj_allele = alleleIndices[jjj];
                    double haplotypeLikelihood = 0.0;
                    for (int r = 0; r < sampleReadCount; r++) {
                        final double value = MathUtils.approximateLog10SumLog10(sampleLikelihoods.get(iii_allele, r),
                                sampleLikelihoods.get(jjj_allele, r)) + MathUtils.LOG_ONE_HALF;
                        haplotypeLikelihood += value;
                    }
                    haplotypeLikelihoodMatrix[iii][jjj] += haplotypeLikelihood;
                }
            }
        }

        // normalize the diploid likelihoods matrix
        return normalize ? normalizeDiploidLikelihoodMatrixFromLog10(haplotypeLikelihoodMatrix) : haplotypeLikelihoodMatrix;
    }

    @Deprecated
    protected static double[][] normalizeDiploidLikelihoodMatrixFromLog10(final double[][] likelihoodMatrix) {
        final int numHaplotypes = likelihoodMatrix.length;
        double[] genotypeLikelihoods = new double[numHaplotypes * (numHaplotypes + 1) / 2];
        int index = 0;
        for (int iii = 0; iii < numHaplotypes; iii++) {
            for (int jjj = 0; jjj <= iii; jjj++) {
                genotypeLikelihoods[index++] = likelihoodMatrix[iii][jjj];
            }
        }
        genotypeLikelihoods = MathUtils.normalizeFromLog10(genotypeLikelihoods, false, true);
        index = 0;
        for (int iii = 0; iii < numHaplotypes; iii++) {
            for (int jjj = 0; jjj <= iii; jjj++) {
                likelihoodMatrix[iii][jjj] = genotypeLikelihoods[index++];
            }
        }
        return likelihoodMatrix;
    }

    // --------------------------------------------------------------------------------
    //
    // Experimental attempts at PCR error rate modeling
    //
    // --------------------------------------------------------------------------------

    protected static final int MAX_STR_UNIT_LENGTH = 8;
    protected static final int MAX_REPEAT_LENGTH = 20;
    protected static final int MIN_ADJUSTED_QSCORE = 10;
    protected static final double INITIAL_QSCORE = 40.0;

    private byte[] pcrIndelErrorModelCache = new byte[MAX_REPEAT_LENGTH * MAX_STR_UNIT_LENGTH + 1];
    private final RepeatCovariate repeatCovariate = new RepeatLengthCovariate();

    private void initializePCRErrorModel() {
        if (pcrErrorModel == PCR_ERROR_MODEL.NONE)
            return;

        repeatCovariate.initialize(MAX_STR_UNIT_LENGTH, MAX_REPEAT_LENGTH);

        pcrIndelErrorModelCache = new byte[MAX_REPEAT_LENGTH + 1];

        final double rateFactor = pcrErrorModel == PCR_ERROR_MODEL.AGGRESSIVE ? 2.0 : 3.0;

        for (int iii = 0; iii <= MAX_REPEAT_LENGTH; iii++)
            pcrIndelErrorModelCache[iii] = getErrorModelAdjustedQual(iii, rateFactor);
    }

    protected static byte getErrorModelAdjustedQual(final int repeatLength, final double rateFactor) {
        return (byte) Math.max(MIN_ADJUSTED_QSCORE, MathUtils.fastRound(INITIAL_QSCORE - Math.exp(((double) repeatLength) / (rateFactor * Math.PI)) + 1.0));
    }

    protected void applyPCRErrorModel(final byte[] readBases, final byte[] readInsQuals, final byte[] readDelQuals) {
        if (pcrErrorModel == PCR_ERROR_MODEL.NONE)
            return;

        for (int iii = 1; iii < readBases.length; iii++) {
            final int repeatLength = repeatCovariate.findTandemRepeatUnits(readBases, iii - 1).getSecond();
            readInsQuals[iii - 1] = (byte) Math.min(0xff & readInsQuals[iii - 1], 0xff & pcrIndelErrorModelCache[repeatLength]);
            readDelQuals[iii - 1] = (byte) Math.min(0xff & readDelQuals[iii - 1], 0xff & pcrIndelErrorModelCache[repeatLength]);
        }
    }
}
