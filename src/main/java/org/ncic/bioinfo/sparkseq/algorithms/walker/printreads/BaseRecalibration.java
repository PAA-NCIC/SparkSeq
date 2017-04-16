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
package org.ncic.bioinfo.sparkseq.algorithms.walker.printreads;

import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import org.apache.log4j.Logger;
import org.ncic.bioinfo.sparkseq.algorithms.utils.EventType;
import org.ncic.bioinfo.sparkseq.algorithms.utils.MathUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.QualityUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.RecalUtils;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.GATKSAMRecord;
import org.ncic.bioinfo.sparkseq.algorithms.utils.reports.GATKReport;
import org.ncic.bioinfo.sparkseq.algorithms.walker.baserecalibrator.QuantizationInfo;
import org.ncic.bioinfo.sparkseq.algorithms.walker.baserecalibrator.ReadCovariates;
import org.ncic.bioinfo.sparkseq.algorithms.walker.baserecalibrator.RecalDatum;
import org.ncic.bioinfo.sparkseq.algorithms.walker.baserecalibrator.RecalibrationTables;
import org.ncic.bioinfo.sparkseq.algorithms.walker.baserecalibrator.covariate.Covariate;
import org.ncic.bioinfo.sparkseq.exceptions.UserException;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Author: wbc
 */
public class BaseRecalibration {
    private static Logger logger = Logger.getLogger(BaseRecalibration.class);
    private final static boolean TEST_CACHING = false;

    private final QuantizationInfo quantizationInfo; // histogram containing the map for qual quantization (calculated after recalibration is done)
    private final RecalibrationTables recalibrationTables;
    private final Covariate[] requestedCovariates; // list of all covariates to be used in this calculation

    private final boolean disableIndelQuals;
    private final int preserveQLessThan;
    private final double globalQScorePrior;
    private final boolean emitOriginalQuals;

    /**
     * Constructor using a GATK Report file
     *
     * @param report             a GATK Report  containing the recalibration information
     * @param quantizationLevels number of bins to quantize the quality scores
     * @param disableIndelQuals  if true, do not emit base indel qualities
     * @param preserveQLessThan  preserve quality scores less than this value
     */
    public BaseRecalibration(final GATKReport report, final int quantizationLevels, final boolean disableIndelQuals, final int preserveQLessThan, final boolean emitOriginalQuals, final double globalQScorePrior) {
        RecalibrationReport recalibrationReport = new RecalibrationReport(report);

        recalibrationTables = recalibrationReport.getRecalibrationTables();
        requestedCovariates = recalibrationReport.getRequestedCovariates();
        quantizationInfo = recalibrationReport.getQuantizationInfo();
        if (quantizationLevels == 0) // quantizationLevels == 0 means no quantization, preserve the quality scores
            quantizationInfo.noQuantization();
        else if (quantizationLevels > 0 && quantizationLevels != quantizationInfo.getQuantizationLevels()) // any other positive value means, we want a different quantization than the one pre-calculated in the recalibration report. Negative values mean the user did not provide a quantization argument, and just wants to use what's in the report.
            quantizationInfo.quantizeQualityScores(quantizationLevels);

        this.disableIndelQuals = disableIndelQuals;
        this.preserveQLessThan = preserveQLessThan;
        this.globalQScorePrior = globalQScorePrior;
        this.emitOriginalQuals = emitOriginalQuals;
    }

    /**
     * Recalibrates the base qualities of a read
     * <p>
     * It updates the base qualities of the read with the new recalibrated qualities (for all event types)
     * <p>
     * Implements a serial recalibration of the reads using the combinational table.
     * First, we perform a positional recalibration, and then a subsequent dinuc correction.
     * <p>
     * Given the full recalibration table, we perform the following preprocessing steps:
     * <p>
     * - calculate the global quality score shift across all data [DeltaQ]
     * - calculate for each of cycle and dinuc the shift of the quality scores relative to the global shift
     * -- i.e., DeltaQ(dinuc) = Sum(pos) Sum(Qual) Qempirical(pos, qual, dinuc) - Qreported(pos, qual, dinuc) / Npos * Nqual
     * - The final shift equation is:
     * <p>
     * Qrecal = Qreported + DeltaQ + DeltaQ(pos) + DeltaQ(dinuc) + DeltaQ( ... any other covariate ... )
     *
     * @param read the read to recalibrate
     */
    public void recalibrateRead(final GATKSAMRecord read) {
        if (emitOriginalQuals && read.getAttribute(SAMTag.OQ.name()) == null) { // Save the old qualities if the tag isn't already taken in the read
            try {
                read.setAttribute(SAMTag.OQ.name(), SAMUtils.phredToFastq(read.getBaseQualities()));
            } catch (IllegalArgumentException e) {
                throw new UserException.MalformedBAM(read, "illegal base quality encountered; " + e.getMessage());
            }
        }

        final ReadCovariates readCovariates = RecalUtils.computeCovariates(read, requestedCovariates);
        final int readLength = read.getReadLength();

        for (final EventType errorModel : EventType.values()) { // recalibrate all three quality strings
            if (disableIndelQuals && errorModel != EventType.BASE_SUBSTITUTION) {
                read.setBaseQualities(null, errorModel);
                continue;
            }

            final byte[] quals = read.getBaseQualities(errorModel);

            // get the keyset for this base using the error model
            final int[][] fullReadKeySet = readCovariates.getKeySet(errorModel);

            // the rg key is constant over the whole read, the global deltaQ is too
            final int rgKey = fullReadKeySet[0][0];
            final RecalDatum empiricalQualRG = recalibrationTables.getReadGroupTable().get(rgKey, errorModel.ordinal());

            if (empiricalQualRG != null) {
                final double epsilon = (globalQScorePrior > 0.0 && errorModel.equals(EventType.BASE_SUBSTITUTION) ? globalQScorePrior : empiricalQualRG.getEstimatedQReported());

                for (int offset = 0; offset < readLength; offset++) { // recalibrate all bases in the read
                    final byte origQual = quals[offset];

                    // only recalibrate usable qualities (the original quality will come from the instrument -- reported quality)
                    if (origQual >= preserveQLessThan) {
                        // get the keyset for this base using the error model
                        final int[] keySet = fullReadKeySet[offset];
                        final RecalDatum empiricalQualQS = recalibrationTables.getQualityScoreTable().get(keySet[0], keySet[1], errorModel.ordinal());
                        final List<RecalDatum> empiricalQualCovs = new ArrayList<RecalDatum>();
                        for (int i = 2; i < requestedCovariates.length; i++) {
                            if (keySet[i] < 0) {
                                continue;
                            }
                            empiricalQualCovs.add(recalibrationTables.getTable(i).get(keySet[0], keySet[1], keySet[i], errorModel.ordinal()));
                        }

                        double recalibratedQualDouble = hierarchicalBayesianQualityEstimate(epsilon, empiricalQualRG, empiricalQualQS, empiricalQualCovs);

                        // recalibrated quality is bound between 1 and MAX_QUAL
                        final byte recalibratedQual = QualityUtils.boundQual(MathUtils.fastRound(recalibratedQualDouble), RecalDatum.MAX_RECALIBRATED_Q_SCORE);

                        // return the quantized version of the recalibrated quality
                        final byte recalibratedQualityScore = quantizationInfo.getQuantizedQuals().get(recalibratedQual);

                        quals[offset] = recalibratedQualityScore;
                    }
                }
            }

            // finally update the base qualities in the read
            read.setBaseQualities(quals, errorModel);
        }
    }

    protected static double hierarchicalBayesianQualityEstimate(final double epsilon, final RecalDatum empiricalQualRG, final RecalDatum empiricalQualQS, final List<RecalDatum> empiricalQualCovs) {
        final double globalDeltaQ = (empiricalQualRG == null ? 0.0 : empiricalQualRG.getEmpiricalQuality(epsilon) - epsilon);
        final double deltaQReported = (empiricalQualQS == null ? 0.0 : empiricalQualQS.getEmpiricalQuality(globalDeltaQ + epsilon) - (globalDeltaQ + epsilon));
        double deltaQCovariates = 0.0;
        for (final RecalDatum empiricalQualCov : empiricalQualCovs) {
            deltaQCovariates += (empiricalQualCov == null ? 0.0 : empiricalQualCov.getEmpiricalQuality(deltaQReported + globalDeltaQ + epsilon) - (deltaQReported + globalDeltaQ + epsilon));
        }

        return epsilon + globalDeltaQ + deltaQReported + deltaQCovariates;
    }
}
