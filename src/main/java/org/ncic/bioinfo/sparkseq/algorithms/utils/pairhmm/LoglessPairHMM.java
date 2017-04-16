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
package org.ncic.bioinfo.sparkseq.algorithms.utils.pairhmm;

import org.ncic.bioinfo.sparkseq.algorithms.utils.QualityUtils;
import static org.ncic.bioinfo.sparkseq.algorithms.utils.pairhmm.PairHMMModel.*;

/**
 * Author: wbc
 */
public class LoglessPairHMM extends N2MemoryPairHMM {
    protected static final double INITIAL_CONDITION = Math.pow(2, 1020);
    protected static final double INITIAL_CONDITION_LOG10 = Math.log10(INITIAL_CONDITION);

    // we divide e by 3 because the observed base could have come from any of the non-observed alleles
    protected static final double TRISTATE_CORRECTION = 3.0;

    /**
     * {@inheritDoc}
     */
    @Override
    public double subComputeReadLikelihoodGivenHaplotypeLog10( final byte[] haplotypeBases,
                                                               final byte[] readBases,
                                                               final byte[] readQuals,
                                                               final byte[] insertionGOP,
                                                               final byte[] deletionGOP,
                                                               final byte[] overallGCP,
                                                               final int hapStartIndex,
                                                               final boolean recacheReadValues,
                                                               final int nextHapStartIndex) {

        if (previousHaplotypeBases == null || previousHaplotypeBases.length != haplotypeBases.length) {
            final double initialValue = INITIAL_CONDITION / haplotypeBases.length;
            // set the initial value (free deletions in the beginning) for the first row in the deletion matrix
            for( int j = 0; j < paddedHaplotypeLength; j++ ) {
                deletionMatrix[0][j] = initialValue;
            }
        }

        if ( ! constantsAreInitialized || recacheReadValues ) {
            initializeProbabilities(transition, insertionGOP, deletionGOP, overallGCP);

            // note that we initialized the constants
            constantsAreInitialized = true;
        }

        initializePriors(haplotypeBases, readBases, readQuals, hapStartIndex);

        for (int i = 1; i < paddedReadLength; i++) {
            // +1 here is because hapStartIndex is 0-based, but our matrices are 1 based
            for (int j = hapStartIndex+1; j < paddedHaplotypeLength; j++) {
                //Inlined the code from updateCell - helps JIT to detect hotspots and produce good native code
                matchMatrix[i][j] = prior[i][j] * ( matchMatrix[i - 1][j - 1] * transition[i][matchToMatch] +
                        insertionMatrix[i - 1][j - 1] * transition[i][indelToMatch] +
                        deletionMatrix[i - 1][j - 1] * transition[i][indelToMatch] );
                insertionMatrix[i][j] = matchMatrix[i - 1][j] * transition[i][matchToInsertion] + insertionMatrix[i - 1][j] * transition[i][insertionToInsertion];
                deletionMatrix[i][j] = matchMatrix[i][j - 1] * transition[i][matchToDeletion] + deletionMatrix[i][j - 1] * transition[i][deletionToDeletion];
            }
        }

        // final probability is the log10 sum of the last element in the Match and Insertion state arrays
        // this way we ignore all paths that ended in deletions! (huge)
        // but we have to sum all the paths ending in the M and I matrices, because they're no longer extended.
        final int endI = paddedReadLength - 1;
        double finalSumProbabilities = 0.0;
        for (int j = 1; j < paddedHaplotypeLength; j++) {
            finalSumProbabilities += matchMatrix[endI][j] + insertionMatrix[endI][j];
        }
        return Math.log10(finalSumProbabilities) - INITIAL_CONDITION_LOG10;
    }

    /**
     * Initializes the matrix that holds all the constants related to the editing
     * distance between the read and the haplotype.
     *
     * @param haplotypeBases the bases of the haplotype
     * @param readBases      the bases of the read
     * @param readQuals      the base quality scores of the read
     * @param startIndex     where to start updating the distanceMatrix (in case this read is similar to the previous read)
     */
    protected void initializePriors(final byte[] haplotypeBases, final byte[] readBases, final byte[] readQuals, final int startIndex) {

        // initialize the pBaseReadLog10 matrix for all combinations of read x haplotype bases
        // Abusing the fact that java initializes arrays with 0.0, so no need to fill in rows and columns below 2.

        for (int i = 0; i < readBases.length; i++) {
            final byte x = readBases[i];
            final byte qual = readQuals[i];
            for (int j = startIndex; j < haplotypeBases.length; j++) {
                final byte y = haplotypeBases[j];
                prior[i+1][j+1] = ( x == y || x == (byte) 'N' || y == (byte) 'N' ?
                        QualityUtils.qualToProb(qual) : (QualityUtils.qualToErrorProb(qual) / (doNotUseTristateCorrection ? 1.0 : TRISTATE_CORRECTION)) );
            }
        }
    }

    /**
     * Initializes the matrix that holds all the constants related to quality scores.
     *
     * @param insertionGOP   insertion quality scores of the read
     * @param deletionGOP    deletion quality scores of the read
     * @param overallGCP     overall gap continuation penalty
     */
    protected static void initializeProbabilities(final double[][] transition, final byte[] insertionGOP, final byte[] deletionGOP, final byte[] overallGCP) {
        PairHMMModel.qualToTransProbs(transition,insertionGOP,deletionGOP,overallGCP);
    }

    /**
     * Updates a cell in the HMM matrix
     *
     * The read and haplotype indices are offset by one because the state arrays have an extra column to hold the
     * initial conditions

     * @param indI             row index in the matrices to update
     * @param indJ             column index in the matrices to update
     * @param prior            the likelihood editing distance matrix for the read x haplotype
     * @param transition        an array with the six transition relevant to this location
     */
    protected void updateCell( final int indI, final int indJ, final double prior, final double[] transition) {

        matchMatrix[indI][indJ] = prior * ( matchMatrix[indI - 1][indJ - 1] * transition[matchToMatch] +
                insertionMatrix[indI - 1][indJ - 1] * transition[indelToMatch] +
                deletionMatrix[indI - 1][indJ - 1] * transition[indelToMatch] );
        insertionMatrix[indI][indJ] = matchMatrix[indI - 1][indJ] * transition[matchToInsertion] + insertionMatrix[indI - 1][indJ] * transition[insertionToInsertion];
        deletionMatrix[indI][indJ] = matchMatrix[indI][indJ - 1] * transition[matchToDeletion] + deletionMatrix[indI][indJ - 1] * transition[deletionToDeletion];
    }
}
