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
package org.ncic.bioinfo.sparkseq.algorithms.utils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.ReferenceContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.AlignmentContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.PileupElement;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.ReadBackedPileup;

/**
 * Author: wbc
 */
public class CGAAlignmentUtils {

    /** Returns the number of mismatches in the pileup element within the given reference context.
     *
     * @param p       the pileup element
     * @param ref     the reference context
     * @param ignoreTargetSite     if true, ignore mismatches at the target locus (i.e. the center of the window)
     * @param qualitySumInsteadOfMismatchCount if true, return the quality score sum of the mismatches rather than the count
     * @return the number of mismatches
     */
    public static int mismatchesInRefWindow(PileupElement p, ReferenceContext ref,
                                            RefContentProvider refContentProvider,
                                            boolean ignoreTargetSite, boolean qualitySumInsteadOfMismatchCount) {
        int sum = 0;

        ReferenceContext referenceContext = refContentProvider.getReferenceContext(ref.getLocus(), 150);

        int windowStart = referenceContext.getLocus().getStart();
        int windowStop = referenceContext.getLocus().getStop();
        byte[] refBases = referenceContext.getBases();
        byte[] readBases = p.getRead().getReadBases();
        byte[] readQualities = p.getRead().getBaseQualities();
        Cigar c = p.getRead().getCigar();

        int readIndex = 0;
        int currentPos = p.getRead().getAlignmentStart();
        int refIndex = Math.max(0, currentPos - windowStart);

        for (int i = 0 ; i < c.numCigarElements() ; i++) {
            CigarElement ce = c.getCigarElement(i);
            int cigarElementLength = ce.getLength();
            switch ( ce.getOperator() ) {
                case M:
                    for (int j = 0; j < cigarElementLength; j++, readIndex++, currentPos++) {
                        // are we past the ref window?
                        if ( currentPos > windowStop )
                            break;

                        // are we before the ref window?
                        if ( currentPos < windowStart )
                            continue;

                        byte refChr = refBases[refIndex++];

                        // do we need to skip the target site?
                        if ( ignoreTargetSite && ref.getLocus().getStart() == currentPos )
                            continue;

                        byte readChr = readBases[readIndex];
                        if ( readChr != refChr )
                            sum += (qualitySumInsteadOfMismatchCount) ? readQualities[readIndex] : 1;
                    }
                    break;
                case I:
                case S:
                    readIndex += cigarElementLength;
                    break;
                case D:
                case N:
                    currentPos += cigarElementLength;
                    if ( currentPos > windowStart )
                        refIndex += Math.min(cigarElementLength, currentPos - windowStart);
                    break;
                case H:
                case P:
                    break;
            }
        }

        return sum;
    }

}
