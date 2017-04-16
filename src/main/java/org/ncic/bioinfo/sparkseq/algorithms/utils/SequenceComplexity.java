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

/**
 * Author: wbc
 */
public class SequenceComplexity {

    /**
     * Indicates what positions in a base sequence is found in homopolymers or STR repeat.
     * <p>
     * <p>
     * The result is an boolean array with as many positions as the input base array.
     * </p>
     * <p>
     * Each entry the result makes reference to the base at the same position in the input, where {@code true}
     * means that it forms part of a repeat.
     * </p>
     *
     * @param bases                  the input bases.
     * @param maxRepeatUnitLength    what is the largest repeat unit to consider.
     * @param minRepeatLengthInUnits what is minimum length of a repeat in units to consider it significantly long. Shorter
     *                               repeats won't be considered as such.
     * @return never {@code null} but an array with the same length as the reference haplotype.
     */
    public static boolean[] findBasesInShortUnitRepeats(final byte[] bases, final int maxRepeatUnitLength,
                                                        final int minRepeatLengthInUnits) {
        final boolean[] result = new boolean[bases.length];
        final int[] repeatAbsoluteLengthCount = new int[maxRepeatUnitLength];
        for (int i = 0; i < maxRepeatUnitLength; i++)
            repeatAbsoluteLengthCount[i] = i + 1;
        for (int i = 0; i < bases.length; i++)
            for (int j = 1; j <= maxRepeatUnitLength; j++) {
                final int prevI = i - j;
                if (prevI < 0) continue;
                if (bases[prevI] == bases[i]) // repeat continuation.
                    repeatAbsoluteLengthCount[j - 1]++;
                else if (minRepeatLengthInUnits <= (repeatAbsoluteLengthCount[j - 1] / j)) {  // end of a long enough repeat.
                    for (int k = i - repeatAbsoluteLengthCount[j - 1]; k < i; k++)
                        result[k] = true;
                    repeatAbsoluteLengthCount[j - 1] = j;
                } else { // end of not long enough repeat.
                    repeatAbsoluteLengthCount[j - 1] = j;
                }
            }
        // Do the marking for the last position in bases.
        for (int j = 1; j <= maxRepeatUnitLength; j++)
            if (minRepeatLengthInUnits <= (repeatAbsoluteLengthCount[j - 1] / j))
                for (int k = bases.length - repeatAbsoluteLengthCount[j - 1]; k < bases.length; k++)
                    result[k] = true;
        return result;
    }

}
