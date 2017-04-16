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

import java.util.Arrays;
import java.util.List;

/**
 * Author: wbc
 */
public class MuTectStats {

    public static double calculateMAD(double[] dd, double median) {
        double[] dev = new double[dd.length];
        for(int i=0; i<dd.length; i++) {
            dev[i] = Math.abs(dd[i] - median);
        }
        return getMedian(dev);

    }

    public static double getMedian(double[] data) {
        Arrays.sort(data);
        Double result;

        if (data.length % 2 == 1) {
            // If the number of entries in the list is not even.

            // Get the middle value.
            // You must floor the result of the division to drop the
            // remainder.
            result = data[(int) Math.floor(data.length/2) ];

        } else {
            // If the number of entries in the list are even.

            // Get the middle two values and average them.
            Double lowerMiddle = data[data.length/2 ];
            Double upperMiddle = data[data.length/2 - 1 ];
            result = (lowerMiddle + upperMiddle) / 2;
        }

        return result;
    }

    public static double[] convertIntegersToDoubles(List<Integer> integers)
    {
        double[] ret = new double[integers.size()];
        for (int i=0; i < ret.length; i++)
        {
            ret[i] = integers.get(i);
        }
        return ret;
    }
}
