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

import java.util.Arrays;

/**
 * Author: wbc
 */
public final class ExactACcounts {
    private final int[] counts;
    private int hashcode = -1;

    public ExactACcounts(final int[] counts) {
        this.counts = counts;
    }

    public int[] getCounts() {
        return counts;
    }

    @Override
    public boolean equals(Object obj) {
        return (obj instanceof ExactACcounts) && Arrays.equals(getCounts(), ((ExactACcounts) obj).getCounts());
    }

    @Override
    public int hashCode() {
        if ( hashcode == -1 )
            hashcode = Arrays.hashCode(getCounts());
        return hashcode;
    }

    @Override
    public String toString() {
        StringBuffer sb = new StringBuffer();
        sb.append(getCounts()[0]);
        for ( int i = 1; i < getCounts().length; i++ ) {
            sb.append("/");
            sb.append(getCounts()[i]);
        }
        return sb.toString();
    }
}
