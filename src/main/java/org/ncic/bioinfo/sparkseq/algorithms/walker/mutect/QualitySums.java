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

import java.util.ArrayList;
import java.util.List;

/**
 * Author: wbc
 */
public class QualitySums {
    private int a = 0;
    private int c = 0;
    private int g = 0;
    private int t = 0;
    private int aCounts = 0;
    private int cCounts = 0;
    private int gCounts = 0;
    private int tCounts = 0;

    // used for tracking individual base quality scores, if requested (expensive)
    private boolean enableQualityScoreTracking = false;
    private List<Integer> aQualityScores = new ArrayList<Integer>();
    private List<Integer> cQualityScores = new ArrayList<Integer>();
    private List<Integer> gQualityScores = new ArrayList<Integer>();
    private List<Integer> tQualityScores = new ArrayList<Integer>();

    public QualitySums(boolean enableQualityScoreTracking) {
        this.enableQualityScoreTracking = enableQualityScoreTracking;
    }

    public int getQualitySum(final char base) {
        if (base == 'a' || base == 'A') { return a; }
        if (base == 'c' || base == 'C') { return c; }
        if (base == 'g' || base == 'G') { return g; }
        if (base == 't' || base == 'T') { return t; }
        throw new RuntimeException("Unknown base: " + base);
    }

    public int getCounts(final char base) {
        if (base == 'a' || base == 'A') { return aCounts; }
        if (base == 'c' || base == 'C') { return cCounts; }
        if (base == 'g' || base == 'G') { return gCounts; }
        if (base == 't' || base == 'T') { return tCounts; }
        throw new RuntimeException("Unknown base: " + base);
    }

    public void incrementSum(final char base, final int count, final int qualSum) {
        if (base == 'a' || base == 'A')      { a += qualSum; aCounts+=count; if (enableQualityScoreTracking) aQualityScores.add(qualSum);}
        else if (base == 'c' || base == 'C') { c += qualSum; cCounts+=count; if (enableQualityScoreTracking) cQualityScores.add(qualSum);}
        else if (base == 'g' || base == 'G') { g += qualSum; gCounts+=count; if (enableQualityScoreTracking) gQualityScores.add(qualSum);}
        else if (base == 't' || base == 'T') { t += qualSum; tCounts+=count; if (enableQualityScoreTracking) tQualityScores.add(qualSum);}
        else throw new RuntimeException("Unknown base: " + base);


    }

    public int getOtherQualities(final char base) {
        int total = a + c + g + t;
        if (base == 'a' || base == 'A') { return total-a; }
        else if (base == 'c' || base == 'C') { return total-c; }
        else if (base == 'g' || base == 'G') { return total-g; }
        else if (base == 't' || base == 'T') { return total-t; }
        else throw new RuntimeException("Unknown base: " + base);
    }

    public List<Integer> getBaseQualityScores(final char base) {
        if (base == 'a' || base == 'A') { return aQualityScores; }
        if (base == 'c' || base == 'C') { return cQualityScores; }
        if (base == 'g' || base == 'G') { return gQualityScores; }
        if (base == 't' || base == 'T') { return tQualityScores; }
        throw new RuntimeException("Unknown base: " + base);
    }

    public void reset() {
        a = 0; c = 0; g = 0; t = 0;
        aCounts = 0; cCounts = 0; gCounts = 0; tCounts = 0;
        aQualityScores.clear();
        cQualityScores.clear();
        gQualityScores.clear();
        tQualityScores.clear();
    }
}