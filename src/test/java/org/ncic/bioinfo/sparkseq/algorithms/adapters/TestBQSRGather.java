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
package org.ncic.bioinfo.sparkseq.algorithms.adapters;

import junit.framework.TestCase;
import org.ncic.bioinfo.sparkseq.algorithms.walker.TestRealignerTargetCreator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

/**
 * Author: wbc
 */
public class TestBQSRGather extends TestCase {

    public void testBQSRGather() {
        List<List<String>> bqsrSources = new ArrayList<>();
        bqsrSources.add(getRecalTableLines("/bqsrtables/0_bqsr.table"));
        bqsrSources.add(getRecalTableLines("/bqsrtables/1_bqsr.table"));
        bqsrSources.add(getRecalTableLines("/bqsrtables/2_bqsr.table"));
        bqsrSources.add(getRecalTableLines("/bqsrtables/3_bqsr.table"));

        BQSRTableGather bqsrTableGather = new BQSRTableGather();

        List<String> mergedTableLines = bqsrTableGather.gatherBQSRTables(bqsrSources);
        List<String> standardTableLines = getRecalTableLines("/bqsrtables/merged.table");
        for (int i = 0; i < mergedTableLines.size(); i++) {
            assertEquals(mergedTableLines.get(i), standardTableLines.get(i));
        }
    }

    public void testParallelBQSRGather() {
        List<List<String>> bqsrSources = new ArrayList<>();
        bqsrSources.add(getRecalTableLines("/bqsrtables/0_bqsr.table"));
        bqsrSources.add(getRecalTableLines("/bqsrtables/1_bqsr.table"));
        bqsrSources.add(getRecalTableLines("/bqsrtables/2_bqsr.table"));
        bqsrSources.add(getRecalTableLines("/bqsrtables/3_bqsr.table"));

        BQSRTableGather bqsrTableGather = new BQSRTableGather();

        List<String> mergedTableLines = bqsrTableGather.gatherBQSRTablesInParallel(bqsrSources, 12, Integer.MAX_VALUE);
        List<String> standardTableLines = getRecalTableLines("/bqsrtables/merged.table");
        for (int i = 0; i < mergedTableLines.size(); i++) {
            assertEquals(mergedTableLines.get(i), standardTableLines.get(i));
        }
    }

    private static List<String> getRecalTableLines(String filePath) {
        java.util.List<String> res = new ArrayList<>();
        String realignedSamPath = TestRealignerTargetCreator.class.getResource(filePath).getFile();
        try (BufferedReader reader = new BufferedReader(new FileReader(new File(realignedSamPath)))) {
            String line = reader.readLine();
            while (line != null) {
                res.add(line);
                line = reader.readLine();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        return res;
    }
}
