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
package org.ncic.bioinfo.sparkseq.algorithms.data.sam;

import junit.framework.TestCase;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLoc;
import org.ncic.bioinfo.sparkseq.data.basic.BasicSamRecord;
import org.ncic.bioinfo.sparkseq.data.common.ReadGroupInfo;
import org.ncic.bioinfo.sparkseq.data.common.RefContigInfo;
import org.ncic.bioinfo.sparkseq.data.common.SamHeaderInfo;
import org.ncic.bioinfo.sparkseq.data.partition.SamRecordPartition;
import org.ncic.bioinfo.sparkseq.fileio.NormalFileLoader;
import scala.collection.immutable.List;


/**
 * Author: wbc
 */
public class TestLocusSamTraverser extends TestCase {
    public void testLocusSamTraverser() {
        String filePath = getClass().getResource("/test.sam").getFile();
        RefContigInfo refContigInfo = RefContigInfo.apply(getClass().getResource("/human_g1k_v37.dict").getFile());
        SamHeaderInfo headerInfo = SamHeaderInfo.sortedHeader(refContigInfo, null);
        headerInfo.addReadGroupInfo(ReadGroupInfo.apply("SRR504516", "sample1"));
        List<BasicSamRecord> samRecords = NormalFileLoader.loadSam(filePath, refContigInfo);

        SamRecordPartition samRecordPartition = new SamRecordPartition(1, 0, samRecords, headerInfo);
        SamContentProvider samContentProvider = new SamContentProvider(samRecordPartition);

        GenomeLoc traverseLoc = new GenomeLoc("1", 0, 69773, 70000);
        LocusSamTraverser traverser = new LocusSamTraverser(samContentProvider, traverseLoc);

        assertEquals(true, traverser.hasNext());
        AlignmentContext ac1 = traverser.next();
        assertEquals(0, ac1.basePileup.getBases().length);
        AlignmentContext ac2 = traverser.next();
        assertEquals('A', ac2.basePileup.getBases()[0]);
        assertEquals(69774, ac2.getPosition());
        int count = 69774;
        while (traverser.hasNext()) {
            AlignmentContext ac = traverser.next();
            count ++;
        }
        assertEquals(70000, count);

        traverser.rewind();
        AlignmentContext ac3 = traverser.next();
        assertEquals(0, ac3.basePileup.getBases().length);
        AlignmentContext ac4 = traverser.next();
        assertEquals('A', ac4.basePileup.getBases()[0]);
    }

    public void testLocusSamTraverser2() {
        String filePath = getClass().getResource("/test.sam").getFile();
        RefContigInfo refContigInfo = RefContigInfo.apply(getClass().getResource("/human_g1k_v37.dict").getFile());
        SamHeaderInfo headerInfo = SamHeaderInfo.sortedHeader(refContigInfo, null);
        headerInfo.addReadGroupInfo(ReadGroupInfo.apply("SRR504516", "sample1"));
        List<BasicSamRecord> samRecords = NormalFileLoader.loadSam(filePath, refContigInfo);

        SamRecordPartition samRecordPartition = new SamRecordPartition(1, 0, samRecords, headerInfo);
        SamContentProvider samContentProvider = new SamContentProvider(samRecordPartition);

        GenomeLoc traverseLoc = new GenomeLoc("1", 0, 2090100, 2090101);
        LocusSamTraverser traverser = new LocusSamTraverser(samContentProvider, traverseLoc);

        assertEquals(true, traverser.hasNext());
        AlignmentContext ac1 = traverser.next();
        assertEquals('C', ac1.basePileup.getBases()[0]);
        assertEquals(2090100, ac1.getPosition());
    }
}
