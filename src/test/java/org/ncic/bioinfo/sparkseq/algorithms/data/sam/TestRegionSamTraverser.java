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
import org.ncic.bioinfo.sparkseq.exceptions.GATKException;
import org.ncic.bioinfo.sparkseq.fileio.NormalFileLoader;
import scala.collection.immutable.List;

/**
 * Author: wbc
 */
public class TestRegionSamTraverser extends TestCase {

    public void testRegionSamTraverser() {
        String filePath = getClass().getResource("/test.sam").getFile();
        RefContigInfo refContigInfo = RefContigInfo.apply(getClass().getResource("/human_g1k_v37.dict").getFile());
        SamHeaderInfo headerInfo = SamHeaderInfo.sortedHeader(refContigInfo, null);
        headerInfo.addReadGroupInfo(ReadGroupInfo.apply("SRR504516", "sample1"));
        List<BasicSamRecord> samRecords = NormalFileLoader.loadSam(filePath, refContigInfo);

        SamRecordPartition samRecordPartition = new SamRecordPartition(1, 0, samRecords, headerInfo);
        SamContentProvider samContentProvider = new SamContentProvider(samRecordPartition);
        RegionSamTraverser traverser = new RegionSamTraverser(samContentProvider);

        GenomeLoc locus = new GenomeLoc("1", 0, 1, 1000);
        java.util.List<GATKSAMRecord> res1 = traverser.getOverlappedReads(locus);
        assertEquals(res1.size(), 0);

        GenomeLoc locus2 = new GenomeLoc("1", 0, 1000, 70000);
        java.util.List<GATKSAMRecord> res2 = traverser.getOverlappedReads(locus2);
        assertEquals(res2.size(), 2);

        GenomeLoc locus3 = new GenomeLoc("1", 0, 1000, 2090293);
        java.util.List<GATKSAMRecord> res3 = traverser.getOverlappedReads(locus3);
        assertEquals(res3.size(), 4);

        GenomeLoc locus4 = new GenomeLoc("1", 0, 2090293, 2090293);
        java.util.List<GATKSAMRecord> res4 = traverser.getOverlappedReads(locus4);
        assertEquals(res4.size(), 1);

        try {
            GenomeLoc locus5 = new GenomeLoc("1", 0, 1000, 2090293);
            java.util.List<GATKSAMRecord> res5 = traverser.getOverlappedReads(locus5);
        } catch (GATKException ex) {
            System.out.println("OK");
        }

    }
}
