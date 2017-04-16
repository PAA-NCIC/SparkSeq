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
package org.ncic.bioinfo.sparkseq.transfer;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import junit.framework.TestCase;
import org.ncic.bioinfo.sparkseq.algorithms.walker.AbstractTestCase;
import org.ncic.bioinfo.sparkseq.algorithms.walker.TestRealignerTargetCreator;
import org.ncic.bioinfo.sparkseq.data.basic.BasicSamRecord;
import org.ncic.bioinfo.sparkseq.data.common.ReadGroupInfo;
import org.ncic.bioinfo.sparkseq.data.common.RefContigInfo;
import org.ncic.bioinfo.sparkseq.data.common.SamHeaderInfo;
import org.ncic.bioinfo.sparkseq.fileio.NormalFileLoader;
import scala.collection.JavaConversions;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Author: wbc
 */
public class TestSAMRecordTransfer extends AbstractTestCase {

    public void testSAMRecordTransfer() {

        RefContigInfo refContigInfo = RefContigInfo.apply(getClass().getResource("/human_g1k_v37.dict").getFile());
        SamHeaderInfo headerInfo = SamHeaderInfo.sortedHeader(refContigInfo, null);
        headerInfo.addReadGroupInfo(ReadGroupInfo.apply("SRR504516", "sample1"));

        SAMFileHeader header = SAMHeaderTransfer.transfer(headerInfo);
        Basic2SAMRecordTransfer basic2SAMRecordTransfer = new Basic2SAMRecordTransfer(header);
        SAMRecord2BasicTransfer samRecord2BasicTransfer = new SAMRecord2BasicTransfer();

        String realPath = TestRealignerTargetCreator.class.getResource("/realigned_reads.sam").getFile();
        List<BasicSamRecord> basicRecords = CollectionConverter.asJavaList(NormalFileLoader.loadSam(realPath, refContigInfo));
        List<SAMRecord> samRecords = basicRecords.stream()
                .map(record -> basic2SAMRecordTransfer.transfer(record))
                .collect(Collectors.toList());

        List<BasicSamRecord> basicRecords2 = samRecords.stream()
                .map(record -> samRecord2BasicTransfer.transfer(record))
                .collect(Collectors.toList());

        for(int i = 0; i < basicRecords.size(); i ++) {
            BasicSamRecord record1 = basicRecords.get(i);
            BasicSamRecord record2 = basicRecords2.get(i);

            assertEquals(record1.toString(), record2.toString());
        }
    }
}
