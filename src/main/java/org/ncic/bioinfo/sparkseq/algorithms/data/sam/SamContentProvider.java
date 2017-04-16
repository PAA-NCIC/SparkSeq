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

import htsjdk.samtools.SAMFileHeader;
import org.ncic.bioinfo.sparkseq.data.basic.BasicSamRecord;
import org.ncic.bioinfo.sparkseq.data.common.SamHeaderInfo;
import org.ncic.bioinfo.sparkseq.data.partition.SamRecordPartition;
import org.ncic.bioinfo.sparkseq.transfer.Basic2SAMRecordTransfer;
import org.ncic.bioinfo.sparkseq.transfer.SAMHeaderTransfer;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Author: wbc
 */
public class SamContentProvider {

    private final static int FAKE_CONTIG_ID = 255;  // 与scala中SamRecordConst中定义冗余

    private final SAMFileHeader samFileHeader;
    private final ArrayList<GATKSAMRecord> gatksamRecords;

    public SamContentProvider(SamRecordPartition samRecordPartition) {
        SamHeaderInfo headerInfo = samRecordPartition.samHeaderInfo();
        SAMFileHeader header = SAMHeaderTransfer.transfer(headerInfo);
        samFileHeader = header;
        gatksamRecords = new ArrayList<>();

        Basic2SAMRecordTransfer transfer = new Basic2SAMRecordTransfer(header);
        scala.collection.Iterable<BasicSamRecord> iterable = samRecordPartition.records();
        scala.collection.Iterator<BasicSamRecord> iter = iterable.iterator();
        while (iter.hasNext()) {
            BasicSamRecord record = iter.next();
            if (record.contigId() == samRecordPartition.contigId()) {   //滤掉unmapped，保证所有read在一个contig上
                gatksamRecords.add(new GATKSAMRecord(transfer.transfer(record)));
            }
        }

        //为samRecord进行排序
        Collections.sort(gatksamRecords,
                (record1, record2) -> record1.getAlignmentStart() - record2.getAlignmentStart());
    }

    public SamContentProvider(List<GATKSAMRecord> gatksamRecords,
                              SAMFileHeader samFileHeader) {
        this(gatksamRecords, samFileHeader, true);
    }

    public SamContentProvider(List<GATKSAMRecord> gatksamRecords,
                              SAMFileHeader samFileHeader, boolean needSort) {
        this.gatksamRecords = new ArrayList<>();
        this.samFileHeader = samFileHeader;

        this.gatksamRecords.addAll(gatksamRecords);

        //为samRecord进行排序
        if (needSort) {
            Collections.sort(gatksamRecords,
                    (record1, record2) -> record1.getAlignmentStart() - record2.getAlignmentStart());
        }
    }

    public SAMFileHeader getSamFileHeader() {
        return samFileHeader;
    }

    public ArrayList<GATKSAMRecord> getGatksamRecords() {
        return gatksamRecords;
    }
}
