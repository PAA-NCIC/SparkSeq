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
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceRecord;
import org.ncic.bioinfo.sparkseq.data.common.ReadGroupInfo;
import org.ncic.bioinfo.sparkseq.data.common.SamHeaderInfo;
import scala.collection.JavaConversions;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Author: wbc
 */
public class SAMHeaderTransfer {

    public static SAMFileHeader transfer(SamHeaderInfo samHeaderInfo){
        SAMFileHeader header = new SAMFileHeader();
        // read groups
        List<ReadGroupInfo> readGroupInfoList = CollectionConverter.asJavaList(samHeaderInfo.getReadGroupInfos());
        List<SAMReadGroupRecord> samReadGroupRecords = readGroupInfoList.stream()
                .map(SAMReadGroupRecordTransfer::transfer)
                .collect(Collectors.toList());
        header.setReadGroups(samReadGroupRecords);
        // Sequence records
        // TODO 获取id的方法不通用 fix me
        int contigCount = samHeaderInfo.getRefContigInfo().getContigIds().size();
        for(int id = 0; id < contigCount; id ++) {
            header.addSequence(new SAMSequenceRecord(
                    samHeaderInfo.getRefContigInfo().getName(id),
                    samHeaderInfo.getRefContigInfo().getLength(id)));
        }
        if(samHeaderInfo.sorted()) {
            header.setAttribute("SO", "coordinate");
        }
        return header;
    }
}
