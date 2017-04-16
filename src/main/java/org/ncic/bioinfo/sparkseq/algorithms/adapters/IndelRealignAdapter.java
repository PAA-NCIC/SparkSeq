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

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.SamContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RODContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RODNames;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLoc;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLocParser;
import org.ncic.bioinfo.sparkseq.algorithms.walker.indelrealigner.IndelRealigner;
import org.ncic.bioinfo.sparkseq.algorithms.walker.realignertargetcreator.RealignerTargetCreator;
import org.ncic.bioinfo.sparkseq.data.basic.BasicSamRecord;
import org.ncic.bioinfo.sparkseq.data.common.RefContigInfo;
import org.ncic.bioinfo.sparkseq.data.common.SamHeaderInfo;
import org.ncic.bioinfo.sparkseq.data.partition.FastaPartition;
import org.ncic.bioinfo.sparkseq.data.partition.SamRecordPartition;
import org.ncic.bioinfo.sparkseq.data.partition.VcfRecordPartition;
import org.ncic.bioinfo.sparkseq.transfer.SAMRecord2BasicTransfer;
import org.ncic.bioinfo.sparkseq.transfer.SAMSequenceDictTransfer;
import scala.collection.JavaConversions;

import java.util.ArrayList;
import java.util.List;

/**
 * Author: wbc
 */
public class IndelRealignAdapter {

    public static List<BasicSamRecord> realign(RefContigInfo refContigInfo,
                                               SamRecordPartition samRecordPartition,
                                               FastaPartition refPartition,
                                               List<VcfRecordPartition> rodPartitions) {
        // 数据准备
        SAMSequenceDictionary samSequenceDictionary = SAMSequenceDictTransfer.transfer(refContigInfo);
        GenomeLocParser parser = new GenomeLocParser(samSequenceDictionary);

        SamContentProvider samContentProvider = new SamContentProvider(samRecordPartition);
        RefContentProvider refContentProvider = new RefContentProvider(samSequenceDictionary, refPartition);

        List<RODContentProvider> rodContentProviders = new java.util.ArrayList<>();
        rodPartitions.forEach(
                rodPartition -> rodContentProviders.add(
                        new RODContentProvider(RODNames.KNOWN_ALLELES + rodPartition.key(), rodPartition, parser))
        );

        // 找target interval
        RealignerTargetCreator realignerTargetCreator = new RealignerTargetCreator(
                parser, refContentProvider, samContentProvider, rodContentProviders);
        realignerTargetCreator.run();

        List<GenomeLoc> targetIntervals = realignerTargetCreator.getTargetIntervals();

        // realign reads
        IndelRealigner indelRealigner = new IndelRealigner(
                parser, refContentProvider, samContentProvider, rodContentProviders, targetIntervals);
        indelRealigner.run();

        // 将结果sam record转化为base sam record
        List<SAMRecord> samRecords = indelRealigner.getResultSam();
        List<BasicSamRecord> basicSamRecords = new ArrayList<>();
        SAMRecord2BasicTransfer transfer = new SAMRecord2BasicTransfer();
        samRecords.forEach(record -> {
            basicSamRecords.add(transfer.transfer(record));
        });
        return basicSamRecords;
    }
}
