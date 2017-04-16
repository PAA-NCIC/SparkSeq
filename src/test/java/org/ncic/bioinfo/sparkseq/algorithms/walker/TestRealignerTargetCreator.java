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
package org.ncic.bioinfo.sparkseq.algorithms.walker;

import htsjdk.samtools.SAMSequenceDictionary;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.SamContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RODContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RODNames;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLoc;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLocParser;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.walker.realignertargetcreator.RealignerTargetCreator;
import org.ncic.bioinfo.sparkseq.data.basic.BasicSamRecord;
import org.ncic.bioinfo.sparkseq.data.common.ReadGroupInfo;
import org.ncic.bioinfo.sparkseq.data.common.RefContigInfo;
import org.ncic.bioinfo.sparkseq.data.common.SamHeaderInfo;
import org.ncic.bioinfo.sparkseq.data.partition.SamRecordPartition;
import org.ncic.bioinfo.sparkseq.data.partition.VcfRecordPartition;
import org.ncic.bioinfo.sparkseq.fileio.NormalFileLoader;
import org.ncic.bioinfo.sparkseq.transfer.SAMSequenceDictTransfer;
import scala.collection.immutable.List;

/**
 * Author: wbc
 */
public class TestRealignerTargetCreator extends AbstractTestCase {

    public void testTargetCreator() {
        RefContigInfo refContigInfo = RefContigInfo.apply(getClass().getResource("/human_g1k_v37.dict").getFile());
        SamHeaderInfo headerInfo = SamHeaderInfo.sortedHeader(refContigInfo, null);
        headerInfo.addReadGroupInfo(ReadGroupInfo.apply("SRR504516", "sample1"));

        SAMSequenceDictionary samSequenceDictionary = SAMSequenceDictTransfer.transfer(refContigInfo);
        GenomeLocParser parser = new GenomeLocParser(samSequenceDictionary);

        String filePath = getClass().getResource("/intervaled.sam").getFile();
        List<BasicSamRecord> samRecords = NormalFileLoader.loadSam(filePath, refContigInfo);
        SamRecordPartition samRecordPartition = new SamRecordPartition(1, 0, samRecords, headerInfo);
        SamContentProvider samContentProvider = new SamContentProvider(samRecordPartition);

        RefContentProvider refContentProvider = getRefContentProvider(samSequenceDictionary);

        VcfRecordPartition vcfRecordPartition = loadVcfPartition("/mills_indel.vcf", "MILLS_INDEL", refContigInfo);
        RODContentProvider rodContentProvider = new RODContentProvider(RODNames.KNOWN_ALLELES, vcfRecordPartition, parser);

        java.util.List<RODContentProvider> rodContentProviders = new java.util.ArrayList<>();
        rodContentProviders.add(rodContentProvider);

        RealignerTargetCreator realignerTargetCreator = new RealignerTargetCreator(
                parser, refContentProvider, samContentProvider, rodContentProviders);

        realignerTargetCreator.run();

        java.util.List<GenomeLoc> result = realignerTargetCreator.getTargetIntervals();

        // 标准答案
        java.util.List<GenomeLoc> standardResult = getRealignerTargetInterval(refContigInfo);
        assertEquals(result.size(), standardResult.size());
        for (int i = 0; i < result.size(); i++) {
            GenomeLoc res1 = result.get(i);
            GenomeLoc res2 = standardResult.get(i);
            assertEquals(res1.getContig(), res2.getContig());
            assertEquals(res1.getStart(), res2.getStart());
            assertEquals(res1.getStop(), res2.getStop());
        }
    }

}
