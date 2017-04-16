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

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.GATKSAMRecord;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.SamContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RODContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RODNames;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.header.StandardWESVCFHeader;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLoc;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLocParser;
import org.ncic.bioinfo.sparkseq.algorithms.walker.mutect.Mutect;
import org.ncic.bioinfo.sparkseq.data.basic.VcfRecord;
import org.ncic.bioinfo.sparkseq.data.common.ReadGroupInfo;
import org.ncic.bioinfo.sparkseq.data.common.RefContigInfo;
import org.ncic.bioinfo.sparkseq.data.common.SamHeaderInfo;
import org.ncic.bioinfo.sparkseq.data.partition.FastaPartition;
import org.ncic.bioinfo.sparkseq.data.partition.VcfRecordPartition;
import org.ncic.bioinfo.sparkseq.transfer.SAMHeaderTransfer;
import org.ncic.bioinfo.sparkseq.transfer.SAMSequenceDictTransfer;
import org.ncic.bioinfo.sparkseq.transfer.VC2VcfRecordTransfer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Author: wbc
 */
public class TestMutect extends AbstractTestCase {

    public void testMutect() {
        RefContigInfo refContigInfo = RefContigInfo.apply(getClass().getResource("/human_g1k_v37.dict").getFile());
        SamHeaderInfo headerInfo = SamHeaderInfo.sortedHeader(refContigInfo, null);
        headerInfo.addReadGroupInfo(ReadGroupInfo.apply("SRR504516", "sample1"));

        SAMSequenceDictionary samSequenceDictionary = SAMSequenceDictTransfer.transfer(refContigInfo);
        GenomeLocParser parser = new GenomeLocParser(samSequenceDictionary);

        List<SAMRecord> normalReads = getNormalReads(headerInfo);
        List<GATKSAMRecord> normalGATKRecords = new ArrayList<>();
        for (SAMRecord record : normalReads) {
            normalGATKRecords.add(new GATKSAMRecord(record));
        }
        SamContentProvider normalContentProvider = new SamContentProvider(normalGATKRecords, SAMHeaderTransfer.transfer(headerInfo));

        List<SAMRecord> tumorReads = getTumorReads(headerInfo);
        List<GATKSAMRecord> tumorGATKRecords = new ArrayList<>();
        for (SAMRecord record : tumorReads) {
            tumorGATKRecords.add(new GATKSAMRecord(record));
        }
        SamContentProvider tumorContentProvider = new SamContentProvider(tumorGATKRecords, SAMHeaderTransfer.transfer(headerInfo));

        RefContentProvider refContentProvider = getMutectRefContentProvider(samSequenceDictionary);

        VcfRecordPartition vcfRecordPartition = loadVcfPartition("/mutect/head_dbsnp.vcf", RODNames.DBSNP, refContigInfo);
        RODContentProvider rodContentProvider = new RODContentProvider(RODNames.DBSNP, vcfRecordPartition, parser);
        java.util.List<RODContentProvider> rodContentProviders = new java.util.ArrayList<>();
        rodContentProviders.add(rodContentProvider);

        List<GenomeLoc> intervals = getIntervals(refContigInfo, "/mutect/0_panel.intervals");

        Mutect mutect = new Mutect(parser, refContentProvider, tumorContentProvider, normalContentProvider, rodContentProviders, intervals);
        // 设置mutect参数
        mutect.MTAC.FRACTION_CONTAMINATION = 0;
        mutect.run();
        List<VariantContext> vcfs = mutect.getResultVCFRecords();
        List<String> outs = mutect.getResultVCFOutInfos();

        VCFHeader header = StandardWESVCFHeader.getHeader();
        VC2VcfRecordTransfer transfer = new VC2VcfRecordTransfer(header, refContigInfo);
        List<VcfRecord> standardResult = getMutectVcf(refContigInfo);

        List<VcfRecord> transferedResult = vcfs.stream()
                .map(record -> transfer.transfer(record))
                .collect(Collectors.toList());
        assertEquals(transferedResult.size(), standardResult.size());
        for (int i = 0; i < transferedResult.size(); i++) {
            VcfRecord record1 = standardResult.get(i);
            VcfRecord record2 = transferedResult.get(i);
            assertEquals(record1.toString(), record2.toString());
        }

        List<String> resultTxt = getMutectOutTxt();
        assertEquals(resultTxt.size(), outs.size());
        for(int i = 0; i < resultTxt.size(); i ++) {
            String line1 = resultTxt.get(i);
            String line2 = outs.get(i);
            assertEquals(line1, line2);
        }
    }

    protected static java.util.List<SAMRecord> getNormalReads(SamHeaderInfo headerInfo) {
        return getReads("/mutect/dedup_Blood.sam", headerInfo);
    }

    private static java.util.List<SAMRecord> getTumorReads(SamHeaderInfo headerInfo) {
        return getReads("/mutect/dedup_Tumor.sam", headerInfo);
    }

    private static RefContentProvider getMutectRefContentProvider(SAMSequenceDictionary samSequenceDictionary) {
        FastaPartition fastaPartition = new FastaPartition(1, 0, "1", readRefContentValue("/mutect/refContent"),
                1, 2501000, 1, 2500000, 1, 2497000);
        return new RefContentProvider(samSequenceDictionary, fastaPartition);
    }

    private static java.util.List<VcfRecord> getMutectVcf(RefContigInfo refContigInfo) {
        return getGvcfValue("/mutect/out.vcf", refContigInfo);
    }

    private static java.util.List<String> getMutectOutTxt() {
        List<String> result = new ArrayList<>();
        String filePath = TestRealignerTargetCreator.class.getResource("/mutect/out.txt").getFile();
        try (BufferedReader reader = new BufferedReader(new FileReader(new File(filePath)))) {
            String line = reader.readLine();
            while (line != null) {
                if (line.length() > 0)
                    result.add(line);
                line = reader.readLine();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return result;
    }

}
