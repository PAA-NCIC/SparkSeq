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
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.GATKSAMRecord;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.SamContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RODContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RODNames;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLocParser;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.ActiveRegionFinder;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.ActiveRegionMapData;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.HaplotypeCaller;
import org.ncic.bioinfo.sparkseq.data.basic.VcfRecord;
import org.ncic.bioinfo.sparkseq.data.common.ReadGroupInfo;
import org.ncic.bioinfo.sparkseq.data.common.RefContigInfo;
import org.ncic.bioinfo.sparkseq.data.common.SamHeaderInfo;
import org.ncic.bioinfo.sparkseq.data.partition.VcfRecordPartition;
import org.ncic.bioinfo.sparkseq.transfer.SAMHeaderTransfer;
import org.ncic.bioinfo.sparkseq.transfer.SAMSequenceDictTransfer;

import java.util.ArrayList;
import java.util.List;

/**
 * Author: wbc
 */
public class TestHaplotypeCaller extends AbstractTestCase {
/*
    public void testHaplotypeCaller() {
        RefContigInfo refContigInfo = RefContigInfo.apply(getClass().getResource("/human_g1k_v37.dict").getFile());
        SamHeaderInfo headerInfo = SamHeaderInfo.sortedHeader(refContigInfo, null);
        headerInfo.addReadGroupInfo(ReadGroupInfo.apply("SRR504516", "sample1"));

        SAMSequenceDictionary samSequenceDictionary = SAMSequenceDictTransfer.transfer(refContigInfo);
        GenomeLocParser parser = new GenomeLocParser(samSequenceDictionary);

        List<SAMRecord> recaledReads = getRecaledReads(headerInfo);
        List<GATKSAMRecord> realignedGATKRecords = new ArrayList<>();
        for (SAMRecord record : recaledReads) {
            realignedGATKRecords.add(new GATKSAMRecord(record));
        }
        SamContentProvider samContentProvider = new SamContentProvider(realignedGATKRecords, SAMHeaderTransfer.transfer(headerInfo));

        RefContentProvider refContentProvider = getRefContentProvider(samSequenceDictionary);

        java.util.List<RODContentProvider> rodContentProviders = new java.util.ArrayList<>();

        ActiveRegionFinder activeRegionFinder = new ActiveRegionFinder(parser, refContentProvider, samContentProvider, rodContentProviders, true);

        activeRegionFinder.run();

        List<ActiveRegionMapData> activeRegionMapDataList = activeRegionFinder.getResultActiveRegions();

        HaplotypeCaller haplotypeCaller = new HaplotypeCaller(
                parser, refContentProvider, samContentProvider,
                rodContentProviders, activeRegionMapDataList, true);
        haplotypeCaller.run();

        List<VariantContext> resultList = haplotypeCaller.getResultVCFRecords();

        List<VcfRecord> gvcfs = getGvcf(refContigInfo);

        for (int i = 0; i < resultList.size(); i++) {
            VariantContext vc = resultList.get(i);
            VcfRecord vr = gvcfs.get(i);
            if (vc.getAlleles().size() == 2 && new String(vc.getAlleles().get(1).getBases()).equals("<NON_REF>")) {
                assertTrue(vr.alt().equals("<NON_REF>"));
                assertEquals(vc.getStart(), vr.position());
            } else {
                assertEquals(new String(vc.getReference().getBases()), vr.ref());
            }
        }
    }

    public void testHaplotypeCallerWithDbsnp() {
        RefContigInfo refContigInfo = RefContigInfo.apply(getClass().getResource("/human_g1k_v37.dict").getFile());
        SamHeaderInfo headerInfo = SamHeaderInfo.sortedHeader(refContigInfo, null);
        headerInfo.addReadGroupInfo(ReadGroupInfo.apply("SRR504516", "sample1"));

        SAMSequenceDictionary samSequenceDictionary = SAMSequenceDictTransfer.transfer(refContigInfo);
        GenomeLocParser parser = new GenomeLocParser(samSequenceDictionary);

        List<SAMRecord> recaledReads = getRecaledReads(headerInfo);
        List<GATKSAMRecord> realignedGATKRecords = new ArrayList<>();
        for (SAMRecord record : recaledReads) {
            realignedGATKRecords.add(new GATKSAMRecord(record));
        }
        SamContentProvider samContentProvider = new SamContentProvider(realignedGATKRecords, SAMHeaderTransfer.transfer(headerInfo));

        RefContentProvider refContentProvider = getRefContentProvider(samSequenceDictionary);

        VcfRecordPartition vcfRecordPartition = loadVcfPartition("/head_dbsnp.vcf", RODNames.DBSNP, refContigInfo);
        RODContentProvider rodContentProvider = new RODContentProvider(RODNames.DBSNP, vcfRecordPartition, parser);
        java.util.List<RODContentProvider> rodContentProviders = new java.util.ArrayList<>();
        rodContentProviders.add(rodContentProvider);

        ActiveRegionFinder activeRegionFinder = new ActiveRegionFinder(parser, refContentProvider, samContentProvider, rodContentProviders, true);

        activeRegionFinder.run();

        List<ActiveRegionMapData> activeRegionMapDataList = activeRegionFinder.getResultActiveRegions();

        HaplotypeCaller haplotypeCaller = new HaplotypeCaller(
                parser, refContentProvider, samContentProvider,
                rodContentProviders, activeRegionMapDataList, true);
        haplotypeCaller.run();

        List<VariantContext> resultList = haplotypeCaller.getResultVCFRecords();

        List<VcfRecord> gvcfs = getGvcf(refContigInfo);

        for (int i = 0; i < resultList.size(); i++) {
            VariantContext vc = resultList.get(i);
            VcfRecord vr = gvcfs.get(i);
            if (vc.getAlleles().size() == 2 && new String(vc.getAlleles().get(1).getBases()).equals("<NON_REF>")) {
                assertTrue(vr.alt().equals("<NON_REF>"));
                assertEquals(vc.getStart(), vr.position());
            } else {
                assertEquals(new String(vc.getReference().getBases()), vr.ref());
            }
        }
    }*/
}
