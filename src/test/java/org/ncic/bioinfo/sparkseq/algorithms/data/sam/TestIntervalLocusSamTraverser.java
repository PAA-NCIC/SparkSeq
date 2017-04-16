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

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.DuplicateReadFilter;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.FailsVendorQualityCheckFilter;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.FilterUtils;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.NotPrimaryAlignmentFilter;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.UnmappedReadFilter;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RODContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RODNames;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLoc;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLocParser;
import org.ncic.bioinfo.sparkseq.algorithms.walker.AbstractTestCase;
import org.ncic.bioinfo.sparkseq.algorithms.walker.TestRealignerTargetCreator;
import org.ncic.bioinfo.sparkseq.data.common.ReadGroupInfo;
import org.ncic.bioinfo.sparkseq.data.common.RefContigInfo;
import org.ncic.bioinfo.sparkseq.data.common.SamHeaderInfo;
import org.ncic.bioinfo.sparkseq.data.partition.FastaPartition;
import org.ncic.bioinfo.sparkseq.transfer.SAMHeaderTransfer;
import org.ncic.bioinfo.sparkseq.transfer.SAMSequenceDictTransfer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Author: wbc
 */
public class TestIntervalLocusSamTraverser extends AbstractTestCase {

    public void testIntervalLocusSamTraverser() {
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

        List<GenomeLoc> intervals = getIntervals(refContigInfo, "/mutect/0_panel.intervals");

        FilterUtils filterUtils = getLocusWalkerFilterUtils();

        List<LocusInfo> locusInfos = getLocusInfo();

        GenomeLoc allLocus = refContentProvider.getLocus();
        IntervalLocusSamTraverser normalTraverser = new IntervalLocusSamTraverser(normalContentProvider, allLocus, intervals, filterUtils);
        normalTraverser.rewind();
        IntervalLocusSamTraverser tumorTraverser = new IntervalLocusSamTraverser(tumorContentProvider, allLocus, intervals, filterUtils);
        normalTraverser.rewind();

        Iterator<LocusInfo> locusInfoIterator = locusInfos.iterator();
        while (normalTraverser.hasNext()) {
            AlignmentContext tumor = tumorTraverser.next();
            AlignmentContext normal = normalTraverser.next();
            LocusInfo info = locusInfoIterator.next();

            assertEquals(info.count, tumor.size() + normal.size());
        }
    }

    protected static java.util.List<SAMRecord> getNormalReads(SamHeaderInfo headerInfo) {
        return getReads("/mutect/dedup_Blood.sam", headerInfo);
    }

    protected static java.util.List<SAMRecord> getTumorReads(SamHeaderInfo headerInfo) {
        return getReads("/mutect/dedup_Tumor.sam", headerInfo);
    }

    protected static RefContentProvider getMutectRefContentProvider(SAMSequenceDictionary samSequenceDictionary) {
        FastaPartition fastaPartition = new FastaPartition(1, 0, "1", readRefContentValue("/mutect/refContent"),
                1, 2501000, 1, 2500000, 1, 2497000);
        return new RefContentProvider(samSequenceDictionary, fastaPartition);
    }

    protected static FilterUtils getLocusWalkerFilterUtils() {
        FilterUtils filterUtils = new FilterUtils();
        filterUtils.addFilter(new UnmappedReadFilter());
        filterUtils.addFilter(new NotPrimaryAlignmentFilter());
        filterUtils.addFilter(new DuplicateReadFilter());
        filterUtils.addFilter(new FailsVendorQualityCheckFilter());
        return filterUtils;
    }

    private static List<LocusInfo> getLocusInfo() {
        List<LocusInfo> infos = new ArrayList<>();
        String locusTxt = TestRealignerTargetCreator.class.getResource("/mutect/LOCUS.txt").getFile();
        try (BufferedReader reader = new BufferedReader(new FileReader(new File(locusTxt)))) {
            String line = reader.readLine();
            while (line != null) {
                if (line.length() > 0) {
                    infos.add(new LocusInfo(line));
                }
                line = reader.readLine();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return infos;
    }

    static class LocusInfo {
        int position;
        int count;

        LocusInfo(String line) {
            String[] strs = line.split(" ");
            count = Integer.parseInt(strs[1]);
            position = Integer.parseInt(strs[0].split(":")[1]);
        }
    }
}
