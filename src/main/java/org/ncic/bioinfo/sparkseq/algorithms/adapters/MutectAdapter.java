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

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.SamContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RODContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.header.StandardWGSVCFHeader;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLoc;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLocParser;
import org.ncic.bioinfo.sparkseq.algorithms.walker.mutect.Mutect;
import org.ncic.bioinfo.sparkseq.data.basic.VcfRecord;
import org.ncic.bioinfo.sparkseq.data.common.Locus;
import org.ncic.bioinfo.sparkseq.data.common.RefContigInfo;
import org.ncic.bioinfo.sparkseq.data.partition.FastaPartition;
import org.ncic.bioinfo.sparkseq.data.partition.SamRecordPartition;
import org.ncic.bioinfo.sparkseq.data.partition.VcfRecordPartition;
import org.ncic.bioinfo.sparkseq.transfer.SAMSequenceDictTransfer;
import org.ncic.bioinfo.sparkseq.transfer.VC2VcfRecordTransfer;
import scala.collection.JavaConversions;

import java.util.ArrayList;
import java.util.List;

/**
 * Author: wbc
 */
public class MutectAdapter {
    public static List<VcfRecord> callVariants(RefContigInfo refContigInfo,
                                               SamRecordPartition tumorSamRecordPartition,
                                               SamRecordPartition normalSamRecordPartition,
                                               FastaPartition refPartition,
                                               List<VcfRecordPartition> rodPartitions,
                                               List<Locus> intervals) {
        SAMSequenceDictionary samSequenceDictionary = SAMSequenceDictTransfer.transfer(refContigInfo);
        GenomeLocParser parser = new GenomeLocParser(samSequenceDictionary);

        SamContentProvider tumorSamContentProvider = new SamContentProvider(tumorSamRecordPartition);
        SamContentProvider normalSamContentProvider = new SamContentProvider(normalSamRecordPartition);
        RefContentProvider refContentProvider = new RefContentProvider(samSequenceDictionary, refPartition);

        List<RODContentProvider> rodContentProviders = new java.util.ArrayList<>();
        rodPartitions.forEach(
                rodPartition -> rodContentProviders.add(
                        new RODContentProvider(rodPartition.key(), rodPartition, parser))
        );

        List<GenomeLoc> intervalLocus = new ArrayList<>();
        GenomeLoc traverseLocus = refContentProvider.getLocus();
        intervals.forEach(
                locus -> {
                    GenomeLoc interval = new GenomeLoc(locus.contigName(), locus.contigId(), locus.start(), locus.stop());
                    if (interval.overlapsP(traverseLocus)) {
                        intervalLocus.add(interval);
                    }
                }
        );

        Mutect mutect = new Mutect(parser, refContentProvider,
                tumorSamContentProvider, normalSamContentProvider, rodContentProviders, intervalLocus);

        List<VariantContext> finalResult = mutect.getResultVCFRecords();

        VCFHeader header = StandardWGSVCFHeader.getHeader();
        VC2VcfRecordTransfer transfer = new VC2VcfRecordTransfer(header, refContigInfo);
        List<VcfRecord> vcfRecords = new ArrayList<>(finalResult.size());
        finalResult.forEach(vc -> vcfRecords.add(transfer.transfer(vc)));
        return vcfRecords;
    }
}
