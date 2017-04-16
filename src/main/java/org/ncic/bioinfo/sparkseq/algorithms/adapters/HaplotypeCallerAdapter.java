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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.SamContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RODContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RODNames;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.VCFHeaderLineIterable;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.header.StandardWGSVCFHeader;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLocParser;
import org.ncic.bioinfo.sparkseq.algorithms.walker.SerializableActiveRegionMapData;
import org.ncic.bioinfo.sparkseq.algorithms.walker.genotypegvcfs.GenotypeGVCFs;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.ActiveRegionFinder;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.ActiveRegionMapData;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.HaplotypeCaller;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.SampleList;
import org.ncic.bioinfo.sparkseq.data.basic.VcfRecord;
import org.ncic.bioinfo.sparkseq.data.common.ReadGroupInfo;
import org.ncic.bioinfo.sparkseq.data.common.RefContigInfo;
import org.ncic.bioinfo.sparkseq.data.common.SamHeaderInfo;
import org.ncic.bioinfo.sparkseq.data.common.VcfHeaderInfo;
import org.ncic.bioinfo.sparkseq.data.partition.FastaPartition;
import org.ncic.bioinfo.sparkseq.data.partition.SamRecordPartition;
import org.ncic.bioinfo.sparkseq.data.partition.VcfRecordPartition;
import org.ncic.bioinfo.sparkseq.transfer.*;
import scala.collection.JavaConversions;

import java.util.ArrayList;
import java.util.List;

/**
 * Author: wbc
 */
public class HaplotypeCallerAdapter {

    public static List<SerializableActiveRegionMapData> getActiveRegions(RefContigInfo refContigInfo,
                                                                         SamRecordPartition samRecordPartition,
                                                                         FastaPartition refPartition,
                                                                         List<VcfRecordPartition> rodPartitions,
                                                                         boolean useGVCF) {

        SAMSequenceDictionary samSequenceDictionary = SAMSequenceDictTransfer.transfer(refContigInfo);
        GenomeLocParser parser = new GenomeLocParser(samSequenceDictionary);

        SamContentProvider samContentProvider = new SamContentProvider(samRecordPartition);
        RefContentProvider refContentProvider = new RefContentProvider(samSequenceDictionary, refPartition);

        List<RODContentProvider> rodContentProviders = new java.util.ArrayList<>();
        rodPartitions.forEach(
                rodPartition -> rodContentProviders.add(
                        new RODContentProvider(rodPartition.key(), rodPartition, parser))
        );

        // 先找出active region的信息
        ActiveRegionFinder activeRegionFinder = new ActiveRegionFinder(
                parser, refContentProvider, samContentProvider, rodContentProviders, useGVCF);
        activeRegionFinder.run();
        List<ActiveRegionMapData> activeRegionMapDatas = activeRegionFinder.getResultActiveRegions();

        SAMRecord2BasicTransfer samTransfer = new SAMRecord2BasicTransfer();
        VCFHeader header = StandardWGSVCFHeader.getHeader();
        VC2VcfRecordTransfer vcfTransfer = new VC2VcfRecordTransfer(header, refContigInfo);

        List<SerializableActiveRegionMapData> serializableActiveRegionMapDataList = new ArrayList<>(activeRegionMapDatas.size());
        for (ActiveRegionMapData mapData : activeRegionMapDatas) {
            serializableActiveRegionMapDataList.add(new SerializableActiveRegionMapData(mapData, samTransfer, vcfTransfer));
        }

        return serializableActiveRegionMapDataList;
    }

    public static class StaticData {
        public HaplotypeCaller haplotypeCaller = null;
        public GenomeLocParser genomeLocParser = null;
        public Basic2SAMRecordTransfer basic2SAMRecordTransfer = null;
        public VCFHeader vcfFileHeader = null;
        public VCFCodec codec = null;
    }

    public static StaticData getStaticDataInstance(RefContigInfo refContigInfo,
                                                   boolean useGVCF,
                                                   SamHeaderInfo samHeaderInfo,
                                                   VcfHeaderInfo vcfHeaderInfo) {
        StaticData data = new StaticData();
        SAMSequenceDictionary samSequenceDictionary = SAMSequenceDictTransfer.transfer(refContigInfo);
        data.genomeLocParser = new GenomeLocParser(samSequenceDictionary);
        samHeaderInfo.addReadGroupInfo(ReadGroupInfo.apply("rg1", "sample1"));
        SAMFileHeader header = SAMHeaderTransfer.transfer(samHeaderInfo);
        List<SAMReadGroupRecord> readGroupInfos = header.getReadGroups();
        List<String> samples = new ArrayList<>();
        for (SAMReadGroupRecord readGroup : readGroupInfos) {
            samples.add(readGroup.getSample());
        }
        data.haplotypeCaller = new HaplotypeCaller(data.genomeLocParser, samples, useGVCF);
        data.basic2SAMRecordTransfer = new Basic2SAMRecordTransfer(header);

        VCFCodec codec = new VCFCodec();
        VCFHeaderLineIterable headerLineIterable = new VCFHeaderLineIterable(vcfHeaderInfo);
        data.vcfFileHeader = (VCFHeader) codec.readActualHeader(headerLineIterable);
        data.codec = codec;
        return data;
    }

    public static List<VcfRecord> callVariants(RefContigInfo refContigInfo,
                                               scala.collection.immutable.List<SerializableActiveRegionMapData> serializableActiveRegionMapDataList,
                                               boolean useGVCF,
                                               SamHeaderInfo samHeaderInfo,
                                               VcfHeaderInfo vcfHeaderInfo) {
        return callVariants(refContigInfo, CollectionConverter.asJavaList(serializableActiveRegionMapDataList),
                useGVCF, samHeaderInfo, vcfHeaderInfo);
    }

    public static List<VcfRecord> callVariants(RefContigInfo refContigInfo,
                                               List<SerializableActiveRegionMapData> serializableActiveRegionMapDataList,
                                               boolean useGVCF,
                                               SamHeaderInfo samHeaderInfo,
                                               VcfHeaderInfo vcfHeaderInfo) {
        StaticData staticData = getStaticDataInstance(refContigInfo, false, samHeaderInfo, vcfHeaderInfo);
        GenomeLocParser parser = staticData.genomeLocParser;
        HaplotypeCaller haplotypeCaller = staticData.haplotypeCaller;
        Basic2SAMRecordTransfer basic2SAMRecordTransfer = staticData.basic2SAMRecordTransfer;
        VCFHeader vcfFileHeader = staticData.vcfFileHeader;
        VCFCodec codec = staticData.codec;

        List<VariantContext> variantContexts = new ArrayList<>();

        for(SerializableActiveRegionMapData mapData: serializableActiveRegionMapDataList) {
            ActiveRegionMapData activeRegionMapData = mapData.toActiveRegionMapData(
                    parser, basic2SAMRecordTransfer, vcfFileHeader, codec);
            variantContexts.addAll(haplotypeCaller.map(activeRegionMapData));
        }

        // 如果使用了gvcf，则需要加一个genotypeGVCFs
        List<VariantContext> finalResult = variantContexts;

        VCFHeader header = StandardWGSVCFHeader.getHeader();
        VC2VcfRecordTransfer transfer = new VC2VcfRecordTransfer(header, refContigInfo);
        List<VcfRecord> vcfRecords = new ArrayList<>(finalResult.size());
        finalResult.forEach(vc -> vcfRecords.add(transfer.transfer(vc)));
        return vcfRecords;
    }

}
