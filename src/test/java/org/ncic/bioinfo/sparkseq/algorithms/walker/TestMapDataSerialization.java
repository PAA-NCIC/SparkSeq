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

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.esotericsoftware.kryo.serializers.JavaSerializer;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import org.ncic.bioinfo.sparkseq.algorithms.adapters.HaplotypeCallerAdapter;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.GATKSAMRecord;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.SamContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RODContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RODNames;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.header.StandardWGSVCFHeader;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLocParser;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.ActiveRegionFinder;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.ActiveRegionMapData;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.HaplotypeCaller;
import org.ncic.bioinfo.sparkseq.data.common.ReadGroupInfo;
import org.ncic.bioinfo.sparkseq.data.common.RefContigInfo;
import org.ncic.bioinfo.sparkseq.data.common.SamHeaderInfo;
import org.ncic.bioinfo.sparkseq.data.common.VcfHeaderInfo;
import org.ncic.bioinfo.sparkseq.data.partition.VcfRecordPartition;
import org.ncic.bioinfo.sparkseq.fileio.NormalFileLoader;
import org.ncic.bioinfo.sparkseq.transfer.*;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * Author: wbc
 */
public class TestMapDataSerialization extends AbstractTestCase {

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

        SAMRecord2BasicTransfer samTransfer = new SAMRecord2BasicTransfer();
        VCFHeader header = StandardWGSVCFHeader.getHeader();
        VC2VcfRecordTransfer vcfTransfer = new VC2VcfRecordTransfer(header, refContigInfo);

        VcfHeaderInfo vcfHeaderInfo = vcfRecordPartition.vcfHeader();
        HaplotypeCallerAdapter.StaticData staticData = HaplotypeCallerAdapter.getStaticDataInstance(refContigInfo, false, headerInfo, vcfHeaderInfo);
        HaplotypeCaller haplotypeCaller = staticData.haplotypeCaller;
        Basic2SAMRecordTransfer basic2SAMRecordTransfer = staticData.basic2SAMRecordTransfer;
        VCFHeader vcfFileHeader = staticData.vcfFileHeader;
        VCFCodec codec = staticData.codec;

        for(ActiveRegionMapData activeRegionMapData : activeRegionMapDataList) {
            SerializableActiveRegionMapData serial = new SerializableActiveRegionMapData(activeRegionMapData, samTransfer, vcfTransfer);
            byte[] data = serializationObject(serial);
            SerializableActiveRegionMapData deser = deserializationObject(data, SerializableActiveRegionMapData.class);
            ActiveRegionMapData back = deser.toActiveRegionMapData(
                    parser, basic2SAMRecordTransfer, vcfFileHeader, codec);

            assertEquals(activeRegionMapData.activeRegion.getReads().size(), back.activeRegion.getReads().size());
            haplotypeCaller.map(back);
        }

    }

    private <T extends Serializable> byte[] serializationObject(T obj) {
        Kryo kryo = new Kryo();
        kryo.setReferences(false);
        kryo.register(obj.getClass(), new JavaSerializer());

        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        Output output = new Output(baos);
        kryo.writeClassAndObject(output, obj);
        output.flush();
        output.close();

        byte[] b = baos.toByteArray();
        try {
            baos.flush();
            baos.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        return b;
    }

    @SuppressWarnings("unchecked")
    private <T extends Serializable> T deserializationObject(byte[] obj,
                                                             Class<T> clazz) {
        Kryo kryo = new Kryo();
        kryo.setReferences(false);
        kryo.register(clazz, new JavaSerializer());

        ByteArrayInputStream bais = new ByteArrayInputStream(
                obj);
        Input input = new Input(bais);
        return (T) kryo.readClassAndObject(input);
    }
}
