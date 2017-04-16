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

import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFormatException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordFactory;
import htsjdk.samtools.TagValueAndUnsignedArrayFlag;
import htsjdk.samtools.TextTagCodec;
import htsjdk.samtools.ValidationStringency;
import org.ncic.bioinfo.sparkseq.compress.BaseCompressTools;
import org.ncic.bioinfo.sparkseq.compress.QualityCompressTools;
import org.ncic.bioinfo.sparkseq.data.basic.BasicSamRecord;

import scala.collection.JavaConversions;
import scala.collection.convert.WrapAsJava;

import java.util.List;
import java.util.Map;

/**
 * Author: wbc
 */
public class Basic2SAMRecordTransfer {

    private final static int FAKE_CONTIG_ID = 255;

    private final SAMFileHeader mFileHeader;
    private final SAMRecordFactory samRecordFactory;
    private final ValidationStringency validationStringency;
    private final TextTagCodec tagCodec = new TextTagCodec();

    public Basic2SAMRecordTransfer(SAMFileHeader header) {
        this.mFileHeader = header;
        this.samRecordFactory = new DefaultSAMRecordFactory();
        this.validationStringency = ValidationStringency.DEFAULT_STRINGENCY;
    }

    public SAMRecord transfer(final BasicSamRecord basicSamRecord) {
        final SAMRecord samRecord =
                samRecordFactory.createSAMRecord(this.mFileHeader);
        samRecord.setValidationStringency(this.validationStringency);
        samRecord.setHeader(this.mFileHeader);
        samRecord.setReadName(basicSamRecord.readName());
        samRecord.setFlags(basicSamRecord.flag());
        if (!basicSamRecord.contigName().equals("*")) {
            samRecord.setReferenceName(basicSamRecord.contigName());
        }
        samRecord.setAlignmentStart(basicSamRecord.position());
        samRecord.setMappingQuality(basicSamRecord.mapQ());

        samRecord.setCigarString(basicSamRecord.cigar());
        if (basicSamRecord.contigId() != FAKE_CONTIG_ID) {
            samRecord.setMateReferenceName(basicSamRecord.mateContigName());
            samRecord.setMateAlignmentStart(basicSamRecord.matePosition());
            samRecord.setInferredInsertSize(basicSamRecord.infferdSize());
        }

        byte[] quality = null;
        if (basicSamRecord.quality().length == 1 && basicSamRecord.quality()[0] == '*') {
            samRecord.setBaseQualities(SAMRecord.NULL_QUALS);
        } else {
            quality = basicSamRecord.quality();
            if (basicSamRecord.compressFlag()) {
                quality = QualityCompressTools.deCompressQual(quality);
            }
        }

        byte[] sequence = null;
        if (basicSamRecord.sequence().length == 1 && basicSamRecord.sequence()[0] == '*') {
            samRecord.setReadBases(SAMRecord.NULL_SEQUENCE);
        } else {
            sequence = basicSamRecord.sequence();
            if (basicSamRecord.compressFlag()) {
                sequence = BaseCompressTools.decompressBase(sequence, quality);
            }
            // 必须在执行之后解压之后再赋值，因为解压可能会导致quality变化。
            samRecord.setBaseQualityString(new String(quality));
            samRecord.setReadString(new String(sequence));
        }

        List<String> tmp = CollectionConverter.asJavaList(basicSamRecord.attributeList());
        for (String attribute : tmp) {
            parseTag(samRecord, attribute);
        }

        return samRecord;
    }

    private void parseTag(final SAMRecord samRecord, final String tag) {
        Map.Entry<String, Object> entry = null;
        try {
            entry = tagCodec.decode(tag);
        } catch (SAMFormatException e) {
            e.printStackTrace();
        }
        if (entry != null) {
            if (entry.getValue() instanceof TagValueAndUnsignedArrayFlag) {
                final TagValueAndUnsignedArrayFlag valueAndFlag =
                        (TagValueAndUnsignedArrayFlag) entry.getValue();
                if (valueAndFlag.isUnsignedArray) {
                    samRecord.setUnsignedArrayAttribute(entry.getKey(),
                            valueAndFlag.value);
                } else {
                    samRecord.setAttribute(entry.getKey(), valueAndFlag.value);
                }
            } else {
                samRecord.setAttribute(entry.getKey(), entry.getValue());
            }
        }
    }

}
