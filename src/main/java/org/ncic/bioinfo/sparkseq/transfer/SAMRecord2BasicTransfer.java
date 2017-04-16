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

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextTagCodec;
import org.ncic.bioinfo.sparkseq.data.basic.BasicSamRecord;
import scala.collection.JavaConversions;

import java.util.ArrayList;
import java.util.List;

/**
 * Author: wbc
 */
public class SAMRecord2BasicTransfer {
    private final static int FAKE_CONTIG_ID = 255;
    private static final String FIELD_SEPARATOR = "\t";

    private final TextTagCodec tagCodec = new TextTagCodec();

    /**
     * Write the record.
     *
     * @param alignment SAMRecord.
     */
    public BasicSamRecord transfer(final SAMRecord alignment) {
        String readName = alignment.getReadName();
        int flags = alignment.getFlags();
        int contigId = alignment.getReferenceIndex();
        if (contigId == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
            contigId = FAKE_CONTIG_ID;
        }
        String contigName = alignment.getReferenceName();
        int position = alignment.getAlignmentStart();
        int mapQ = alignment.getMappingQuality();
        String cigar = alignment.getCigarString();
        int mateContigId = alignment.getMateReferenceIndex();
        if (mateContigId == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
            mateContigId = FAKE_CONTIG_ID;
        }
        String mateContigName = alignment.getMateReferenceName();
        int matePosition = alignment.getMateAlignmentStart();
        int inferredSize = alignment.getInferredInsertSize();
        byte[] sequence = alignment.getReadString().getBytes();
        byte[] quality = alignment.getBaseQualityString().getBytes();
        List<SAMRecord.SAMTagAndValue> attributes = alignment.getAttributes();
        List<String> encodedTags = new ArrayList<>(attributes.size());
        for (SAMRecord.SAMTagAndValue attribute : attributes) {
            encodedTags.add(tagCodec.encode(attribute.tag, attribute.value));
        }

        // 在这里统一不压缩，压缩操作是由scala代码中的执行引擎统一调配的。
        return BasicSamRecord.apply(false, readName, flags, contigId, contigName, position, mapQ,
                cigar, mateContigId, mateContigName, matePosition, inferredSize,
                sequence, quality, encodedTags);
    }
}
