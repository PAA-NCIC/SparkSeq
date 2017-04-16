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

import com.github.lindenb.jbwa.jni.BwaIndex;
import com.github.lindenb.jbwa.jni.BwaMem;
import com.github.lindenb.jbwa.jni.ShortRead;
import org.apache.commons.lang3.StringUtils;
import org.ncic.bioinfo.sparkseq.data.basic.BasicSamRecord;
import org.ncic.bioinfo.sparkseq.data.basic.FastqPairRecord;
import org.ncic.bioinfo.sparkseq.data.common.ReadGroupInfo;
import org.ncic.bioinfo.sparkseq.data.common.RefContigInfo;
import org.ncic.bioinfo.sparkseq.exceptions.PipelineException;
import org.ncic.bioinfo.sparkseq.transfer.FastqRecord2ShortReadTransfer;
import scala.collection.JavaConversions;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * 该类用单例模式控制了每一个进程只加载一份index
 * <p>
 * Author: wbc
 */
public class BwaMemAdapter {

    private static volatile BwaIndex bwaIndexInstance = null;

    private static BwaIndex getBwaIndexInstance(String bwaJNILibPath, String referencePath) throws IOException {
        if (bwaIndexInstance == null) {
            synchronized (BwaMemAdapter.class) {
                if (bwaIndexInstance == null) {
                    System.load(bwaJNILibPath);
                    bwaIndexInstance = new BwaIndex(new File(referencePath));
                }
            }
        }
        return bwaIndexInstance;
    }


    public static List<BasicSamRecord> pairAlign(String bwaJNILibPath,
                                                 String referencePath,
                                                 ReadGroupInfo readGroupInfo,
                                                 RefContigInfo refContigInfo,
                                                 List<FastqPairRecord> fastqPairRecords) {

        List<BasicSamRecord> results = new ArrayList<>();
        BwaIndex bwaIndex = null;
        try {
            bwaIndex = getBwaIndexInstance(bwaJNILibPath, referencePath);
        } catch (IOException e) {
            e.printStackTrace();
            throw new PipelineException("Error when load index in JNI bwa");
        }

        int chunkSize = 1000;
        List<ShortRead> reads1 = new ArrayList<>(chunkSize);
        List<ShortRead> reads2 = new ArrayList<>(chunkSize);

        for(FastqPairRecord pairRecord : fastqPairRecords) {
            if(reads1.size() == chunkSize) {
                alignAndAddIntoResult(bwaIndex, results, readGroupInfo, refContigInfo, reads1, reads2);
                reads1.clear();
                reads2.clear();
            }

            reads1.add(FastqRecord2ShortReadTransfer.transferRead1(pairRecord));
            reads2.add(FastqRecord2ShortReadTransfer.transferRead2(pairRecord));
        }
        if(reads1.size() > 0) {
            alignAndAddIntoResult(bwaIndex, results, readGroupInfo, refContigInfo, reads1, reads2);
        }
        return results;
    }

    private static void alignAndAddIntoResult(BwaIndex bwaIndex, List<BasicSamRecord> results,
                                              ReadGroupInfo readGroupInfo, RefContigInfo refContigInfo,
                                              List<ShortRead> reads1, List<ShortRead> reads2){
        try {
            BwaMem mem = new BwaMem(bwaIndex);
            String[] rawResults = mem.align(reads1, reads2);
            String rgInfo = "RG:Z:" + readGroupInfo.id();
            String rgInfoWithTab = "\t" + rgInfo;
            for (String rawResult : rawResults) {
                // 这个库有时候会一次出现多条，中间以\n间隔
                String[] rawReads = StringUtils.split(rawResult, '\n');
                for (String rawRead : rawReads) {
                    String finalReadString = null;
                    if (rawRead.endsWith("\t")) {
                        finalReadString = rawRead + rgInfo;
                    } else {
                        finalReadString = rawRead + rgInfoWithTab;
                    }
                    results.add(BasicSamRecord.apply(finalReadString, refContigInfo, false));
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
            throw new PipelineException("Error when call JNI bwa");
        }
    }

}
