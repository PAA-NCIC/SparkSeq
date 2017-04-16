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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMLineParser;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import junit.framework.TestCase;
import org.apache.commons.lang3.StringUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLoc;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefContentProvider;
import org.ncic.bioinfo.sparkseq.data.basic.VcfRecord;
import org.ncic.bioinfo.sparkseq.data.common.RefContigInfo;
import org.ncic.bioinfo.sparkseq.data.common.SamHeaderInfo;
import org.ncic.bioinfo.sparkseq.data.common.VcfHeaderInfo;
import org.ncic.bioinfo.sparkseq.data.partition.FastaPartition;
import org.ncic.bioinfo.sparkseq.data.partition.VcfRecordPartition;
import org.ncic.bioinfo.sparkseq.fileio.NormalFileLoader;
import org.ncic.bioinfo.sparkseq.transfer.SAMHeaderTransfer;
import scala.collection.JavaConversions;
import scala.collection.immutable.List;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 * Author: wbc
 */
public abstract class AbstractTestCase extends TestCase {

    protected static String readRefContentValue(String filePath) {
        StringBuilder builder = new StringBuilder();
        String refContentFilePath = TestRealignerTargetCreator.class.getResource(filePath).getFile();
        try (BufferedReader reader = new BufferedReader(new FileReader(new File(refContentFilePath)))) {
            String line = reader.readLine();
            while (line != null) {
                builder.append(line);
                line = reader.readLine();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return builder.toString();
    }

    protected static String readRefContent() {
        return readRefContentValue("/refheadContent");
    }

    protected static String readLargeRefContent() {
        return readRefContentValue("/large/refheadContent");
    }

    protected static RefContentProvider getRefContentProvider(SAMSequenceDictionary samSequenceDictionary) {
        FastaPartition fastaPartition = new FastaPartition(1, 0, "1", readRefContent(),
                1, 100500, 1, 100000, 1, 97000);
        return new RefContentProvider(samSequenceDictionary, fastaPartition);
    }

    protected static RefContentProvider getLargeRefContentProvider(SAMSequenceDictionary samSequenceDictionary) {
        FastaPartition fastaPartition = new FastaPartition(1, 0, "1", readLargeRefContent(),
                99000, 10001000, 100000, 10000000, 103000, 9997000);
        return new RefContentProvider(samSequenceDictionary, fastaPartition);
    }

    protected static VcfRecordPartition loadVcfPartition(String filePath, String name, RefContigInfo refContigInfo) {
        String vcfFilePath = TestRealignerTargetCreator.class.getResource(filePath).getFile();
        VcfHeaderInfo vcfHeaderInfo = NormalFileLoader.loadVcfHeader(vcfFilePath);
        List<VcfRecord> vcfRecordList = NormalFileLoader.loadVcf(vcfFilePath, refContigInfo);
        return new VcfRecordPartition(1, name, 0, vcfHeaderInfo, vcfRecordList);
    }

    protected static java.util.List<GenomeLoc> getRealignerTargetIntervalValue(String filePath, RefContigInfo refContigInfo) {
        String intervalFilePath = TestRealignerTargetCreator.class.getResource(filePath).getFile();
        java.util.List<GenomeLoc> result = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(new File(intervalFilePath)))) {
            String line = reader.readLine();
            while (line != null) {
                if (line.length() > 0) {
                    String[] split1 = line.split(":");
                    String contigName = split1[0];
                    int contigId = refContigInfo.getId(contigName);
                    String[] split2 = split1[1].split("-");
                    int start = Integer.parseInt(split2[0]);
                    int end = (split2.length == 1) ? start : Integer.parseInt(split2[1]);
                    result.add(new GenomeLoc(contigName, contigId, start, end));
                }
                line = reader.readLine();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        return result;
    }

    protected static java.util.List<GenomeLoc> getRealignerTargetInterval(RefContigInfo refContigInfo) {
        return getRealignerTargetIntervalValue("/target_interval.list", refContigInfo);
    }

    protected static java.util.List<GenomeLoc> getLargeRealignerTargetInterval(RefContigInfo refContigInfo) {
        return getRealignerTargetIntervalValue("/large/target_interval.list", refContigInfo);
    }

    protected static java.util.List<SAMRecord> getOriginReads(SamHeaderInfo headerInfo) {
        return getReads("/intervaled.sam", headerInfo);
    }

    protected static java.util.List<SAMRecord> getRealignedReads(SamHeaderInfo headerInfo) {
        return getReads("/realigned_reads.sam", headerInfo);
    }

    protected static java.util.List<SAMRecord> getRecaledReads(SamHeaderInfo headerInfo) {
        return getReads("/recal_reads.sam", headerInfo);
    }

    protected static java.util.List<SAMRecord> getLargeOriginReads(SamHeaderInfo headerInfo) {
        return getReads("/large/target.sam", headerInfo);
    }

    protected static java.util.List<SAMRecord> getLargeRealignedReads(SamHeaderInfo headerInfo) {
        return getReads("/large/realigned_reads.sam", headerInfo);
    }

    protected static java.util.List<SAMRecord> getLargeRecaledReads(SamHeaderInfo headerInfo) {
        return getReads("/large/recal_reads.sam", headerInfo);
    }

    protected static java.util.List<SAMRecord> getReads(String filePath, SamHeaderInfo headerInfo) {

        String realPath = TestRealignerTargetCreator.class.getResource(filePath).getFile();

        SAMFileHeader header = SAMHeaderTransfer.transfer(headerInfo);
        SAMLineParser parser = new SAMLineParser(header);
        java.util.List<SAMRecord> result = new java.util.ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(new File(realPath)))) {
            String line = reader.readLine();
            while (line != null) {
                if (line.length() > 0 && !line.startsWith("@")) {
                    result.add(parser.parseLine(line));
                }
                line = reader.readLine();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        return result;
    }

    /**
     * Map的key是readName + flag
     *
     * @param recordList
     * @return
     */
    protected static Map<String, SAMRecord> readList2Map(java.util.List<SAMRecord> recordList) {
        Map<String, SAMRecord> resultMap = new HashMap<>();
        recordList.forEach(record -> resultMap.put(record.getReadName() + record.getFlags(), record));
        return resultMap;
    }

    private static java.util.List<String> getRecalTableLinesValue(String filePath) {
        java.util.List<String> res = new ArrayList<>();
        String realignedSamPath = TestRealignerTargetCreator.class.getResource(filePath).getFile();
        try (BufferedReader reader = new BufferedReader(new FileReader(new File(realignedSamPath)))) {
            String line = reader.readLine();
            while (line != null) {
                res.add(line);
                line = reader.readLine();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        return res;
    }

    protected static java.util.List<String> getRecalTableLines() {
        return getRecalTableLinesValue("/recal_data.table");
    }

    protected static java.util.List<String> getLargeRecalTableLines() {
        return getRecalTableLinesValue("/large/recal_data.table");
    }

    protected static java.util.List<VcfRecord> getGvcfValue(String filePath, RefContigInfo refContigInfo) {
        String realignedSamPath = TestRealignerTargetCreator.class.getResource(filePath).getFile();
        List<VcfRecord> vcfRecordList = NormalFileLoader.loadVcf(realignedSamPath, refContigInfo);
        return JavaConversions.asJavaList(vcfRecordList);
    }

    protected static java.util.List<VcfRecord> getGvcf(RefContigInfo refContigInfo) {
        return getGvcfValue("/gvcf.vcf", refContigInfo);
    }

    protected static java.util.List<VcfRecord> getLargeGvcf(RefContigInfo refContigInfo) {
        return getGvcfValue("/large/gvcf.vcf", refContigInfo);
    }

    protected static java.util.List<VcfRecord> getVcf(RefContigInfo refContigInfo) {
        return getGvcfValue("/raw.snps.indels.vcf", refContigInfo);
    }

    protected static java.util.List<VcfRecord> getLargeVcf(RefContigInfo refContigInfo) {
        return getGvcfValue("/large/raw.snps.indels.vcf", refContigInfo);
    }

    protected static java.util.List<GenomeLoc> getIntervals(RefContigInfo refContigInfo, String filePath) {
        String intervalPath = TestRealignerTargetCreator.class.getResource(filePath).getFile();
        java.util.List<GenomeLoc> intervals = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(new File(intervalPath)))) {
            String line = reader.readLine();
            while (line != null) {
                if(line.length() > 0 && !line.startsWith("@")) {
                    String[] split = StringUtils.split(line, '\t');
                    String contigName = split[0];
                    int contigId = refContigInfo.getId(contigName);
                    int start = Integer.parseInt(split[1]);
                    int stop = Integer.parseInt(split[2]);
                    intervals.add(new GenomeLoc(contigName, contigId, start, stop));
                }
                line = reader.readLine();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        return intervals;
    }
}
