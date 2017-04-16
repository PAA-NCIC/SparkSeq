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

import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.FilterUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLoc;
import org.ncic.bioinfo.sparkseq.exceptions.GATKException;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Author: wbc
 */
public class RegionSamTraverser {

    private final ArrayList<GATKSAMRecord> samList;
    private int startIdx = 0;
    private int endIdx = 0;
    private int lastLocStart = 0;
    private int lastLocStop = 0;
    private final int readCount;

    public RegionSamTraverser(SamContentProvider samContentProvider) {
        this(samContentProvider, new FilterUtils());
    }

    public RegionSamTraverser(SamContentProvider samContentProvider, FilterUtils filterUtils) {

        samList = new ArrayList<>();
        samContentProvider.getGatksamRecords().forEach(record -> {
            if (filterUtils.filter(record)) {
                samList.add(record);
            }
        });

        readCount = samList.size();
        rewind();
    }


    /**
     * 将遍历器调整到初始状态，从头开始遍历
     */
    public void rewind() {
        startIdx = 0;
        endIdx = 0;
        lastLocStart = 0;
        lastLocStop = 0;
    }

    /**
     * 强烈建议每次访问的loc必须是顺序的，这样能减少查询的次数（目前版本是线性的）
     * 选出的region是startIdx和endIdx的左闭右开区间
     *
     * @param loc
     * @return
     */
    public List<GATKSAMRecord> getOverlappedReads(GenomeLoc loc) {
        if(loc.getStart() < lastLocStart || loc.getStop() < lastLocStop) {
            throw new GATKException("Locus to traverse is not in sequence");
        }
        lastLocStart = loc.getStart();
        lastLocStop = loc.getStop();
        while (startIdx < readCount && samList.get(startIdx).getAlignmentEnd() < lastLocStart) {
            startIdx++;
        }
        while (endIdx < readCount && samList.get(endIdx).getAlignmentStart() <= lastLocStop) {
            endIdx++;
        }
        if(startIdx < readCount) {
            List<GATKSAMRecord> resultRecords = new ArrayList<>(endIdx - startIdx);
            for(int i = startIdx; i < endIdx; i ++) {
                GATKSAMRecord record = samList.get(i);
                if(record.getAlignmentEnd() >= lastLocStart)    // 需要再次判断，因为有更短的read可能会没有overlap
                resultRecords.add(samList.get(i));
            }
            return resultRecords;
        } else {
            return Collections.EMPTY_LIST;
        }
    }
}
