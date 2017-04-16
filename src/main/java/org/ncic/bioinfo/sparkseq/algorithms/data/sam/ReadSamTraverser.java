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
import org.ncic.bioinfo.sparkseq.exceptions.GATKException;

import java.util.ArrayList;

/**
 * Author: wbc
 */
public class ReadSamTraverser {

    private final ArrayList<GATKSAMRecord> samList;
    private int index = 0;
    private final int readCount;

    public ReadSamTraverser(SamContentProvider samContentProvider) {
        this(samContentProvider, new FilterUtils());
    }

    public ReadSamTraverser(SamContentProvider samContentProvider, FilterUtils filterUtils) {

        samList = new ArrayList<>();
        samContentProvider.getGatksamRecords().forEach(record -> {
            if(filterUtils.filter(record)) {
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
        index = 0;
    }

    public boolean hasNext() {
        return index < readCount;
    }

    public GATKSAMRecord next() {
        if (!hasNext()) {
            throw new GATKException("Out of bound when traverse sam content provider");
        }
        return samList.get(index++);
    }
}
