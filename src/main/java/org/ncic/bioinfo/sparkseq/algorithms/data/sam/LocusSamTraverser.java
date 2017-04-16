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
import org.ncic.bioinfo.sparkseq.algorithms.utils.ReadUtils;
import org.ncic.bioinfo.sparkseq.exceptions.GATKException;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLoc;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Author: wbc
 */
public class LocusSamTraverser {

    protected final List<GATKSAMRecord> samList;
    protected final String contigName;
    protected final int contigId;
    protected final int startCoordinate;
    protected final int endCoordinate;
    protected final int readCount;

    // 遍历过程中的变量
    protected int curCoordinate = 0;
    protected LinkedList<PileupTraker> bufferList = new LinkedList<>();
    // 当前指向的samList中的下标，以及当前代表下标指向read的tracker
    protected int curIndex = 0;
    protected PileupTraker curTracker = null;

    public LocusSamTraverser(SamContentProvider samContentProvider, GenomeLoc traverseInterval) {
        this(samContentProvider, traverseInterval, new FilterUtils());
    }

    public LocusSamTraverser(SamContentProvider samContentProvider, GenomeLoc traverseInterval,
                             FilterUtils filterUtils) {
        samList = samContentProvider.getGatksamRecords().stream()
                .filter(record -> filterUtils.filter(record))
                .collect(Collectors.toList());

        readCount = samList.size();
        contigName = traverseInterval.getContig();
        contigId = traverseInterval.getContigIndex();
        startCoordinate = traverseInterval.getStart();
        endCoordinate = traverseInterval.getStop();

        rewind();
    }

    public void rewind() {
        curCoordinate = startCoordinate;
        bufferList.clear();
        curIndex = 0;
        if (readCount > 0) {
            curTracker = new PileupTraker(samList.get(0));
        }
    }

    public boolean hasNext() {
        return curCoordinate <= endCoordinate;
    }

    public AlignmentContext next() {
        if (!hasNext()) {
            throw new GATKException("Out of bound when traverse sam content provider");
        }

        // 删除bufferList中已经超出的头部部分。通过判断尾部是否已经越过当前coordinate来删除
        // 因为read是不等长的，所以有可能更早的read没有超出，而晚一点的更短read已经超出，
        // 所以buffer中并不是所有的read都可以生成pileup。
        while (!bufferList.isEmpty()
                && !bufferList.getFirst().isBeforeEnd(curCoordinate)) {
            bufferList.removeFirst();
        }

        // 加入新的read
        // 需要判断：1.read尾在coordinate之前，直接跳过
        // 2.read头在coordinate之前，read尾在coordinate之后，加入buffer
        // 3.read头在coordinate之后，break
        while (curIndex < readCount) {
            if (!curTracker.isBeforeEnd(curCoordinate)) {
                curIndex++;
                if (curIndex < readCount) {
                    curTracker = new PileupTraker(samList.get(curIndex));
                }
                continue;
            }
            if (curTracker.isAfterStart(curCoordinate)) {
                // 如果遍历的位点截断了一些read，需要先做一个调整，让machine移动到当前位点
                int readStart = curTracker.stateMachine.getRead().getAlignmentStart();
                for (int i = 0; i < curCoordinate - readStart; i++) {
                    curTracker.stepForwardOnGenome();
                }
                bufferList.addLast(curTracker);
                curIndex++;
                if (curIndex < readCount) {
                    curTracker = new PileupTraker(samList.get(curIndex));
                }
            } else {
                break;
            }
        }

        //Get pileup to a arrayList
        List<PileupElement> pileupReads = new ArrayList<>();
        bufferList.forEach(tracker -> {
            if (tracker.isBeforeEnd(curCoordinate)) {
                tracker.stepForwardOnGenome();

                if (!ReadUtils.isBaseInsideAdaptor(tracker.stateMachine.getRead(), curCoordinate)) {
                    PileupElement element = tracker.stateMachine.makePileupElement();
                    pileupReads.add(element);
                }
            }
        });

        GenomeLoc locus = new GenomeLoc(contigName, contigId, curCoordinate, curCoordinate);
        ReadBackedPileup readBackedPileup = new ReadBackedPileupImpl(locus, pileupReads);

        curCoordinate++;
        return new AlignmentContext(locus, readBackedPileup);
    }

    class PileupTraker {
        private int start = 0;
        private int end = 0;
        private AlignmentStateMachine stateMachine;

        PileupTraker(GATKSAMRecord read) {
            this.stateMachine = new AlignmentStateMachine(read);
            this.start = read.getAlignmentStart();
            this.end = read.getAlignmentEnd();
        }

        void stepForwardOnGenome() {
            stateMachine.stepForwardOnGenome();
        }

        boolean isBeforeEnd(int coordinate) {
            return (coordinate <= end);
        }

        boolean isAfterStart(int coordinate) {
            return (coordinate >= start);
        }
    }

    /**
     * 设置新的coordinate，但是不能向后倒退
     *
     * @param newCoordinate
     */
    public void forwardToCoordinate(int newCoordinate) {
        if(this.curCoordinate > newCoordinate) {
            throw new GATKException("Can't go back when traverse in locus mode");
        }
        this.curCoordinate = newCoordinate;
    }

}
