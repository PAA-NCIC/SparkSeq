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
package org.ncic.bioinfo.sparkseq.algorithms.engine;

import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLocParser;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefMetaDataTracker;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.SamContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RODContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.ActiveRegion;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.ActiveRegionMapData;

import java.util.List;

/**
 * Author: wbc
 */
public abstract class ActiveRegionWalker extends Walker {

    List<ActiveRegionMapData> activeRegionMapDataList = null;

    public ActiveRegionWalker(GenomeLocParser genomeLocParser,
                              RefContentProvider refContentProvider,
                              SamContentProvider samContentProvider,
                              List<RODContentProvider> rodContentProviderList,
                              List<ActiveRegionMapData> activeRegionMapDataList) {
        super(genomeLocParser, refContentProvider, samContentProvider, rodContentProviderList);
        this.activeRegionMapDataList = activeRegionMapDataList;
    }

    public void run() {

        for (ActiveRegionMapData activeRegionMapData : activeRegionMapDataList) {
            map(activeRegionMapData);
        }

        onTraversalDone();
    }

    protected abstract void map(ActiveRegionMapData activeRegionMapData);

    protected abstract void onTraversalDone();

}
