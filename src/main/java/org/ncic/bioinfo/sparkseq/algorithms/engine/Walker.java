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

import htsjdk.samtools.SAMReadGroupRecord;
import org.apache.log4j.Logger;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLocParser;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.SamContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RODContentProvider;

import java.util.List;

/**
 * Author: wbc
 */
public abstract class Walker {

    final protected static Logger logger = Logger.getLogger(Walker.class);

    protected GenomeLocParser genomeLocParser;
    protected RefContentProvider refContentProvider;
    protected SamContentProvider samContentProvider;
    protected List<RODContentProvider> rodContentProviderList;

    public Walker(GenomeLocParser genomeLocParser,
                  RefContentProvider refContentProvider,
                  SamContentProvider samContentProvider,
                  List<RODContentProvider> rodContentProviderList) {
        this.genomeLocParser = genomeLocParser;
        this.refContentProvider = refContentProvider;
        this.samContentProvider = samContentProvider;
        this.rodContentProviderList = rodContentProviderList;
    }

    protected abstract void initialize();

    public abstract void run();

    public String getSampleName() {
        List<SAMReadGroupRecord> readGroupInfos = samContentProvider.getSamFileHeader().getReadGroups();
        return (readGroupInfos.size() > 0) ? readGroupInfos.get(0).getSample() : "Sample1";
    }

}
