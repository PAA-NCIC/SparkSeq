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
package org.ncic.bioinfo.sparkseq.algorithms.utils.transformers;

import org.ncic.bioinfo.sparkseq.algorithms.engine.Walker;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.GATKSAMRecord;
import org.ncic.bioinfo.sparkseq.algorithms.utils.reports.GATKReport;
import org.ncic.bioinfo.sparkseq.algorithms.walker.printreads.BaseRecalibration;

import java.io.File;
import java.util.Map;

/**
 * Author: wbc
 */
public class BQSRReadTransformer extends ReadTransformer {
    private boolean enabled;
    private BaseRecalibration bqsr = null;

    private GATKReport recalTable = null;

    public BQSRReadTransformer(GATKReport recalTable) {
        this.recalTable = recalTable;
    }

    @Override
    public OrderingConstraint getOrderingConstraint() {
        return OrderingConstraint.MUST_BE_FIRST;
    }

    @Override
    public void initializeSub(Map<String, Object> args, final Walker walker) {
        int quantizationLevels = (Integer)args.get("quantizationLevels");
        boolean disableIndelQuals = (Boolean)args.get("disableIndelQuals");
        boolean emitOriginalQuals = (Boolean)args.get("emitOriginalQuals");
        int PRESERVE_QSCORES_LESS_THAN = (Integer)args.get("PRESERVE_QSCORES_LESS_THAN");
        double globalQScorePrior = (Double)args.get("globalQScorePrior");
        // TODO -- See important note below about applying BQSR to a reduced BAM file:
        // If it is important to make sure that BQSR is not applied (as opposed to having the covariates computed) against a reduced bam file,
        // we need to figure out how to make this work.  The problem is that the ReadTransformers are initialized before the ReadDataSource
        // inside the GenomeAnalysisEngine, so we generate a NPE when trying to retrieve the SAMFileHeaders.  Ultimately, I don't think this is
        // a necessary check anyways since we disallow running BaseRecalibrator on reduced bams (so we can't generate the recal tables to use here).
        // Although we could add this check to the apply() method below, it's kind of ugly and inefficient.
        // The call here would be: RecalUtils.checkForInvalidRecalBams(engine.getSAMFileHeaders(), engine.getArguments().ALLOW_BQSR_ON_REDUCED_BAMS);
        this.bqsr = new BaseRecalibration(recalTable,
                quantizationLevels, disableIndelQuals, PRESERVE_QSCORES_LESS_THAN, emitOriginalQuals, globalQScorePrior);
    }

    @Override
    public boolean enabled() {
        return enabled;
    }

    /**
     * initialize a new BQSRReadTransformer that applies BQSR on the fly to incoming reads.
     */
    @Override
    public GATKSAMRecord apply(GATKSAMRecord read) {
        bqsr.recalibrateRead(read);
        return read;
    }
}
