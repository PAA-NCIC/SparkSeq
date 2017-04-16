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

import org.ncic.bioinfo.sparkseq.algorithms.utils.reports.GATKReport;
import org.ncic.bioinfo.sparkseq.algorithms.walker.printreads.RecalibrationReport;
import org.ncic.bioinfo.sparkseq.exceptions.GATKException;
import org.ncic.bioinfo.sparkseq.exceptions.PipelineException;
import org.ncic.bioinfo.sparkseq.transfer.GATKReportTransfer;
import scala.collection.JavaConversions;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.RecursiveTask;
import java.util.concurrent.TimeUnit;

/**
 * Author: wbc
 */
public class BQSRTableGather {

    public static List<String> gatherBQSRTables(List<List<String>> reportLinesGroup) {

        if (reportLinesGroup.size() == 0) {
            throw new PipelineException("No bqsr tables to merge");
        }

        final RecalibrationReport finalReport = new RecalibrationReport(
                GATKReportTransfer.lines2Report(reportLinesGroup.get(0)));

        int count = reportLinesGroup.size();
        for (int i = 1; i < count; i++) {
            List<String> reportLines = reportLinesGroup.get(i);
            GATKReport gatkReport = GATKReportTransfer.lines2Report(reportLines);
            if (gatkReport.getTable("RecalTable0").getNumRows() != 0) {
                finalReport.combine(new RecalibrationReport(gatkReport));
            }
        }

        finalReport.calculateQuantizedQualities();
        GATKReport gatkReport = finalReport.createGATKReport();
        return GATKReportTransfer.report2Lines(gatkReport);
    }

    public static List<String> gatherBQSRTablesInParallel(final List<List<String>> rawReportLinesGroup, int threadCount, int tableCount) {
        ExecutorService pool = Executors.newFixedThreadPool(threadCount);

        int tableSize = (rawReportLinesGroup.size() < tableCount) ? rawReportLinesGroup.size() : tableCount;
        List<List<String>> reportLinesGroup = rawReportLinesGroup.subList(0, tableSize);

        RecalibrationReport[] reports = new RecalibrationReport[reportLinesGroup.size()];

        for (int i = 0; i < reportLinesGroup.size(); i++) {
            final int reportId = i;
            pool.execute(new Runnable() {
                @Override
                public void run() {
                    List<String> lines = reportLinesGroup.get(reportId);
                    GATKReport gatkReport = GATKReportTransfer.lines2Report(lines);
                    reports[reportId] = new RecalibrationReport(gatkReport);
                }
            });
        }

        pool.shutdown();

        try {
            while (!pool.awaitTermination(1, TimeUnit.SECONDS)) ;
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        RecalibrationReport finalReport = reports[0];
        for (int i = 1; i < reports.length; i++) {
            if (reports[i].getRecalibrationTables().getTable(0).getDimensions()[0] != 0) {
                finalReport.combine(reports[i]);
            }
        }
        finalReport.calculateQuantizedQualities();
        GATKReport gatkReport = finalReport.createGATKReport();
        return GATKReportTransfer.report2Lines(gatkReport);

        /*ForkJoinPool forkJoinPool = new ForkJoinPool();
        Future<RecalibrationReport> futureResult = forkJoinPool.submit(new ForkJoinMergeTask(reports, 0, reports.length - 1));

        try {
            RecalibrationReport finalReport = futureResult.get();
            finalReport.calculateQuantizedQualities();
            GATKReport gatkReport = finalReport.createGATKReport();
            return GATKReportTransfer.report2Lines(gatkReport);
        } catch (Exception e) {
            throw new GATKException("Error when combine BQSR report:" + e.getMessage());
        }*/
    }

    /**
     * 暂时弃用，因为结果会不一致
     */
    private static class ForkJoinMergeTask extends RecursiveTask<RecalibrationReport> {
        private int start;
        private int stop;
        RecalibrationReport[] reports;

        public ForkJoinMergeTask(RecalibrationReport[] reports, int start, int end) {
            this.reports = reports;
            this.start = start;
            this.stop = end;
        }

        public RecalibrationReport compute() {
            if (start == stop) {
                return reports[start];
            } else if (start == stop - 1) {
                reports[start].combine(reports[stop]);
                return reports[start];
            } else {
                int middle = (start + stop) / 2;
                ForkJoinMergeTask leftTask = new ForkJoinMergeTask(reports, start, middle);
                ForkJoinMergeTask rightTask = new ForkJoinMergeTask(reports, middle + 1, stop);
                leftTask.fork();
                rightTask.fork();
                RecalibrationReport leftRes = leftTask.join();
                RecalibrationReport rightRes = rightTask.join();
                leftRes.combine(rightRes);
                return leftRes;
            }
        }
    }
}
