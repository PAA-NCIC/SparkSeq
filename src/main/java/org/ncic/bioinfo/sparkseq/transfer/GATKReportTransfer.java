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

import org.ncic.bioinfo.sparkseq.algorithms.utils.reports.GATKReport;
import org.ncic.bioinfo.sparkseq.algorithms.utils.reports.GATKReportTable;
import org.ncic.bioinfo.sparkseq.algorithms.utils.reports.GATKReportVersion;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Author: wbc
 */
public class GATKReportTransfer {

    public static GATKReport lines2Report(List<String> lines) {
        Iterator<String> lineIter = lines.iterator();
        GATKReport report = new GATKReport();

        // version line
        String reportHeaderLine = lineIter.next();
        GATKReportVersion version = GATKReportVersion.fromHeader(reportHeaderLine);
        report.setVersion(version);

        // tables
        int nTables = Integer.parseInt(reportHeaderLine.split(":")[2]);
        // Read each table according ot the number of tables
        for (int i = 0; i < nTables; i++) {
            report.addTable(new GATKReportTable(lineIter, version));
        }

        return report;
    }

    public static List<String> report2Lines(GATKReport report) {
        List<String> lines = new ArrayList<>();
        // version line
        lines.add(GATKReport.GATKREPORT_HEADER_PREFIX + report.getVersion().toString()
                + GATKReport.SEPARATOR + report.getTables().size());

        for (GATKReportTable table : report.getTables()) {
            String[] tableLines = table.transIntoLines();
            for(String str : tableLines) {
                lines.add(str);
            }
            lines.add("");
        }
        return lines;
    }

}
