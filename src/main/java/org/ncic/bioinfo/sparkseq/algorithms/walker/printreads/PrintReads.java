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
package org.ncic.bioinfo.sparkseq.algorithms.walker.printreads;

import htsjdk.samtools.SAMReadGroupRecord;
import org.ncic.bioinfo.sparkseq.algorithms.engine.ReadWalker;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLocParser;
import org.ncic.bioinfo.sparkseq.algorithms.utils.RandomGenerator;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefMetaDataTracker;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.ReferenceContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.GATKSAMRecord;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.SamContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.Filter;
import org.ncic.bioinfo.sparkseq.algorithms.utils.reports.GATKReport;
import org.ncic.bioinfo.sparkseq.algorithms.utils.transformers.BQSRReadTransformer;
import org.ncic.bioinfo.sparkseq.algorithms.utils.transformers.MisencodedBaseQualityReadTransformer;
import org.ncic.bioinfo.sparkseq.algorithms.utils.transformers.ReadTransformer;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RODContentProvider;
import org.ncic.bioinfo.sparkseq.exceptions.ReviewedGATKException;
import org.ncic.bioinfo.sparkseq.exceptions.UserException;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;

/**
 * Author: wbc
 */
public class PrintReads extends ReadWalker {

    private GATKReport recalTable = null;

    private String readGroup = null;

    /**
     * For example, --platform ILLUMINA or --platform 454.
     */
    private String platform = null;

    /**
     * Only prints the first n reads of the file
     */
    int nReadsToPrint = -1;

    /**
     * Only reads from samples listed in the provided file(s) will be included in the output.
     */
    public Set<File> sampleFile = new TreeSet<>();

    /**
     * Only reads from the sample(s) will be included in the output.
     */
    public Set<String> sampleNames = new TreeSet<>();

    /**
     * Erase all extra attributes in the read but keep the read group information
     */
    public boolean simplifyReads = false;

    public boolean NO_PG_TAG = false;

    private Set<String> readGroupsToKeep = Collections.emptySet();

    public static final String PROGRAM_RECORD_NAME = "GATK PrintReads";   // The name that will go in the @PG tag

    private Random random;

    private List<ReadTransformer> readTransformers = null;

    // result
    private List<GATKSAMRecord> result = new ArrayList<>();

    public PrintReads(GenomeLocParser genomeLocParser,
                      RefContentProvider refContentProvider,
                      SamContentProvider samContentProvider,
                      List<RODContentProvider> rodContentProviderList,
                      GATKReport recalTable) {
        super(genomeLocParser, refContentProvider, samContentProvider, rodContentProviderList);
        this.recalTable = recalTable;
    }

    @Override
    public void initialize() {
        List<ReadTransformer> activeTransformers = new ArrayList<>();
        BQSRReadTransformer bqsrReadTransformer = new BQSRReadTransformer(recalTable);
        Map<String, Object> bqsrReadTransformerArgs = buildBQSRTransformerArgsMap();
        bqsrReadTransformer.initialize(ReadTransformer.ApplicationTime.HANDLED_IN_WALKER, bqsrReadTransformerArgs, this);
        activeTransformers.add(bqsrReadTransformer);

        MisencodedBaseQualityReadTransformer misencodedBaseQualityReadTransformer =
                new MisencodedBaseQualityReadTransformer();
        Map<String, Object> baseQualityTransformerArgs = buildMisencodedBaseQualityReadTransformerArgsMap();
        misencodedBaseQualityReadTransformer.initialize(ReadTransformer.ApplicationTime.HANDLED_IN_WALKER, baseQualityTransformerArgs, this);
        activeTransformers.add(misencodedBaseQualityReadTransformer);

        // 给transformer赋值
        setReadTransformers(activeTransformers);

        random = RandomGenerator.getRandomGenerator();
    }

    protected void setReadTransformers(final List<ReadTransformer> readTransformers) {
        if (readTransformers == null)
            throw new ReviewedGATKException("read transformers cannot be null");

        // sort them in priority order
        Collections.sort(readTransformers, new ReadTransformer.ReadTransformerComparator());

        // make sure we don't have an invalid set of active read transformers
        checkActiveReadTransformers(readTransformers);

        this.readTransformers = readTransformers;
    }

    protected void checkActiveReadTransformers(final List<ReadTransformer> readTransformers) {
        if (readTransformers == null)
            throw new IllegalArgumentException("read transformers cannot be null");

        ReadTransformer sawMustBeFirst = null;
        ReadTransformer sawMustBeLast = null;

        for (final ReadTransformer r : readTransformers) {
            if (r.getOrderingConstraint() == ReadTransformer.OrderingConstraint.MUST_BE_FIRST) {
                if (sawMustBeFirst != null)
                    throw new UserException.IncompatibleReadFiltersException(sawMustBeFirst.toString(), r.toString());
                sawMustBeFirst = r;
            } else if (r.getOrderingConstraint() == ReadTransformer.OrderingConstraint.MUST_BE_LAST) {
                if (sawMustBeLast != null)
                    throw new UserException.IncompatibleReadFiltersException(sawMustBeLast.toString(), r.toString());
                sawMustBeLast = r;
            }
        }
    }

    private Map<String, Object> buildBQSRTransformerArgsMap() {
        Map<String, Object> args = new HashMap<>();
        args.put("quantizationLevels", 0);
        args.put("disableIndelQuals", false);
        args.put("emitOriginalQuals", false);
        args.put("PRESERVE_QSCORES_LESS_THAN", 6);
        args.put("globalQScorePrior", -1.0);
        return args;
    }

    private Map<String, Object> buildMisencodedBaseQualityReadTransformerArgsMap() {
        Map<String, Object> args = new HashMap<>();
        args.put("fixQuals", false);
        args.put("ALLOW_POTENTIALLY_MISENCODED_QUALS", false);
        return args;
    }

    @Override
    protected List<Filter> getFilter() {
        List<Filter> filters = new ArrayList<>();
        return filters;
    }

    public boolean filter(ReferenceContext ref, GATKSAMRecord read) {
        // check that the read belongs to an RG that we need to keep
        if (!readGroupsToKeep.isEmpty()) {
            final SAMReadGroupRecord readGroup = read.getReadGroup();
            if (!readGroupsToKeep.contains(readGroup.getReadGroupId()))
                return false;
        }

        // check if we've reached the output limit
        if (nReadsToPrint == 0) {
            return false;          // n == 0 means we've printed all we needed.
        } else if (nReadsToPrint > 0) {
            nReadsToPrint--;       // n > 0 means there are still reads to be printed.
        }

        return true;
    }

    @Override
    protected void map(final ReferenceContext ref,
                       final GATKSAMRecord originalRead,
                       final RefMetaDataTracker metaDataTracker) {
        if (!filter(ref, originalRead)) {
            return;
        }

        GATKSAMRecord workingRead = originalRead;

        for (final ReadTransformer transformer : readTransformers) {
            workingRead = transformer.apply(workingRead);
        }

        if (simplifyReads) workingRead = workingRead.simplify();
        result.add(workingRead);
    }

    @Override
    protected void onTraversalDone() {

    }

    public List<GATKSAMRecord> getResultRecords() {
        return result;
    }
}
