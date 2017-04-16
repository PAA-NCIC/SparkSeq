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
package org.ncic.bioinfo.sparkseq.algorithms.walker.baserecalibrator;

import htsjdk.tribble.Feature;
import org.ncic.bioinfo.sparkseq.algorithms.utils.RecalUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.reports.GATKReportTable;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RodBinding;

import java.io.PrintStream;
import java.util.Collections;
import java.util.List;

/**
 * A collection of the arguments that are used for BQSR. Used to be common to both CovariateCounterWalker and TableRecalibrationWalker.
 * This set of arguments will also be passed to the constructor of every Covariate when it is instantiated.
 * <p>
 * Author: wbc
 */
public class RecalibrationArgumentCollection {

    /**
     * This algorithm treats every reference mismatch as an indication of error. However, real genetic variation is expected to mismatch the reference,
     * so it is critical that a database of known polymorphic sites is given to the tool in order to skip over those sites. This tool accepts any number of RodBindings (VCF, Bed, etc.)
     * for use as this database. For users wishing to exclude an interval list of known variation simply use -XL my.interval.list to skip over processing those sites.
     * Please note however that the statistics reported by the tool will not accurately reflected those sites skipped by the -XL argument.
     */
    public List<RodBinding<Feature>> knownSites = Collections.emptyList();

    /**
     * Note that the ReadGroup and QualityScore covariates are required and do not need to be specified.
     * Also, unless --no_standard_covs is specified, the Cycle and Context covariates are standard and are included by default.
     * Use the --list argument to see the available covariates.
     */
    public String[] COVARIATES = null;

    /*
 * The Cycle and Context covariates are standard and are included by default unless this argument is provided.
 * Note that the ReadGroup and QualityScore covariates are required and cannot be excluded.
 */
    public boolean DO_NOT_USE_STANDARD_COVARIATES = false;

    /**
     * This calculation is critically dependent on being able to skip over known polymorphic sites. Please be sure that you know what you are doing if you use this option.
     */
    public boolean RUN_WITHOUT_DBSNP = false;

    /**
     * BaseRecalibrator accepts a --solid_recal_mode <MODE> flag which governs how the recalibrator handles the
     * reads which have had the reference inserted because of color space inconsistencies.
     */
    public RecalUtils.SOLID_RECAL_MODE SOLID_RECAL_MODE = RecalUtils.SOLID_RECAL_MODE.SET_Q_ZERO;

    /**
     * BaseRecalibrator accepts a --solid_nocall_strategy <MODE> flag which governs how the recalibrator handles
     * no calls in the color space tag. Unfortunately because of the reference inserted bases mentioned above, reads with no calls in
     * their color space tag can not be recalibrated.
     */
    public RecalUtils.SOLID_NOCALL_STRATEGY SOLID_NOCALL_STRATEGY = RecalUtils.SOLID_NOCALL_STRATEGY.THROW_EXCEPTION;

    /**
     * The context covariate will use a context of this size to calculate its covariate value for base mismatches. Must be between 1 and 13 (inclusive). Note that higher values will increase runtime and required java heap size.
     */
    public int MISMATCHES_CONTEXT_SIZE = 2;

    /**
     * The context covariate will use a context of this size to calculate its covariate value for base insertions and deletions. Must be between 1 and 13 (inclusive). Note that higher values will increase runtime and required java heap size.
     */
    public int INDELS_CONTEXT_SIZE = 3;

    /**
     * The cycle covariate will generate an error if it encounters a cycle greater than this value.
     * This argument is ignored if the Cycle covariate is not used.
     */
    public int MAXIMUM_CYCLE_VALUE = 500;

    /**
     * A default base qualities to use as a prior (reported quality) in the mismatch covariate model. This value will replace all base qualities in the read for this default value. Negative value turns it off. [default is off]
     */
    public byte MISMATCHES_DEFAULT_QUALITY = -1;

    /**
     * A default base qualities to use as a prior (reported quality) in the insertion covariate model. This parameter is used for all reads without insertion quality scores for each base. [default is on]
     */
    public byte INSERTIONS_DEFAULT_QUALITY = 45;

    /**
     * A default base qualities to use as a prior (reported quality) in the mismatch covariate model. This value will replace all base qualities in the read for this default value. Negative value turns it off. [default is on]
     */
    public byte DELETIONS_DEFAULT_QUALITY = 45;

    /**
     * Reads with low quality bases on either tail (beginning or end) will not be considered in the context. This parameter defines the quality below which (inclusive) a tail is considered low quality
     */
    public byte LOW_QUAL_TAIL = 2;

    /**
     * BQSR generates a quantization table for quick quantization later by subsequent tools. BQSR does not quantize the base qualities, this is done by the engine with the -qq or -BQSR options.
     * This parameter tells BQSR the number of levels of quantization to use to build the quantization table.
     */
    public int QUANTIZING_LEVELS = 16;

    /**
     * The tag name for the binary tag covariate (if using it)
     */
    public String BINARY_TAG_NAME = null;

    /*
     * whether GATK report tables should have rows in sorted order, starting from leftmost column
     */
    public Boolean SORT_BY_ALL_COLUMNS = false;

    public String DEFAULT_PLATFORM = null;

    public String FORCE_PLATFORM = null;

    public String FORCE_READGROUP = null;

    public PrintStream RECAL_TABLE_UPDATE_LOG = null;

    /**
     * The repeat covariate will use a context of this size to calculate it's covariate value for base insertions and deletions
     */
    public int MAX_STR_UNIT_LENGTH = 8;

    public int MAX_REPEAT_LENGTH = 20;

    public GATKReportTable generateReportTable(final String covariateNames) {
        GATKReportTable argumentsTable;
        if (SORT_BY_ALL_COLUMNS) {
            argumentsTable = new GATKReportTable("Arguments", "Recalibration argument collection values used in this run", 2, GATKReportTable.TableSortingWay.SORT_BY_COLUMN);
        } else {
            argumentsTable = new GATKReportTable("Arguments", "Recalibration argument collection values used in this run", 2);
        }
        argumentsTable.addColumn("Argument");
        argumentsTable.addColumn(RecalUtils.ARGUMENT_VALUE_COLUMN_NAME);
        argumentsTable.addRowID("covariate", true);
        argumentsTable.set("covariate", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, covariateNames);
        argumentsTable.addRowID("no_standard_covs", true);
        argumentsTable.set("no_standard_covs", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, DO_NOT_USE_STANDARD_COVARIATES);
        argumentsTable.addRowID("run_without_dbsnp", true);
        argumentsTable.set("run_without_dbsnp", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, RUN_WITHOUT_DBSNP);
        argumentsTable.addRowID("solid_recal_mode", true);
        argumentsTable.set("solid_recal_mode", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, SOLID_RECAL_MODE);
        argumentsTable.addRowID("solid_nocall_strategy", true);
        argumentsTable.set("solid_nocall_strategy", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, SOLID_NOCALL_STRATEGY);
        argumentsTable.addRowID("mismatches_context_size", true);
        argumentsTable.set("mismatches_context_size", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, MISMATCHES_CONTEXT_SIZE);
        argumentsTable.addRowID("indels_context_size", true);
        argumentsTable.set("indels_context_size", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, INDELS_CONTEXT_SIZE);
        argumentsTable.addRowID("mismatches_default_quality", true);
        argumentsTable.set("mismatches_default_quality", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, MISMATCHES_DEFAULT_QUALITY);
        argumentsTable.addRowID("deletions_default_quality", true);
        argumentsTable.set("deletions_default_quality", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, DELETIONS_DEFAULT_QUALITY);
        argumentsTable.addRowID("insertions_default_quality", true);
        argumentsTable.set("insertions_default_quality", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, INSERTIONS_DEFAULT_QUALITY);
        argumentsTable.addRowID("maximum_cycle_value", true);
        argumentsTable.set("maximum_cycle_value", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, MAXIMUM_CYCLE_VALUE);
        argumentsTable.addRowID("low_quality_tail", true);
        argumentsTable.set("low_quality_tail", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, LOW_QUAL_TAIL);
        argumentsTable.addRowID("default_platform", true);
        argumentsTable.set("default_platform", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, DEFAULT_PLATFORM);
        argumentsTable.addRowID("force_platform", true);
        argumentsTable.set("force_platform", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, FORCE_PLATFORM);
        argumentsTable.addRowID("quantizing_levels", true);
        argumentsTable.set("quantizing_levels", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, QUANTIZING_LEVELS);
        argumentsTable.addRowID("recalibration_report", true);
        argumentsTable.set("recalibration_report", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, "null");
        argumentsTable.addRowID("binary_tag_name", true);
        argumentsTable.set("binary_tag_name", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, BINARY_TAG_NAME == null ? "null" : BINARY_TAG_NAME);
        return argumentsTable;
    }
}
