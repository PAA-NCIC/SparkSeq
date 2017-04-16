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
package org.ncic.bioinfo.sparkseq.algorithms.walker.mutect;

/**
 * Author: wbc
 */
public class MuTectArgumentCollection {
    public boolean NOOP = false;

    public boolean ENABLE_EXTENDED_OUTPUT = false;

    public boolean ENABLE_QSCORE_OUTPUT = false;

    public boolean ARTIFACT_DETECTION_MODE = false;

    public String TUMOR_SAMPLE_NAME = "Tumor";

    public String BAM_TUMOR_SAMPLE_NAME = null;

    public String NORMAL_SAMPLE_NAME = "Normal";

    public boolean FORCE_OUTPUT = false;

    public boolean FORCE_ALLELES = false;

    public boolean ONLY_PASSING_CALLS = false;

    public float INITIAL_TUMOR_LOD_THRESHOLD = 4.0f;

    public float TUMOR_LOD_THRESHOLD = 6.3f;

    public float FRACTION_CONTAMINATION = 0.02f;

    public float MINIMUM_MUTATION_CELL_FRACTION = 0.00f;

    public float NORMAL_LOD_THRESHOLD = 2.2f;

    public float NORMAL_ARTIFACT_LOD_THRESHOLD = 1.0f;

    public float STRAND_ARTIFACT_LOD_THRESHOLD = 2.0f;

    public float STRAND_ARTIFACT_POWER_THRESHOLD = 0.9f;

    public float NORMAL_DBSNP_LOD_THRESHOLD = 5.5f;

    public float MINIMUM_NORMAL_ALLELE_FRACTION = 0.00f;

    public float TUMOR_F_PRETEST = 0.005f;

    public int MIN_QSCORE = 5;

    public int GAP_EVENTS_THRESHOLD = 3;

    public float HEAVILY_CLIPPED_READ_FRACTION = 0.30f;

    public float FRACTION_MAPQ0_THRESHOLD = 0.5f;

    public double PIR_MEDIAN_THRESHOLD = 10;

    public double PIR_MAD_THRESHOLD = 3;

    public int REQUIRED_MAXIMUM_ALT_ALLELE_MAPPING_QUALITY_SCORE = 20;

    /**
     * Parameters for ALT ALLELE IN NORMAL filter
     **/
    public int MAX_ALT_ALLELES_IN_NORMAL_COUNT = 2;

    public int MAX_ALT_ALLELES_IN_NORMAL_QSCORE_SUM = 20;

    public double MAX_ALT_ALLELE_IN_NORMAL_FRACTION = 0.03;

    public int POWER_CONSTANT_QSCORE = 30;

    public double POWER_CONSTANT_AF = 0.3f;

}
