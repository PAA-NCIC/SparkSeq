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
package org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.readlikelihood;

import org.ncic.bioinfo.sparkseq.algorithms.utils.haplotype.Haplotype;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.GATKSAMRecord;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.AssemblyResultSet;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.SampleList;

import java.util.List;
import java.util.Map;

/**
 * Author: wbc
 */
public interface ReadLikelihoodCalculationEngine {

    enum Implementation {
        /**
         * Classic full pair-hmm all haplotypes vs all reads.
         */
        PairHMM,

        /**
         * Graph-base likelihoods.
         */
        GraphBased,

        /**
         * Random likelihoods, used to establish a baseline benchmark for other meaningful implementations.
         */
        Random
    }


    /**
     * Calculates the likelihood of reads across many samples evaluated against haplotypes resulting from the
     * active region assembly process.
     *
     * @param assemblyResultSet the input assembly results.
     * @param samples the list of targeted samples.
     * @param perSampleReadList the input read sets stratified per sample.
     *
     * @throws NullPointerException if either parameter is {@code null}.
     *
     * @return never {@code null}, and with at least one entry for input sample (keys in {@code perSampleReadList}.
     *    The value maps can be potentially empty though.
     */
    public ReadLikelihoods<Haplotype> computeReadLikelihoods(AssemblyResultSet assemblyResultSet, SampleList samples,
                                                             Map<String, List<GATKSAMRecord>> perSampleReadList);

    public void close();
}
