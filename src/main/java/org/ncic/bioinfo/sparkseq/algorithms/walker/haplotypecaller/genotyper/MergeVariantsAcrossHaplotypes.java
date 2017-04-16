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
package org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.genotyper;

import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLoc;
import org.ncic.bioinfo.sparkseq.algorithms.utils.haplotype.Haplotype;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.readlikelihood.ReadLikelihoods;

import java.util.List;
import java.util.TreeSet;

/**
 * Author: wbc
 */
public class MergeVariantsAcrossHaplotypes {
    /**
     * Merge variants across the haplotypes, updating the haplotype event maps and startPos set as appropriate
     *
     * @param haplotypes a list of haplotypes whose events we want to merge
     * @param readLikelihoods map from sample name -> read likelihoods for each haplotype
     * @param startPosKeySet a set of starting positions of all events among the haplotypes
     * @param ref the reference bases
     * @param refLoc the span of the reference bases
     * @return true if anything was merged
     */
    public boolean merge( final List<Haplotype> haplotypes,
                          final ReadLikelihoods<Haplotype> readLikelihoods,
                          final TreeSet<Integer> startPosKeySet,
                          final byte[] ref,
                          final GenomeLoc refLoc ) {
        return false;
    }
}
