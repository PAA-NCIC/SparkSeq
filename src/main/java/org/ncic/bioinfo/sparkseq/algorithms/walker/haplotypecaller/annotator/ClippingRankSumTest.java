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
package org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator;

import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.ncic.bioinfo.sparkseq.algorithms.utils.AlignmentUtils;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.GATKSAMRecord;

import java.util.*;

/**
 * Rank Sum Test for hard-clipped bases on REF vs. ALT reads
 *
 * <p>This variant-level annotation tests whether the data supporting the reference allele shows more or less base clipping (hard clips) than those supporting the alternate allele. The ideal result is a value close to zero, which indicates there is little to no difference.  A negative value indicates that the reads supporting the alternate allele have more hard-clipped bases than those supporting the reference allele. Conversely, a positive value indicates that the reads supporting the alternate allele have fewer hard-clipped bases than those supporting the reference allele. Finding a statistically significant difference either way suggests that the sequencing and/or mapping process may have been biased or affected by an artifact.</p>
 *
 * <h3>Statistical notes</h3>
 * <p>The value output for this annotation is the u-based z-approximation from the Mann-Whitney-Wilcoxon Rank Sum Test applied to base clips (number of hard-clipped bases on reads supporting REF vs. number of hard-clipped bases on reads supporting ALT). See the <a href="http://www.broadinstitute.org/gatk/guide/article?id=">method document on statistical tests</a> for a more detailed explanation of the ranksum test.</p>
 *
 * <h3>Caveat</h3>
 * <p>The clipping rank sum test can not be calculated for sites without a mixture of reads showing both the reference and alternate alleles.</p>
 *
 * Author: wbc
 */
public class ClippingRankSumTest extends RankSumTest {
    @Override
    public List<String> getKeyNames() { return Arrays.asList("ClippingRankSum"); }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("ClippingRankSum", 1, VCFHeaderLineType.Float, "Z-score From Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases")); }

    @Override
    protected Double getElementForRead(final GATKSAMRecord read, final int refLoc) {
        return (double) AlignmentUtils.getNumHardClippedBases(read);
    }
 }
