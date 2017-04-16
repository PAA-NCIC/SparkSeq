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
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.GATKSAMRecord;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.PileupElement;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.interfaces.StandardAnnotation;

import java.util.*;


/**
 * Rank Sum Test for mapping qualities of REF vs. ALT reads
 *
 * <p>This variant-level annotation compares the mapping qualities of the reads supporting the reference allele with those supporting the alternate allele. The ideal result is a value close to zero, which indicates there is little to no difference. A negative value indicates that the reads supporting the alternate allele have lower mapping quality scores than those supporting the reference allele. Conversely, a positive value indicates that the reads supporting the alternate allele have higher mapping quality scores than those supporting the reference allele.</p>
 * <p>This annotation can be used to evaluate confidence in a variant call and is a recommended covariate for variant recalibration (VQSR). Finding a statistically significant difference in quality either way suggests that the sequencing and/or mapping process may have been biased or affected by an artifact. In practice, we only filter out low negative values when evaluating variant quality because the idea is to filter out variants for which the quality of the data supporting the alternate allele is comparatively low. The reverse case, where it is the quality of data supporting the reference allele that is lower (resulting in positive ranksum scores), is not really informative for filtering variants.
 *
 * <h3>Statistical notes</h3>
 * <p>The value output for this annotation is the u-based z-approximation from the Mann-Whitney-Wilcoxon Rank Sum Test for mapping qualities (MAPQ of reads supporting REF vs. MAPQ of reads supporting ALT). See the <a href="http://www.broadinstitute.org/gatk/guide/article?id=4732">method document on statistical tests</a> for a more detailed explanation of the ranksum test.</p>
 *
 * <h3>Caveat</h3>
 * <p>The mapping quality rank sum test can not be calculated for sites without a mixture of reads showing both the reference and alternate alleles.</p>
 *
 * <h3>Related annotations</h3>
 * <ul>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_RMSMappingQuality.php">RMSMappingQuality</a></b> gives an estimation of the overal read mapping quality supporting a variant call.</li>
 * </ul>
 *
 */
public class MappingQualityRankSumTest extends RankSumTest implements StandardAnnotation {
    @Override
    public List<String> getKeyNames() { return Arrays.asList("MQRankSum"); }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("MQRankSum", 1, VCFHeaderLineType.Float, "Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities")); }

    @Override
    protected Double getElementForRead(final GATKSAMRecord read, final int refLoc) {
        return (double)read.getMappingQuality();
    }

    @Override
    protected Double getElementForPileupElement(final PileupElement p) {
        return (double)p.getRead().getMappingQuality();
    }
}