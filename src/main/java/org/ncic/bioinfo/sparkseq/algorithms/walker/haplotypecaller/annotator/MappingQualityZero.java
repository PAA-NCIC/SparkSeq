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

import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import htsjdk.variant.variantcontext.VariantContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefMetaDataTracker;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.ReferenceContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.AlignmentContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.GATKSAMRecord;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.PileupElement;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.ReadBackedPileup;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.PerReadAlleleLikelihoodMap;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.interfaces.ActiveRegionBasedAnnotation;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.interfaces.AnnotatorCompatible;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.interfaces.InfoFieldAnnotation;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.interfaces.StandardAnnotation;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * Count of all reads with MAPQ = 0 across all samples
 * <p>
 * <p>This anotation gives you the count of all reads that have MAPQ = 0 across all samples. The count of reads with MAPQ0 can be used for quality control; high counts typically indicate regions where it is difficult to make confident calls.</p>
 * <p>
 * <h3>Related annotations</h3>
 * <ul>
 * <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_MappingQualityZeroBySample.php">MappingQualityZeroBySample</a></b> gives the count of reads with MAPQ=0 for each individual sample.</li>
 * <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_LowMQ.php">LowMQ</a></b> gives the proportion of reads with low mapping quality (MAPQ below 10, including 0).</li>
 * </ul>
 */
public class MappingQualityZero extends InfoFieldAnnotation implements StandardAnnotation, ActiveRegionBasedAnnotation {

    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {
        if ((vc.isSNP() || !vc.isVariant()) && stratifiedContexts != null)
            return annotatePileup(ref, stratifiedContexts, vc);
        else if (stratifiedPerReadAlleleLikelihoodMap != null && vc.isVariant())
            return annotateWithLikelihoods(stratifiedPerReadAlleleLikelihoodMap, vc);
        else
            return null;
    }

    private Map<String, Object> annotatePileup(final ReferenceContext ref,
                                               final Map<String, AlignmentContext> stratifiedContexts,
                                               final VariantContext vc) {
        if (stratifiedContexts.size() == 0)
            return null;

        int mq0 = 0;
        for (Map.Entry<String, AlignmentContext> sample : stratifiedContexts.entrySet()) {
            final AlignmentContext context = sample.getValue();
            final ReadBackedPileup pileup = context.getBasePileup();
            for (PileupElement p : pileup) {
                if (p.getMappingQual() == 0)
                    mq0++;
            }
        }
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%d", mq0));
        return map;
    }

    private Map<String, Object> annotateWithLikelihoods(final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap,
                                                        final VariantContext vc) {
        if (stratifiedPerReadAlleleLikelihoodMap == null)
            return null;

        int mq0 = 0;
        for (PerReadAlleleLikelihoodMap likelihoodMap : stratifiedPerReadAlleleLikelihoodMap.values()) {
            for (GATKSAMRecord read : likelihoodMap.getLikelihoodReadMap().keySet()) {

                if (read.getMappingQuality() == 0)
                    mq0++;
            }
        }
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%d", mq0));
        return map;
    }


    public List<String> getKeyNames() {
        return Arrays.asList(VCFConstants.MAPPING_QUALITY_ZERO_KEY);
    }

    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(VCFStandardHeaderLines.getInfoLine(getKeyNames().get(0)));
    }
}