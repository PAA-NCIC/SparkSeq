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
package org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.model;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.log4j.Logger;
import org.ncic.bioinfo.sparkseq.algorithms.utils.AlignmentContextUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.BaseUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLocParser;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefMetaDataTracker;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.ReferenceContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.AlignmentContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.PileupElement;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.ReadBackedPileup;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.PerReadAlleleLikelihoodMap;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.argcollection.UnifiedArgumentCollection;
import org.ncic.bioinfo.sparkseq.exceptions.ReviewedGATKException;

import java.util.List;
import java.util.Map;

/**
 * Author: wbc
 */
public abstract class GenotypeLikelihoodsCalculationModel {

    public static final String DUMMY_LANE = "Lane1";
    public static final String DUMMY_SAMPLE_NAME = "DummySample1";

    public enum Model {
        SNP,
        INDEL,
        GENERALPLOIDYSNP,
        GENERALPLOIDYINDEL,
        BOTH;
    }

    protected final UnifiedArgumentCollection UAC;
    protected Logger logger;

    /**
     * Create a new object
     * @param logger        logger
     * @param UAC           unified arg collection
     */
    protected GenotypeLikelihoodsCalculationModel(UnifiedArgumentCollection UAC, Logger logger) {
        if ( logger == null || UAC == null ) throw new ReviewedGATKException("Bad arguments");
        this.UAC = UAC;
        this.logger = logger;
    }

    /**
     * Can be overridden by concrete subclasses
     *
     * @param tracker               rod data
     * @param ref                   reference context
     * @param contexts              stratified alignment contexts
     * @param contextType           stratified context type
     * @param allAllelesToUse the alternate allele to use, null if not set
     * @param useBAQedPileup        should we use the BAQed pileup or the raw one?
     * @param locParser             Genome Loc Parser
     * @return variant context where genotypes are no-called but with GLs
     */
    public abstract VariantContext getLikelihoods(final RefMetaDataTracker tracker,
                                                  final ReferenceContext ref,
                                                  final Map<String, AlignmentContext> contexts,
                                                  final AlignmentContextUtils.ReadOrientation contextType,
                                                  final List<Allele> allAllelesToUse,
                                                  final boolean useBAQedPileup,
                                                  final GenomeLocParser locParser,
                                                  final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap);


    protected int getFilteredDepth(ReadBackedPileup pileup) {
        int count = 0;
        for ( PileupElement p : pileup ) {
            if ( BaseUtils.isRegularBase( p.getBase() ) )
                count++;
        }

        return count;
    }

}