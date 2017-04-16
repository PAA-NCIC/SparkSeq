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
import org.apache.log4j.Logger;
import org.ncic.bioinfo.sparkseq.algorithms.utils.AlignmentContextUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLocParser;
import org.ncic.bioinfo.sparkseq.algorithms.utils.haplotype.Haplotype;
import org.ncic.bioinfo.sparkseq.algorithms.utils.pairhmm.PairHMMIndelErrorModel;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefMetaDataTracker;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.ReferenceContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.AlignmentContext;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.GeneralPloidyGenotypeLikelihoods;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.PerReadAlleleLikelihoodMap;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.argcollection.UnifiedArgumentCollection;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Author: wbc
 */
public class GeneralPloidyIndelGenotypeLikelihoodsCalculationModel extends GeneralPloidyGenotypeLikelihoodsCalculationModel {
    private static final int MAX_NUM_ALLELES_TO_GENOTYPE = 4;

    private PairHMMIndelErrorModel pairModel;
 /*
    private static ThreadLocal<HashMap<PileupElement, LinkedHashMap<Allele, Double>>> indelLikelihoodMap =
            new ThreadLocal<HashMap<PileupElement, LinkedHashMap<Allele, Double>>>() {
                protected synchronized HashMap<PileupElement, LinkedHashMap<Allele, Double>> initialValue() {
                    return new HashMap<PileupElement, LinkedHashMap<Allele, Double>>();
                }
            };
   */

    private LinkedHashMap<Allele, Haplotype> haplotypeMap;

     /*
    static {
        indelLikelihoodMap.set(new HashMap<PileupElement, LinkedHashMap<Allele, Double>>());
    }
       */

    protected GeneralPloidyIndelGenotypeLikelihoodsCalculationModel(final UnifiedArgumentCollection UAC, final Logger logger) {
        super(UAC, logger);


        pairModel = new PairHMMIndelErrorModel(UAC.INDEL_GAP_OPEN_PENALTY, UAC.INDEL_GAP_CONTINUATION_PENALTY,
                UAC.OUTPUT_DEBUG_INDEL_INFO, UAC.pairHMM);
        haplotypeMap = new LinkedHashMap<Allele, Haplotype>();
    }


    protected GeneralPloidyGenotypeLikelihoods getPoolGenotypeLikelihoodObject(final List<Allele> alleles,
                                                                               final double[] logLikelihoods,
                                                                               final int ploidy,
                                                                               final HashMap<String, ErrorModel> perLaneErrorModels,
                                                                               final boolean useBQAedPileup,
                                                                               final ReferenceContext ref,
                                                                               final boolean ignoreLaneInformation,
                                                                               final PerReadAlleleLikelihoodMap perReadAlleleLikelihoodMap){
        return new GeneralPloidyIndelGenotypeLikelihoods(alleles, logLikelihoods, ploidy,perLaneErrorModels,ignoreLaneInformation, pairModel, haplotypeMap, ref, perReadAlleleLikelihoodMap);
    }

    protected List<Allele> getInitialAllelesToUse(final RefMetaDataTracker tracker,
                                                  final ReferenceContext ref,
                                                  final Map<String, AlignmentContext> contexts,
                                                  final AlignmentContextUtils.ReadOrientation contextType,
                                                  final GenomeLocParser locParser,
                                                  final List<Allele> allAllelesToUse){


        List<Allele> alleles = IndelGenotypeLikelihoodsCalculationModel.getInitialAlleleList(tracker, ref, contexts, contextType, UAC,true);

        if (alleles.size() > MAX_NUM_ALLELES_TO_GENOTYPE)
            alleles = alleles.subList(0,MAX_NUM_ALLELES_TO_GENOTYPE);
        if (contextType == AlignmentContextUtils.ReadOrientation.COMPLETE) {
            haplotypeMap.clear();
        }
        IndelGenotypeLikelihoodsCalculationModel.getHaplotypeMapFromAlleles(alleles, ref, ref.getLocus(), haplotypeMap);

        // sanity check: if haplotype map couldn't be created, clear allele list
        if (haplotypeMap.isEmpty())
            alleles.clear();
        return alleles;

    }

    protected List<Allele> getFinalAllelesToUse(final RefMetaDataTracker tracker,
                                                final ReferenceContext ref,
                                                final List<Allele> allAllelesToUse,
                                                final ArrayList<PoolGenotypeData> GLs) {

        // find the alternate allele(s) that we should be using
        final List<Allele> alleles = new ArrayList<Allele>();
        if ( allAllelesToUse != null )
            alleles.addAll(allAllelesToUse);
        else if (!GLs.isEmpty())
            alleles.addAll(GLs.get(0).alleles);
        return alleles;

    }

    protected int getEndLocation(final RefMetaDataTracker tracker,
                                 final ReferenceContext ref,
                                 final List<Allele> allelesToUse) {
        return ref.getLocus().getStart() + allelesToUse.get(0).length() - 1;
    }
}