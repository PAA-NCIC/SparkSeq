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
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.log4j.Logger;
import org.ncic.bioinfo.sparkseq.algorithms.utils.AlignmentContextUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.BaseUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GATKVariantContextUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLoc;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLocParser;
import org.ncic.bioinfo.sparkseq.algorithms.utils.haplotype.Haplotype;
import org.ncic.bioinfo.sparkseq.algorithms.utils.pairhmm.PairHMMIndelErrorModel;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefMetaDataTracker;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.ReferenceContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.AlignmentContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.PileupElement;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.ReadBackedPileup;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.ConsensusAlleleCounter;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.GenotypingOutputMode;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.PerReadAlleleLikelihoodMap;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.argcollection.UnifiedArgumentCollection;

import java.util.ArrayList;
import java.util.EnumSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Author: wbc
 */
public class IndelGenotypeLikelihoodsCalculationModel extends GenotypeLikelihoodsCalculationModel {
    private static final int HAPLOTYPE_SIZE = 80;

    private boolean DEBUG = false;
    private boolean ignoreSNPAllelesWhenGenotypingIndels = false;
    private PairHMMIndelErrorModel pairModel;


    private LinkedHashMap<Allele, Haplotype> haplotypeMap;

    private List<Allele> alleleList = new ArrayList<Allele>();


    protected IndelGenotypeLikelihoodsCalculationModel(final UnifiedArgumentCollection UAC,
                                                       final Logger logger) {
        super(UAC, logger);
        pairModel = new PairHMMIndelErrorModel(UAC.INDEL_GAP_OPEN_PENALTY, UAC.INDEL_GAP_CONTINUATION_PENALTY,
                UAC.OUTPUT_DEBUG_INDEL_INFO, UAC.pairHMM);
        DEBUG = UAC.OUTPUT_DEBUG_INDEL_INFO;
        haplotypeMap = new LinkedHashMap<Allele, Haplotype>();
        ignoreSNPAllelesWhenGenotypingIndels = UAC.IGNORE_SNP_ALLELES;
    }

    protected static List<Allele> computeConsensusAlleles(final ReferenceContext ref,
                                                          final Map<String, AlignmentContext> contexts,
                                                          final AlignmentContextUtils.ReadOrientation contextType,
                                                          final UnifiedArgumentCollection UAC) {
        ConsensusAlleleCounter counter = new ConsensusAlleleCounter(true, UAC.MIN_INDEL_COUNT_FOR_GENOTYPING, UAC.MIN_INDEL_FRACTION_PER_SAMPLE);
        return counter.computeConsensusAlleles(ref, contexts, contextType);
    }

    private final static EnumSet<VariantContext.Type> allowableTypes = EnumSet.of(VariantContext.Type.INDEL, VariantContext.Type.MIXED);


    public VariantContext getLikelihoods(final RefMetaDataTracker tracker,
                                         final ReferenceContext ref,
                                         final Map<String, AlignmentContext> contexts,
                                         final AlignmentContextUtils.ReadOrientation contextType,
                                         final List<Allele> allAllelesToUse,
                                         final boolean useBAQedPileup,
                                         final GenomeLocParser locParser,
                                         final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {

        GenomeLoc loc = ref.getLocus();
//        if (!ref.getLocus().equals(lastSiteVisited)) {
        if (contextType == AlignmentContextUtils.ReadOrientation.COMPLETE) {
            // starting a new site: clear allele list
            haplotypeMap.clear();
            perReadAlleleLikelihoodMap.clear(); // clean mapping sample-> per read, per allele likelihoods
            alleleList = getInitialAlleleList(tracker, ref, contexts, contextType, UAC, ignoreSNPAllelesWhenGenotypingIndels);
            if (alleleList.isEmpty())
                return null;
        }

        getHaplotypeMapFromAlleles(alleleList, ref, loc, haplotypeMap); // will update haplotypeMap adding elements
        if (haplotypeMap == null || haplotypeMap.isEmpty())
            return null;

        // start making the VariantContext
        // For all non-snp VC types, VC end location is just startLocation + length of ref allele including padding base.
        final int endLoc = loc.getStart() + alleleList.get(0).length() - 1;
        final int eventLength = getEventLength(alleleList);

        final VariantContextBuilder builder = new VariantContextBuilder("UG_call", loc.getContig(), loc.getStart(), endLoc, alleleList);

        // create the genotypes; no-call everyone for now
        GenotypesContext genotypes = GenotypesContext.create();
        final int ploidy = UAC.genotypeArgs.samplePloidy;
        final List<Allele> noCall = GATKVariantContextUtils.noCallAlleles(ploidy);

        // For each sample, get genotype likelihoods based on pileup
        // compute prior likelihoods on haplotypes, and initialize haplotype likelihood matrix with them.

        for (Map.Entry<String, AlignmentContext> sample : contexts.entrySet()) {
            AlignmentContext context = AlignmentContextUtils.stratify(sample.getValue(), contextType);

            if (!perReadAlleleLikelihoodMap.containsKey(sample.getKey())){
                // no likelihoods have been computed for this sample at this site
                perReadAlleleLikelihoodMap.put(sample.getKey(), new PerReadAlleleLikelihoodMap());
            }
            final ReadBackedPileup pileup = context.getBasePileup();
            if (pileup != null) {
                final GenotypeBuilder b = new GenotypeBuilder(sample.getKey());
                final double[] genotypeLikelihoods = pairModel.computeDiploidReadHaplotypeLikelihoods(pileup, haplotypeMap, ref, eventLength, perReadAlleleLikelihoodMap.get(sample.getKey()), UAC.getSampleContamination().get(sample.getKey()));
                b.PL(genotypeLikelihoods);
                b.alleles(noCall);
                b.DP(getFilteredDepth(pileup));
                genotypes.add(b.make());

                if (DEBUG) {
                    System.out.format("Sample:%s Alleles:%s GL:", sample.getKey(), alleleList.toString());
                    for (int k = 0; k < genotypeLikelihoods.length; k++)
                        System.out.format("%1.4f ", genotypeLikelihoods[k]);
                    System.out.println();
                }
            }
        }

        return builder.genotypes(genotypes).make();
    }

    public static void getHaplotypeMapFromAlleles(final List<Allele> alleleList,
                                                  final ReferenceContext ref,
                                                  final GenomeLoc loc,
                                                  final LinkedHashMap<Allele, Haplotype> haplotypeMap) {
        // protect against having an indel too close to the edge of a contig
        if (loc.getStart() <= HAPLOTYPE_SIZE)
            haplotypeMap.clear();
            // check if there is enough reference window to create haplotypes (can be an issue at end of contigs)
        else if (ref.getWindow().getStop() < loc.getStop() + HAPLOTYPE_SIZE)
            haplotypeMap.clear();
        else if (alleleList.isEmpty())
            haplotypeMap.clear();
        else {
            final int eventLength = getEventLength(alleleList);
            final int hsize = ref.getWindow().size() - Math.abs(eventLength) - 1;
            final int numPrefBases = ref.getLocus().getStart() - ref.getWindow().getStart() + 1;

            if (hsize <= 0)  // protect against event lengths larger than ref window sizes
                haplotypeMap.clear();
            else
                haplotypeMap.putAll(Haplotype.makeHaplotypeListFromAlleles(alleleList, loc.getStart(),
                        ref, hsize, numPrefBases));
        }
    }

    public static int getEventLength(List<Allele> alleleList) {
        Allele refAllele = alleleList.get(0);
        Allele altAllele = alleleList.get(1);
        // look for alt allele that has biggest length distance to ref allele
        int maxLenDiff = 0;
        for (Allele a : alleleList) {
            if (a.isNonReference()) {
                int lenDiff = Math.abs(a.getBaseString().length() - refAllele.getBaseString().length());
                if (lenDiff > maxLenDiff) {
                    maxLenDiff = lenDiff;
                    altAllele = a;
                }
            }
        }

        return altAllele.getBaseString().length() - refAllele.getBaseString().length();

    }

    public static List<Allele> getInitialAlleleList(final RefMetaDataTracker tracker,
                                                    final ReferenceContext ref,
                                                    final Map<String, AlignmentContext> contexts,
                                                    final AlignmentContextUtils.ReadOrientation contextType,
                                                    final UnifiedArgumentCollection UAC,
                                                    final boolean ignoreSNPAllelesWhenGenotypingIndels) {

        List<Allele> alleles = new ArrayList<Allele>();
        if (UAC.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES) {
            VariantContext vc = null;
            for (final VariantContext vc_input : tracker.getValues(UAC.alleles, ref.getLocus())) {
                if (vc_input != null &&
                        allowableTypes.contains(vc_input.getType()) &&
                        ref.getLocus().getStart() == vc_input.getStart()) {
                    vc = vc_input;
                    break;
                }
            }
            // ignore places where we don't have a variant
            if (vc == null)
                return alleles;

            if (ignoreSNPAllelesWhenGenotypingIndels) {
                // if there's an allele that has same length as the reference (i.e. a SNP or MNP), ignore it and don't genotype it
                for (Allele a : vc.getAlleles())
                    if (a.isNonReference() && a.getBases().length == vc.getReference().getBases().length)
                        continue;
                    else
                        alleles.add(a);

            } else {
                alleles.addAll(vc.getAlleles());
            }

        } else {
            alleles = computeConsensusAlleles(ref, contexts, contextType, UAC);
        }
        return alleles;
    }

    // Overload function in GenotypeLikelihoodsCalculationModel so that, for an indel case, we consider a deletion as part of the pileup,
    // so that per-sample DP will include deletions covering the event.
    protected int getFilteredDepth(ReadBackedPileup pileup) {
        int count = 0;
        for (PileupElement p : pileup) {
            if (p.isDeletion() || BaseUtils.isRegularBase(p.getBase()))
                count++;
        }

        return count;
    }

}
