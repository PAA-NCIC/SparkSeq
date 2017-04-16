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

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.log4j.Logger;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLocParser;
import org.ncic.bioinfo.sparkseq.algorithms.utils.QualityUtils;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefMetaDataTracker;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.ReferenceContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.AlignmentContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.GATKSAMRecord;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.PileupElement;
import org.ncic.bioinfo.sparkseq.algorithms.utils.Utils;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.MostLikelyAllele;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.PerReadAlleleLikelihoodMap;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.interfaces.AnnotatorCompatible;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.interfaces.InfoFieldAnnotation;

import java.util.ArrayList;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;

/**
 * Author: wbc
 */
public abstract class StrandBiasTest extends InfoFieldAnnotation {
    private final static Logger logger = Logger.getLogger(StrandBiasTest.class);

    @Override
    public void initialize(final AnnotatorCompatible walker, final GenomeLocParser parser, final Set<VCFHeaderLine> headerLines) {
        boolean hasSBBSannotation = false;
        for (final VCFHeaderLine line : headerLines) {
            if (line instanceof VCFFormatHeaderLine) {
                final VCFFormatHeaderLine formatline = (VCFFormatHeaderLine) line;
                if (formatline.getID().equals(StrandBiasBySample.STRAND_BIAS_BY_SAMPLE_KEY_NAME)) {
                    hasSBBSannotation = true;
                    break;
                }
            }
        }

        if (hasSBBSannotation) {
            logger.info("StrandBiasBySample annotation exists in input VCF header. Attempting to use StrandBiasBySample " +
                    "values to calculate strand bias annotation values. If no sample has the SB genotype annotation, annotation may still fail.");
            return;
        }

    }

    @Override
    //template method for calculating strand bias annotations using the three different methods
    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {
        if (!vc.isVariant())
            return null;

        if (vc.hasGenotypes()) {
            boolean hasSB = false;
            for (final Genotype g : vc.getGenotypes()) {
                if (g.hasAnyAttribute(StrandBiasBySample.STRAND_BIAS_BY_SAMPLE_KEY_NAME)) {
                    hasSB = true;
                    break;
                }
            }
            if (hasSB)
                return calculateAnnotationFromGTfield(vc.getGenotypes());
        }

        //stratifiedContexts can come come from VariantAnnotator, but will be size 0 if no reads were provided
        if (vc.isSNP() && stratifiedContexts != null && stratifiedContexts.size() > 0) {
            return calculateAnnotationFromStratifiedContexts(stratifiedContexts, vc);
        }

        //stratifiedPerReadAllelelikelihoodMap can come from HaplotypeCaller call to VariantAnnotatorEngine
        else if (stratifiedPerReadAlleleLikelihoodMap != null) {
            return calculateAnnotationFromLikelihoodMap(stratifiedPerReadAlleleLikelihoodMap, vc);
        } else
            // for non-snp variants, we  need per-read likelihoods.
            // for snps, we can get same result from simple pileup
            return null;
    }

    protected abstract Map<String, Object> calculateAnnotationFromGTfield(final GenotypesContext genotypes);

    protected abstract Map<String, Object> calculateAnnotationFromStratifiedContexts(final Map<String, AlignmentContext> stratifiedContexts,
                                                                                     final VariantContext vc);

    protected abstract Map<String, Object> calculateAnnotationFromLikelihoodMap(final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap,
                                                                                final VariantContext vc);

    /**
     * Create the contingency table by retrieving the per-sample strand bias annotation and adding them together
     *
     * @param genotypes the genotypes from which to pull out the per-sample strand bias annotation
     * @param minCount  minimum threshold for the sample strand bias counts for each ref and alt.
     *                  If both ref and alt counts are above minCount the whole sample strand bias is added to the resulting table
     * @return the table used for several strand bias tests, will be null if none of the genotypes contain the per-sample SB annotation
     */
    protected int[][] getTableFromSamples(final GenotypesContext genotypes, final int minCount) {
        if (genotypes == null) {
            throw new IllegalArgumentException("Genotypes cannot be null.");
        }

        final int[] sbArray = {0, 0, 0, 0}; // reference-forward-reverse -by- alternate-forward-reverse
        boolean foundData = false;

        for (final Genotype g : genotypes) {
            if (g.isNoCall() || !g.hasAnyAttribute(StrandBiasBySample.STRAND_BIAS_BY_SAMPLE_KEY_NAME))
                continue;

            foundData = true;
            final ArrayList<Integer> sbbsString = (ArrayList<Integer>) g.getAnyAttribute(StrandBiasBySample.STRAND_BIAS_BY_SAMPLE_KEY_NAME);
            final int[] data = Utils.list2Array(sbbsString);
            if (passesMinimumThreshold(data, minCount)) {
                for (int index = 0; index < sbArray.length; index++) {
                    sbArray[index] += data[index];
                }
            }
        }

        return (foundData ? decodeSBBS(sbArray) : null);
    }

    /**
     * Allocate and fill a 2x2 strand contingency table.  In the end, it'll look something like this:
     * fw      rc
     * allele1   #       #
     * allele2   #       #
     *
     * @return a 2x2 contingency table
     */
    protected static int[][] getSNPContingencyTable(final Map<String, AlignmentContext> stratifiedContexts,
                                                    final Allele ref,
                                                    final Allele alt,
                                                    final int minQScoreToConsider,
                                                    final int minCount) {
        int[][] table = new int[2][2];

        for (final Map.Entry<String, AlignmentContext> sample : stratifiedContexts.entrySet()) {
            final int[] myTable = new int[4];
            for (final PileupElement p : sample.getValue().getBasePileup()) {

                if (!isUsableBase(p)) // ignore deletions and bad MQ
                    continue;

                if (p.getQual() < minQScoreToConsider || p.getMappingQual() < minQScoreToConsider)
                    continue;

                updateTable(myTable, Allele.create(p.getBase(), false), p.getRead(), ref, alt);
            }

            if (passesMinimumThreshold(myTable, minCount))
                copyToMainTable(myTable, table);
        }

        return table;
    }

    /**
     * Allocate and fill a 2x2 strand contingency table.  In the end, it'll look something like this:
     * fw      rc
     * allele1   #       #
     * allele2   #       #
     *
     * @return a 2x2 contingency table
     */
    public static int[][] getContingencyTable(final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap,
                                              final VariantContext vc,
                                              final int minCount) {
        if (stratifiedPerReadAlleleLikelihoodMap == null) {
            throw new IllegalArgumentException("stratifiedPerReadAlleleLikelihoodMap cannot be null");
        }
        if (vc == null) {
            throw new IllegalArgumentException("input vc cannot be null");
        }

        final Allele ref = vc.getReference();
        final Allele alt = vc.getAltAlleleWithHighestAlleleCount();
        final int[][] table = new int[2][2];

        for (final PerReadAlleleLikelihoodMap maps : stratifiedPerReadAlleleLikelihoodMap.values()) {
            final int[] myTable = new int[4];
            for (final Map.Entry<GATKSAMRecord, Map<Allele, Double>> el : maps.getLikelihoodReadMap().entrySet()) {
                final MostLikelyAllele mostLikelyAllele = PerReadAlleleLikelihoodMap.getMostLikelyAllele(el.getValue());
                final GATKSAMRecord read = el.getKey();
                updateTable(myTable, mostLikelyAllele.getAlleleIfInformative(), read, ref, alt);
            }
            if (passesMinimumThreshold(myTable, minCount))
                copyToMainTable(myTable, table);
        }

        return table;
    }

    /**
     * Helper method to copy the per-sample table to the main table
     *
     * @param perSampleTable per-sample table (single dimension)
     * @param mainTable      main table (two dimensions)
     */
    private static void copyToMainTable(final int[] perSampleTable, final int[][] mainTable) {
        mainTable[0][0] += perSampleTable[0];
        mainTable[0][1] += perSampleTable[1];
        mainTable[1][0] += perSampleTable[2];
        mainTable[1][1] += perSampleTable[3];
    }


    /**
     * Can the base in this pileup element be used in comparative tests?
     *
     * @param p the pileup element to consider
     * @return true if this base is part of a meaningful read for comparison, false otherwise
     */
    private static boolean isUsableBase(final PileupElement p) {
        return !(p.isDeletion() ||
                p.getMappingQual() == 0 ||
                p.getMappingQual() == QualityUtils.MAPPING_QUALITY_UNAVAILABLE ||
                ((int) p.getQual()) < QualityUtils.MIN_USABLE_Q_SCORE);
    }

    private static void updateTable(final int[] table, final Allele allele, final GATKSAMRecord read, final Allele ref, final Allele alt) {

        final boolean matchesRef = allele.equals(ref, true);
        final boolean matchesAlt = allele.equals(alt, true);

        if (matchesRef || matchesAlt) {
            final int offset = matchesRef ? 0 : 2;

            if (read.isStrandless()) {
                // a strandless read counts as observations on both strand, at 50% weight, with a minimum of 1
                // (the 1 is to ensure that a strandless read always counts as an observation on both strands, even
                // if the read is only seen once, because it's a merged read or other)
                table[offset]++;
                table[offset + 1]++;
            } else {
                // a normal read with an actual strand
                final boolean isFW = !read.getReadNegativeStrandFlag();
                table[offset + (isFW ? 0 : 1)]++;
            }
        }
    }

    /**
     * Does this strand data array pass the minimum threshold for inclusion?
     *
     * @param data the array
     * @return true if it passes the minimum threshold, false otherwise
     * @minCount The minimum threshold of counts in the array
     */
    protected static boolean passesMinimumThreshold(final int[] data, final int minCount) {
        // the ref and alt totals must be greater than MIN_COUNT
        return data[0] + data[1] + data[2] + data[3] > minCount;
    }

    /**
     * Helper function to parse the genotype annotation into the SB annotation array
     *
     * @param string the string that is returned by genotype.getAnnotation("SB")
     * @return the array used by the per-sample Strand Bias annotation
     */
    private static int[] encodeSBBS(final String string) {
        final int[] array = new int[4];
        final StringTokenizer tokenizer = new StringTokenizer(string, ",", false);
        for (int index = 0; index < 4; index++) {
            array[index] = Integer.parseInt(tokenizer.nextToken());
        }
        return array;
    }

    /**
     * Helper function to turn the  SB annotation array into a contingency table
     *
     * @param array the array used by the per-sample Strand Bias annotation
     * @return the table used by the StrandOddsRatio annotation
     */
    private static int[][] decodeSBBS(final int[] array) {
        if (array.length != 4) {
            throw new IllegalArgumentException("Expecting a length = 4 strand bias array.");
        }
        final int[][] table = new int[2][2];
        table[0][0] = array[0];
        table[0][1] = array[1];
        table[1][0] = array[2];
        table[1][1] = array[3];
        return table;
    }
}
