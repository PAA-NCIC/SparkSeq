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
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.ncic.bioinfo.sparkseq.algorithms.utils.AlignmentContextUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.AlignmentUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.BaseUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.MathUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.QualityUtils;
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
import org.ncic.bioinfo.sparkseq.exceptions.ReviewedGATKException;

import java.io.Serializable;
import java.util.*;

/**
 * Consistency of the site with strictly two segregating haplotypes
 *
 * <p>For diploid organisms, barring chromosomal abnormalities, we expect that any given sample has no more than 2 segregating haplotypes at a given site. If there is evidence for more
 * than 2 segregating haplotypes, the read data should be considered suspect and the evidence artifactual. Higher scores are indicative of regions with bad alignments, typically leading to artifactual SNP and indel calls.</p>
 *
 * <h3>Caveats</h3>
 * <p>HaplotypeCaller does not output this annotation because it already evaluates haplotype segregation internally. This annotation is only informative (and available) for variants called by Unified Genotyper.</p>
 * Author: wbc
 */
public class HaplotypeScore extends InfoFieldAnnotation implements StandardAnnotation, ActiveRegionBasedAnnotation {
    private final static boolean DEBUG = false;
    private final static int MIN_CONTEXT_WING_SIZE = 10;
    private final static int MAX_CONSENSUS_HAPLOTYPES_TO_CONSIDER = 50;
    private final static char REGEXP_WILDCARD = '.';

    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {
        if (vc.isSNP() && stratifiedContexts != null)
            return annotatePileup(ref, stratifiedContexts, vc);
        else
            return null;
    }

    private Map<String, Object> annotatePileup(final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc) {

        if (stratifiedContexts.size() == 0) // size 0 means that call was made by someone else and we have no data here
            return null;

        final AlignmentContext context = AlignmentContextUtils.joinContexts(stratifiedContexts.values());

        final int contextWingSize = Math.min((ref.getWindow().size() - 1) / 2, MIN_CONTEXT_WING_SIZE);
        final int contextSize = contextWingSize * 2 + 1;

        final int locus = ref.getLocus().getStart() + (ref.getLocus().getStop() - ref.getLocus().getStart()) / 2;

        final ReadBackedPileup pileup = context.getBasePileup();

        // Compute all haplotypes consistent with the current read pileup
        final List<Haplotype> haplotypes = computeHaplotypes(pileup, contextSize, locus, vc);

        final MathUtils.RunningAverage scoreRA = new MathUtils.RunningAverage();
        if (haplotypes != null) {
            for (final Genotype genotype : vc.getGenotypes()) {
                final AlignmentContext thisContext = stratifiedContexts.get(genotype.getSampleName());
                if (thisContext != null) {
                    final ReadBackedPileup thisPileup = thisContext.getBasePileup();
                    scoreRA.add(scoreReadsAgainstHaplotypes(haplotypes, thisPileup, contextSize, locus)); // Taking the simple average of all sample's score since the score can be negative and the RMS doesn't make sense
                }
            }
        }

        // annotate the score in the info field
        final Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%.4f", scoreRA.mean()));
        return map;
    }

    private static class HaplotypeComparator implements Comparator<Haplotype>, Serializable {

        public int compare(Haplotype a, Haplotype b) {
            if (a.getQualitySum() < b.getQualitySum())
                return 1;
            if (a.getQualitySum() > b.getQualitySum()) {
                return -1;
            }
            return 0;
        }
    }

    private List<Haplotype> computeHaplotypes(final ReadBackedPileup pileup, final int contextSize, final int locus, final VariantContext vc) {
        // Compute all possible haplotypes consistent with current pileup

        int haplotypesToCompute = vc.getAlternateAlleles().size() + 1;

        final PriorityQueue<Haplotype> candidateHaplotypeQueue = new PriorityQueue<Haplotype>(100, new HaplotypeComparator());
        final PriorityQueue<Haplotype> consensusHaplotypeQueue = new PriorityQueue<Haplotype>(MAX_CONSENSUS_HAPLOTYPES_TO_CONSIDER, new HaplotypeComparator());

        for (final PileupElement p : pileup) {
            final Haplotype haplotypeFromRead = getHaplotypeFromRead(p, contextSize, locus);
            if ( haplotypeFromRead != null )
                candidateHaplotypeQueue.add(haplotypeFromRead);
        }

        // Now that priority queue has been built with all reads at context, we need to merge and find possible segregating haplotypes
        Haplotype elem;
        while ((elem = candidateHaplotypeQueue.poll()) != null) {
            boolean foundHaplotypeMatch = false;
            Haplotype lastCheckedHaplotype = null;
            for (final Haplotype haplotypeFromList : consensusHaplotypeQueue) {
                final Haplotype consensusHaplotype = getConsensusHaplotype(elem, haplotypeFromList);
                if (consensusHaplotype != null) {
                    foundHaplotypeMatch = true;
                    if (consensusHaplotype.getQualitySum() > haplotypeFromList.getQualitySum()) {
                        consensusHaplotypeQueue.remove(haplotypeFromList);
                        consensusHaplotypeQueue.add(consensusHaplotype);
                    }
                    break;
                } else {
                    lastCheckedHaplotype = haplotypeFromList;
                }
            }

            if (!foundHaplotypeMatch && consensusHaplotypeQueue.size() < MAX_CONSENSUS_HAPLOTYPES_TO_CONSIDER) {
                consensusHaplotypeQueue.add(elem);
            } else if (!foundHaplotypeMatch && lastCheckedHaplotype != null && elem.getQualitySum() > lastCheckedHaplotype.getQualitySum()) {
                consensusHaplotypeQueue.remove(lastCheckedHaplotype);
                consensusHaplotypeQueue.add(elem);
            }
        }

        // Now retrieve the N most popular haplotypes
        if (consensusHaplotypeQueue.size() > 0) {
            // The consensus haplotypes are in a quality-ordered priority queue, so the best haplotypes are just the ones at the front of the queue
            final Haplotype haplotype1 = consensusHaplotypeQueue.poll();

            List<Haplotype> hlist = new ArrayList<Haplotype>();
            hlist.add(new Haplotype(haplotype1.getBases(), 60));

            for (int k = 1; k < haplotypesToCompute; k++) {
                Haplotype haplotype2 = consensusHaplotypeQueue.poll();
                if (haplotype2 == null) {
                    haplotype2 = haplotype1;
                } // Sometimes only the reference haplotype can be found
                hlist.add(new Haplotype(haplotype2.getBases(), 20));
            }
            return hlist;
        } else
            return null;
    }

    /**
     * Return a haplotype object constructed from the read or null if read's cigar is null
     *
     * @param p                pileup element representing the read
     * @param contextSize      the context size to use
     * @param locus            the position
     * @return possibly null Haplotype object constructed from the read
     */
    private Haplotype getHaplotypeFromRead(final PileupElement p, final int contextSize, final int locus) {
        final GATKSAMRecord read = p.getRead();
        if ( read.getCigar() == null )
            return null;

        final byte[] haplotypeBases = new byte[contextSize];
        Arrays.fill(haplotypeBases, (byte) REGEXP_WILDCARD);
        final byte[] baseQualities = new byte[contextSize];
        Arrays.fill(baseQualities, (byte)0);

        byte[] readBases = read.getReadBases();
        readBases = AlignmentUtils.readToAlignmentByteArray(read.getCigar(), readBases); // Adjust the read bases based on the Cigar string
        byte[] readQuals = read.getBaseQualities();
        readQuals = AlignmentUtils.readToAlignmentByteArray(read.getCigar(), readQuals); // Shift the location of the qual scores based on the Cigar string

        final int readOffsetFromPileup = AlignmentUtils.calcAlignmentByteArrayOffset(read.getCigar(), p, read.getAlignmentStart(), locus);
        final int baseOffsetStart = readOffsetFromPileup - (contextSize - 1) / 2;

        for (int i = 0; i < contextSize; i++) {
            final int baseOffset = i + baseOffsetStart;
            if (baseOffset < 0) {
                continue;
            }
            if (baseOffset >= readBases.length) {
                break;
            }
            if (readQuals[baseOffset] == PileupElement.DELETION_BASE) {
                readQuals[baseOffset] = PileupElement.DELETION_QUAL;
            }
            if (!BaseUtils.isRegularBase(readBases[baseOffset])) {
                readBases[baseOffset] = (byte) REGEXP_WILDCARD;
                readQuals[baseOffset] = (byte) 0;
            } // N's shouldn't be treated as distinct bases
            readQuals[baseOffset] = (byte) Math.min((int) readQuals[baseOffset], p.getMappingQual());
            if (((int) readQuals[baseOffset]) < 5) {
                readQuals[baseOffset] = (byte) 0;
            } // quals less than 5 are used as codes and don't have actual probabilistic meaning behind them
            haplotypeBases[i] = readBases[baseOffset];
            baseQualities[i] = readQuals[baseOffset];
        }

        return new Haplotype(haplotypeBases, baseQualities);
    }

    private Haplotype getConsensusHaplotype(final Haplotype haplotypeA, final Haplotype haplotypeB) {
        final byte[] a = haplotypeA.getBases();
        final byte[] b = haplotypeB.getBases();

        if (a.length != b.length) {
            throw new ReviewedGATKException("Haplotypes a and b must be of same length");
        }

        byte chA, chB;
        final byte wc = (byte) REGEXP_WILDCARD;

        final int length = a.length;
        final byte[] consensusChars = new byte[length];
        final int[] consensusQuals = new int[length];

        final int[] qualsA = haplotypeA.getQuals();
        final int[] qualsB = haplotypeB.getQuals();

        for (int i = 0; i < length; i++) {
            chA = a[i];
            chB = b[i];

            if ((chA != chB) && (chA != wc) && (chB != wc))
                return null;

            if ((chA == wc) && (chB == wc)) {
                consensusChars[i] = wc;
                consensusQuals[i] = 0;
            } else if ((chA == wc)) {
                consensusChars[i] = chB;
                consensusQuals[i] = qualsB[i];
            } else if ((chB == wc)) {
                consensusChars[i] = chA;
                consensusQuals[i] = qualsA[i];
            } else {
                consensusChars[i] = chA;
                consensusQuals[i] = qualsA[i] + qualsB[i];
            }
        }

        return new Haplotype(consensusChars, consensusQuals);
    }

    // calculate the haplotype scores by walking over all reads and comparing them to the haplotypes
    private double scoreReadsAgainstHaplotypes(final List<Haplotype> haplotypes, final ReadBackedPileup pileup, final int contextSize, final int locus) {
        if (DEBUG) System.out.printf("HAP1: %s%n", haplotypes.get(0));
        if (DEBUG) System.out.printf("HAP2: %s%n", haplotypes.get(1));

        final ArrayList<double[]> haplotypeScores = new ArrayList<double[]>();
        for (final PileupElement p : pileup) {
            // Score all the reads in the pileup, even the filtered ones
            final double[] scores = new double[haplotypes.size()];
            for (int i = 0; i < haplotypes.size(); i++) {
                final Haplotype haplotype = haplotypes.get(i);
                final double score = scoreReadAgainstHaplotype(p, contextSize, haplotype, locus);
                scores[i] = score;
                if (DEBUG) {
                    System.out.printf("  vs. haplotype %d = %f%n", i, score);
                }
            }
            haplotypeScores.add(scores);
        }

        double overallScore = 0.0;
        for (final double[] readHaplotypeScores : haplotypeScores) {
            overallScore += MathUtils.arrayMin(readHaplotypeScores);
        }

        return overallScore;
    }

    private double scoreReadAgainstHaplotype(final PileupElement p, final int contextSize, final Haplotype haplotype, final int locus) {
        double expected = 0.0;
        double mismatches = 0.0;

        final GATKSAMRecord read = p.getRead();
        if ( read.getCigar() == null )
            return 0.0;

        // What's the expected mismatch rate under the model that this read is actually sampled from
        // this haplotype?  Let's assume the consensus base c is a random choice one of A, C, G, or T, and that
        // the observed base is actually from a c with an error rate e.  Since e is the rate at which we'd
        // see a miscalled c, the expected mismatch rate is really e.  So the expected number of mismatches
        // is just sum_i e_i for i from 1..n for n sites
        //
        // Now, what's the probabilistic sum of mismatches?  Suppose that the base b is equal to c.  Well, it could
        // actually be a miscall in a matching direction, which would happen at a e / 3 rate.  If b != c, then
        // the chance that it is actually a mismatch is 1 - e, since any of the other 3 options would be a mismatch.
        // so the probability-weighted mismatch rate is sum_i ( matched ? e_i / 3 : 1 - e_i ) for i = 1 ... n
        final byte[] haplotypeBases = haplotype.getBases();
        byte[] readBases = read.getReadBases();

        readBases = AlignmentUtils.readToAlignmentByteArray(p.getRead().getCigar(), readBases); // Adjust the read bases based on the Cigar string
        byte[] readQuals = read.getBaseQualities();
        readQuals = AlignmentUtils.readToAlignmentByteArray(p.getRead().getCigar(), readQuals); // Shift the location of the qual scores based on the Cigar string
        int readOffsetFromPileup = AlignmentUtils.calcAlignmentByteArrayOffset(p.getRead().getCigar(), p, read.getAlignmentStart(), locus);
        final int baseOffsetStart = readOffsetFromPileup - (contextSize - 1) / 2;

        for (int i = 0; i < contextSize; i++) {
            final int baseOffset = i + baseOffsetStart;
            if (baseOffset < 0) {
                continue;
            }
            if (baseOffset >= readBases.length) {
                break;
            }

            final byte haplotypeBase = haplotypeBases[i];
            final byte readBase = readBases[baseOffset];

            final boolean matched = (readBase == haplotypeBase || haplotypeBase == (byte) REGEXP_WILDCARD);
            byte qual = readQuals[baseOffset];
            if (qual == PileupElement.DELETION_BASE) {
                qual = PileupElement.DELETION_QUAL;
            } // calcAlignmentByteArrayOffset fills the readQuals array with DELETION_BASE at deletions
            qual = (byte) Math.min((int) qual, p.getMappingQual());
            if (((int) qual) >= 5) { // quals less than 5 are used as codes and don't have actual probabilistic meaning behind them
                final double e = QualityUtils.qualToErrorProb(qual);
                expected += e;
                mismatches += matched ? e : 1.0 - e / 3.0;
            }

            // a more sophisticated calculation would include the reference quality, but it's nice to actually penalize
            // the mismatching of poorly determined regions of the consensus
        }

        return mismatches - expected;
    }

    public List<String> getKeyNames() {
        return Arrays.asList("HaplotypeScore");
    }

    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFInfoHeaderLine("HaplotypeScore", 1, VCFHeaderLineType.Float, "Consistency of the site with at most two segregating haplotypes"));
    }

    private static class Haplotype  {
        private final byte[] bases;
        private final int[] quals;
        private int qualitySum = -1;

        public Haplotype( final byte[] bases, final int[] quals ) {
            this.bases = bases;
            this.quals = quals;
        }

        public Haplotype( final byte[] bases, final int qual ) {
            this.bases = bases;
            quals = new int[bases.length];
            Arrays.fill(quals, qual);
        }

        public Haplotype( final byte[] bases, final byte[] quals ) {
            this.bases = bases;
            this.quals = new int[quals.length];
            for ( int i = 0 ; i < quals.length; i++ )
                this.quals[i] = (int)quals[i];
        }

        public double getQualitySum() {
            if ( qualitySum == -1 ) {
                qualitySum = 0;
                for ( final int qual : quals ) {
                    qualitySum += qual;
                }
            }
            return qualitySum;
        }

        public int[] getQuals() {
            return quals.clone();
        }

        public byte[] getBases() {
            return bases.clone();
        }
    }
}
