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
package org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.log4j.Logger;
import org.ncic.bioinfo.sparkseq.algorithms.utils.AlignmentContextUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GATKVariantContextUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLoc;
import org.ncic.bioinfo.sparkseq.algorithms.utils.ReadUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.clip.ReadClipper;
import org.ncic.bioinfo.sparkseq.algorithms.data.basic.Pair;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.ReferenceContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.AlignmentContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.GATKSAMRecord;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.PileupElement;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.ReadBackedPileup;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Author: wbc
 */
public class ConsensusAlleleCounter {
    final protected static Logger logger = Logger.getLogger(ConsensusAlleleCounter.class);
    private final int minIndelCountForGenotyping;
    private final boolean doMultiAllelicCalls;
    private final double minFractionInOneSample;

    public ConsensusAlleleCounter(final boolean doMultiAllelicCalls,
                                  final int minIndelCountForGenotyping,
                                  final double minFractionInOneSample) {
        this.minIndelCountForGenotyping = minIndelCountForGenotyping;
        this.doMultiAllelicCalls = doMultiAllelicCalls;
        this.minFractionInOneSample = minFractionInOneSample;
    }

    /**
     * Returns a list of Alleles at this locus that may be segregating
     *
     * @param ref
     * @param contexts
     * @param contextType
     * @return
     */
    public List<Allele> computeConsensusAlleles(ReferenceContext ref,
                                                Map<String, AlignmentContext> contexts,
                                                AlignmentContextUtils.ReadOrientation contextType) {
        final Map<String, Integer> consensusIndelStrings = countConsensusAlleles(ref, contexts, contextType);
        return consensusCountsToAlleles(ref, consensusIndelStrings);
    }

    //
    // TODO -- WARNING DOESN'T WORK WITH REDUCED READS
    //
    private Map<String, Integer> countConsensusAlleles(ReferenceContext ref,
                                                       Map<String, AlignmentContext> contexts,
                                                       AlignmentContextUtils.ReadOrientation contextType) {
        final GenomeLoc loc = ref.getLocus();
        HashMap<String, Integer> consensusIndelStrings = new HashMap<String, Integer>();

        int insCount = 0, delCount = 0;
        // quick check of total number of indels in pileup
        for ( Map.Entry<String, AlignmentContext> sample : contexts.entrySet() ) {
            final AlignmentContext context = AlignmentContextUtils.stratify(sample.getValue(), contextType);

            final ReadBackedPileup indelPileup = context.getBasePileup();
            insCount += indelPileup.getNumberOfInsertionsAfterThisElement();
            delCount += indelPileup.getNumberOfDeletionsAfterThisElement();
        }

        if ( insCount < minIndelCountForGenotyping && delCount < minIndelCountForGenotyping )
            return Collections.emptyMap();

        for (Map.Entry<String, AlignmentContext> sample : contexts.entrySet()) {
            // todo -- warning, can be duplicating expensive partition here
            AlignmentContext context = AlignmentContextUtils.stratify(sample.getValue(), contextType);

            final ReadBackedPileup indelPileup = context.getBasePileup();

            final int nIndelReads = indelPileup.getNumberOfInsertionsAfterThisElement() + indelPileup.getNumberOfDeletionsAfterThisElement();
            final int nReadsOverall = indelPileup.getNumberOfElements();

            if ( nIndelReads == 0 || (nIndelReads / (1.0 * nReadsOverall)) < minFractionInOneSample ) {
                continue;
            }

            for (PileupElement p : indelPileup) {
                final GATKSAMRecord read = ReadClipper.hardClipAdaptorSequence(p.getRead());
                if (read == null)
                    continue;
                if (ReadUtils.is454Read(read)) {
                    continue;
                }

                if ( p.isBeforeInsertion() ) {
                    final String insertionBases = p.getBasesOfImmediatelyFollowingInsertion();
                    // edge case: ignore a deletion immediately preceding an insertion as p.getBasesOfImmediatelyFollowingInsertion() returns null [EB]
                    if ( insertionBases == null )
                        continue;

                    boolean foundKey = false;
                    // copy of hashmap into temp arrayList
                    ArrayList<Pair<String,Integer>> cList = new ArrayList<Pair<String,Integer>>();
                    for (Map.Entry<String, Integer> s : consensusIndelStrings.entrySet()) {
                        cList.add(new Pair<String, Integer>(s.getKey(), s.getValue()));
                    }

                    if (read.getAlignmentEnd() == loc.getStart()) {
                        // first corner condition: a read has an insertion at the end, and we're right at the insertion.
                        // In this case, the read could have any of the inserted bases and we need to build a consensus

                        for (int k=0; k < cList.size(); k++) {
                            String s = cList.get(k).getFirst();
                            int cnt = cList.get(k).getSecond();
                            // case 1: current insertion is prefix of indel in hash map
                            if (s.startsWith(insertionBases)) {
                                cList.set(k,new Pair<String, Integer>(s,cnt+1));
                                foundKey = true;
                            }
                            else if (insertionBases.startsWith(s)) {
                                // case 2: indel stored in hash table is prefix of current insertion
                                // In this case, new bases are new key.
                                foundKey = true;
                                cList.set(k,new Pair<String, Integer>(insertionBases,cnt+1));
                            }
                        }
                        if (!foundKey)
                            // none of the above: event bases not supported by previous table, so add new key
                            cList.add(new Pair<String, Integer>(insertionBases,1));

                    }
                    else if (read.getAlignmentStart() == loc.getStart()+1) {
                        // opposite corner condition: read will start at current locus with an insertion
                        for (int k=0; k < cList.size(); k++) {
                            String s = cList.get(k).getFirst();
                            int cnt = cList.get(k).getSecond();
                            if (s.endsWith(insertionBases)) {
                                // case 1: current insertion (indelString) is suffix of indel in hash map (s)
                                cList.set(k,new Pair<String, Integer>(s,cnt+1));
                                foundKey = true;
                            }
                            else if (insertionBases.endsWith(s)) {
                                // case 2: indel stored in hash table is prefix of current insertion
                                // In this case, new bases are new key.
                                foundKey = true;
                                cList.set(k,new Pair<String, Integer>(insertionBases,cnt+1));
                            }
                        }
                        if (!foundKey)
                            // none of the above: event bases not supported by previous table, so add new key
                            cList.add(new Pair<String, Integer>(insertionBases,1));


                    }
                    else {
                        // normal case: insertion somewhere in the middle of a read: add count to arrayList
                        int cnt = consensusIndelStrings.containsKey(insertionBases)? consensusIndelStrings.get(insertionBases):0;
                        cList.add(new Pair<String, Integer>(insertionBases,cnt+1));
                    }

                    // copy back arrayList into hashMap
                    consensusIndelStrings.clear();
                    for (Pair<String,Integer> pair : cList) {
                        consensusIndelStrings.put(pair.getFirst(),pair.getSecond());
                    }

                }
                else if ( p.isBeforeDeletionStart() ) {
                    final String deletionString = String.format("D%d",p.getLengthOfImmediatelyFollowingIndel());
                    int cnt = consensusIndelStrings.containsKey(deletionString)? consensusIndelStrings.get(deletionString):0;
                    consensusIndelStrings.put(deletionString,cnt+1);
                }
            }
        }

        return consensusIndelStrings;
    }

    private List<Allele> consensusCountsToAlleles(final ReferenceContext ref,
                                                  final Map<String, Integer> consensusIndelStrings) {
        final GenomeLoc loc = ref.getLocus();
        final Collection<VariantContext> vcs = new ArrayList<VariantContext>();
        int maxAlleleCnt = 0;
        Allele refAllele, altAllele;

        for (final Map.Entry<String, Integer> elt : consensusIndelStrings.entrySet()) {
            final String s = elt.getKey();
            final int curCnt = elt.getValue();
            int stop = 0;

            // if observed count if above minimum threshold, we will genotype this allele
            if (curCnt < minIndelCountForGenotyping)
                continue;

            if (s.startsWith("D")) {
                // get deletion length
                final int dLen = Integer.valueOf(s.substring(1));
                // get ref bases of accurate deletion
                final int startIdxInReference = 1 + loc.getStart() - ref.getWindow().getStart();
                stop = loc.getStart() + dLen;
                final byte[] refBases = Arrays.copyOfRange(ref.getBases(), startIdxInReference - 1, startIdxInReference + dLen);   // add reference padding

                if (Allele.acceptableAlleleBases(refBases, false)) {
                    refAllele = Allele.create(refBases, true);
                    altAllele = Allele.create(ref.getBase(), false);
                }
                else continue; // don't go on with this allele if refBases are non-standard
            } else {
                // insertion case
                final String insertionBases = (char)ref.getBase() + s;  // add reference padding
                if (Allele.acceptableAlleleBases(insertionBases, false)) { // don't allow N's in insertions
                    refAllele = Allele.create(ref.getBase(), true);
                    altAllele = Allele.create(insertionBases, false);
                    stop = loc.getStart();
                }
                else continue; // go on to next allele if consensus insertion has any non-standard base.
            }


            final VariantContextBuilder builder = new VariantContextBuilder().source("");
            builder.loc(loc.getContig(), loc.getStart(), stop);
            builder.alleles(Arrays.asList(refAllele, altAllele));
            builder.noGenotypes();
            if (doMultiAllelicCalls) {
                vcs.add(builder.make());
                if (vcs.size() >= GenotypeLikelihoods.MAX_ALT_ALLELES_THAT_CAN_BE_GENOTYPED)
                    break;
            } else if (curCnt > maxAlleleCnt) {
                maxAlleleCnt = curCnt;
                vcs.clear();
                vcs.add(builder.make());
            }
        }

        if (vcs.isEmpty())
            return Collections.emptyList(); // nothing else to do, no alleles passed minimum count criterion

        final VariantContext mergedVC = GATKVariantContextUtils.simpleMerge(vcs, null, GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED, GATKVariantContextUtils.GenotypeMergeType.UNSORTED, false, false, null, false, false);
        return mergedVC.getAlleles();
    }
}
