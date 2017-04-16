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

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.commons.lang.ArrayUtils;
import org.apache.log4j.Logger;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLoc;
import org.ncic.bioinfo.sparkseq.algorithms.utils.haplotype.Haplotype;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.readlikelihood.ReadLikelihoods;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;

/**
 * Author: wbc
 */
public class LDMerger extends MergeVariantsAcrossHaplotypes {
    private final static Logger logger = Logger.getLogger(LDMerger.class);

    private final boolean DEBUG;
    private final int minSamplesToMergeSNPs;
    private final int minSamplesToMergeOtherEvents;

    public LDMerger(boolean DEBUG, int minSamplesToMergeSNPs, int minSamplesToMergeOtherEvents) {
        super();
        this.DEBUG = DEBUG;
        this.minSamplesToMergeSNPs = minSamplesToMergeSNPs;
        this.minSamplesToMergeOtherEvents = minSamplesToMergeOtherEvents;
    }

    protected LDMerger() {
        this(false, 1, 1);
    }

    // TODO -- should be class arguments and static variables in HC
    protected final static int MAX_DISTANCE_BETWEEN_SNPS_TO_MERGE = 6;
    protected final static int MAX_DISTANCE_BETWEEN_OTHER_EVENTS_TO_MERGE = 25;

    /**
     * We require 99% confidence that only the phased haplotypes exist in the population to merge the records
     */
    protected final static double MERGE_EVENTS_PROB_PHASED_THRESHOLD = 0.99;

    /**
     * Merge as many events among the haplotypes as possible based on pairwise LD among variants
     *
     * @param haplotypes a list of haplotypes whose events we want to merge
     * @param readLikelihoods map from sample name -> read likelihoods for each haplotype
     * @param startPosKeySet a set of starting positions of all events among the haplotypes
     * @param ref the reference bases
     * @param refLoc the span of the reference bases
     */
    @Override
    public boolean merge( final List<Haplotype> haplotypes,
                          final ReadLikelihoods<Haplotype> readLikelihoods,
                          final TreeSet<Integer> startPosKeySet,
                          final byte[] ref,
                          final GenomeLoc refLoc ) {
        if ( haplotypes == null ) throw new IllegalArgumentException("haplotypes cannot be null");
        if ( readLikelihoods == null ) throw new IllegalArgumentException("readLikelihoods cannot be null");
        if ( startPosKeySet == null ) throw new IllegalArgumentException("startPosKeySet cannot be null");
        if ( ref == null ) throw new IllegalArgumentException("ref cannot be null");
        if ( refLoc == null ) throw new IllegalArgumentException("refLoc cannot be null");
        if ( refLoc.size() != ref.length ) throw new IllegalArgumentException("refLoc size " + refLoc.size() + " != ref.length " + ref.length + " at " + refLoc);

        if( startPosKeySet.size() <= 1 ) { return false; }

        final int nSamples = readLikelihoods.sampleCount();
        final HaplotypeLDCalculator r2Calculator = new HaplotypeLDCalculator(haplotypes, readLikelihoods);
        boolean somethingWasMerged = false;
        boolean mapWasUpdated = true;
        while( mapWasUpdated ) {
            mapWasUpdated = mergeConsecutiveEventsBasedOnLDOnce(haplotypes, r2Calculator, nSamples, startPosKeySet, ref, refLoc);
            somethingWasMerged |= mapWasUpdated;
        }
        return somethingWasMerged;
    }

    /**
     * Merge the next pair of events, if possible
     *
     * @param haplotypes a list of haplotypes whose events we want to merge
     * @param ldCalculator calculates R^2 for pairs of events on demand
     * @param startPosKeySet a set of starting positions of all events among the haplotypes
     * @param ref the reference bases
     * @param refLoc the span of the reference bases
     * @return true if something was merged, false otherwise
     */
    protected boolean mergeConsecutiveEventsBasedOnLDOnce( final List<Haplotype> haplotypes,
                                                           final HaplotypeLDCalculator ldCalculator,
                                                           final int nSamples,
                                                           final TreeSet<Integer> startPosKeySet,
                                                           final byte[] ref,
                                                           final GenomeLoc refLoc ) {
        // loop over the set of start locations and consider pairs that start near each other
        final Iterator<Integer> iter = startPosKeySet.iterator();
        int thisStart = iter.next();
        while( iter.hasNext() ) {
            final int nextStart = iter.next();
            final LDMergeData toMerge = getPairOfEventsToMerge(haplotypes, thisStart, nextStart);

            if ( toMerge.canBeMerged(nSamples) ) {
                final double pPhased = ldCalculator.computeProbOfBeingPhased(toMerge.firstVC, toMerge.secondVC);

                if( DEBUG ) {
                    logger.info("Found consecutive biallelic events with R^2 = " + String.format("%.4f", pPhased));
                    logger.info("-- " + toMerge.firstVC);
                    logger.info("-- " + toMerge.secondVC);
                }

                if( pPhased > MERGE_EVENTS_PROB_PHASED_THRESHOLD) {
                    final VariantContext mergedVC = createMergedVariantContext(toMerge.firstVC, toMerge.secondVC, ref, refLoc);
                    // if for some reason the merging resulting in a bad allele, mergedVC will be null, and we will just remove first and second
                    replaceVariantContextsInMap(haplotypes, startPosKeySet, mergedVC, toMerge.firstVC, toMerge.secondVC);
                    return true; // break out of tree set iteration since it was just updated, start over from the beginning and keep merging events
                }
            }

            thisStart = nextStart;
        }

        return false;
    }

    /**
     * Info about potential LD merge of two variant contexts
     */
    private class LDMergeData {
        VariantContext firstVC = null, secondVC = null;
        boolean canBeMerged = true;

        /** Tell this object that it cant be merged for some reason */
        public LDMergeData cantBeMerged() {
            canBeMerged = false;
            return this;
        }

        /**
         * Can these two events be merged
         * @param nSamples the number of samples we're considering
         * @return true if we can merge our two variant contexts
         */
        public boolean canBeMerged(final int nSamples) {
            if ( ! canBeMerged || firstVC == null || secondVC == null )
                return false;

            final int distance = secondVC.getStart() - firstVC.getEnd();
            if ( firstVC.isSNP() && secondVC.isSNP() ) {
                return nSamples >= minSamplesToMergeSNPs && distance <= MAX_DISTANCE_BETWEEN_SNPS_TO_MERGE;
            } else {
                return nSamples >= minSamplesToMergeOtherEvents && distance <= MAX_DISTANCE_BETWEEN_OTHER_EVENTS_TO_MERGE;
            }
        }
    }

    /**
     * Get the information about the potential merge of two events starting at thisStart and nextStart
     * @param haplotypes our haplotypes
     * @param thisStart the starting position of the first event to merge
     * @param nextStart the starting position of the next event to merge
     * @return never {@code null}.
     */
    private LDMergeData getPairOfEventsToMerge(final List<Haplotype> haplotypes, final int thisStart, final int nextStart) {
        final LDMergeData mergeData = new LDMergeData();

        for( final Haplotype h : haplotypes ) {
            // only make complex substitutions out of consecutive biallelic sites
            final VariantContext thisHapVC = h.getEventMap().get(thisStart);
            if( thisHapVC != null && !thisHapVC.isSymbolic() ) { // something was found at this location on this haplotype
                if( mergeData.firstVC == null ) {
                    mergeData.firstVC = thisHapVC;
                } else if( !thisHapVC.hasSameAllelesAs( mergeData.firstVC) ) {
                    return mergeData.cantBeMerged();
                }
            }
            final VariantContext nextHapVC = h.getEventMap().get(nextStart);
            if( nextHapVC != null && !nextHapVC.isSymbolic() ) { // something was found at the next location on this haplotype
                if( mergeData.secondVC == null ) {
                    mergeData.secondVC = nextHapVC;
                } else if( !nextHapVC.hasSameAllelesAs( mergeData.secondVC) ) {
                    return mergeData.cantBeMerged();
                }
            }
        }

        // don't try to merge overlapping events
        if ( mergeData.firstVC != null && mergeData.secondVC != null && mergeData.firstVC.getEnd() >= mergeData.secondVC.getStart() )
            return mergeData.cantBeMerged();

        return mergeData;
    }

    // BUGBUG: make this merge function more general
    protected VariantContext createMergedVariantContext( final VariantContext thisVC, final VariantContext nextVC, final byte[] ref, final GenomeLoc refLoc ) {
        final int thisStart = thisVC.getStart();
        final int nextStart = nextVC.getStart();
        byte[] refBases = new byte[]{};
        byte[] altBases = new byte[]{};
        refBases = ArrayUtils.addAll(refBases, thisVC.getReference().getBases());
        altBases = ArrayUtils.addAll(altBases, thisVC.getAlternateAllele(0).getBases());
        int locus;
        for( locus = thisStart + refBases.length; locus < nextStart; locus++ ) {
            final byte refByte = ref[locus - refLoc.getStart()];
            refBases = ArrayUtils.add(refBases, refByte);
            altBases = ArrayUtils.add(altBases, refByte);
        }
        refBases = ArrayUtils.addAll(refBases, ArrayUtils.subarray(nextVC.getReference().getBases(), locus > nextStart ? 1 : 0, nextVC.getReference().getBases().length)); // special case of deletion including the padding base of consecutive indel
        altBases = ArrayUtils.addAll(altBases, nextVC.getAlternateAllele(0).getBases());

        int iii = 0;
        if( refBases.length == altBases.length ) { // insertion + deletion of same length creates an MNP --> trim common prefix bases off the beginning of the allele
            while( iii < refBases.length && refBases[iii] == altBases[iii] ) { iii++; }
            if ( iii == refBases.length ) {
                // we've become a null allele, such as with CA/C + A/AA -> CA/CA => after trimming there's nothing left
                // so return a null variant context so we can eliminate the variants from consideration
                return null;
            }
        }


        final Allele refAllele = Allele.create( ArrayUtils.subarray(refBases, iii, refBases.length), true );
        final Allele altAllele =  Allele.create( ArrayUtils.subarray(altBases, iii, altBases.length), false );
        return new VariantContextBuilder("merged", thisVC.getChr(), thisVC.getStart() + iii, nextVC.getEnd(), Arrays.asList(refAllele, altAllele)).make();
    }

    /**
     * Update the event maps in all haplotypes to replace a replacement of update1 and 2 with replacement
     *
     * @param haplotypes the haplotypes whose event maps we need to update
     * @param startPosKeySet a sorted set of start positions that we must update
     * @param replacement a VariantContext to replace update1 and update2 with.  Can be null, indicating that we just want to remove update1 and update2
     * @param update1 the first VC we want to update
     * @param update2 the second VC we want to update
     */
    private void replaceVariantContextsInMap(final List<Haplotype> haplotypes,
                                             final TreeSet<Integer> startPosKeySet,
                                             final VariantContext replacement,
                                             final VariantContext update1, final VariantContext update2) {
        // remove the old event from the eventMap on every haplotype and the start pos key set, replace with merged event
        for( final Haplotype h : haplotypes ) {
            // if we had both events, add replacement.  In some cases the haplotype may not have both
            // events but they were still merged because the haplotype isn't a particularly informative
            // haplotype in any case.  The order of operations here is important because we are modifying the map
            final boolean shouldAdd = h.getEventMap().containsKey(update1.getStart()) && h.getEventMap().containsKey(update2.getStart());
            h.getEventMap().remove(update1.getStart());
            h.getEventMap().remove(update2.getStart());
            if ( shouldAdd && replacement != null ) {
                h.getEventMap().addVC(replacement, false); // cannot merge we other events at the same position
            }
        }

        startPosKeySet.remove(update1.getStart());
        startPosKeySet.remove(update2.getStart());
        if ( replacement != null ) startPosKeySet.add(replacement.getStart());
    }
}
