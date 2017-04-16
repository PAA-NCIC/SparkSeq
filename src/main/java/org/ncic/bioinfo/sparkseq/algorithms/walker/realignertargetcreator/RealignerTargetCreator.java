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
package org.ncic.bioinfo.sparkseq.algorithms.walker.realignertargetcreator;

import htsjdk.variant.variantcontext.VariantContext;
import org.ncic.bioinfo.sparkseq.algorithms.engine.LocusWalker;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLoc;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLocParser;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefMetaDataTracker;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.ReferenceContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.AlignmentContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.PileupElement;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.ReadBackedPileup;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.SamContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.BadCigarFilter;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.BadMateFilter;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.Filter;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.MappingQualityUnavailableFilter;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.MappingQualityZeroFilter;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RODContentProvider;

import java.util.ArrayList;
import java.util.List;

/**
 * Author: wbc
 */
public class RealignerTargetCreator extends LocusWalker {

    public static final int windowSize = 10;
    public static final int maxIntervalSize = 500;

    private EventPair sum;

    private List<GenomeLoc> result = new ArrayList<>();

    public RealignerTargetCreator(GenomeLocParser genomeLocParser,
                                  RefContentProvider refContentProvider,
                                  SamContentProvider samContentProvider,
                                  List<RODContentProvider> rodContentProviderList) {
        super(genomeLocParser, refContentProvider, samContentProvider, rodContentProviderList);
        initialize();
    }

    @Override
    protected void initialize() {
        sum = new EventPair(null, null);
    }

    @Override
    protected List<Filter> getFilter() {
        List<Filter> filters = new ArrayList<>();
        filters.add(new MappingQualityZeroFilter());
        filters.add(new MappingQualityUnavailableFilter());
        filters.add(new BadMateFilter());
        filters.add(new BadCigarFilter());
        return filters;
    }

    protected void map(RefMetaDataTracker tracker,
                       ReferenceContext ref,
                       AlignmentContext context) {
        boolean hasIndel = false;
        boolean hasInsertion = false;
        boolean hasPointEvent = false;

        int furthestStopPos = -1;

        // look at the rods for indels or SNPs
        if (tracker != null) {
            for (VariantContext vc : tracker.getValues(VariantContext.class)) {
                switch (vc.getType()) {
                    case INDEL:
                        hasIndel = true;
                        if (vc.isSimpleInsertion())
                            hasInsertion = true;
                        break;
                    case SNP:
                        hasPointEvent = true;
                        break;
                    case MIXED:
                        hasPointEvent = true;
                        hasIndel = true;
                        if (vc.isSimpleInsertion())
                            hasInsertion = true;
                        break;
                    default:
                        break;
                }
                if (hasIndel)
                    furthestStopPos = vc.getEnd();
            }
        }

        // look at the normal context to get deletions and positions with high entropy
        final ReadBackedPileup pileup = context.getBasePileup();

        for (PileupElement p : pileup) {

            // check the ends of the reads to see how far they extend
            furthestStopPos = Math.max(furthestStopPos, p.getRead().getAlignmentEnd());

            // is it a deletion or insertion?
            if (p.isDeletion() || p.isBeforeInsertion()) {
                hasIndel = true;
                if (p.isBeforeInsertion())
                    hasInsertion = true;
            }
        }

        // return null if no event occurred
        if (!hasIndel && !hasPointEvent)
            return;

        // return null if we didn't find any usable reads/rods associated with the event
        if (furthestStopPos == -1)
            return;

        GenomeLoc eventLoc = context.getLocation();
        if (hasInsertion)
            eventLoc = genomeLocParser.createGenomeLoc(eventLoc.getContig(), eventLoc.getStart(), eventLoc.getStart() + 1);

        EventType eventType = (hasIndel ? (hasPointEvent ? EventType.BOTH : EventType.INDEL_EVENT) : EventType.POINT_EVENT);

        Event event = new Event(eventLoc, furthestStopPos, eventType);
        reduce(event, sum);
    }

    private EventPair reduce(Event value, EventPair sum) {
        if (value == null) {
            ; // do nothing
        } else if (sum.left == null) {
            sum.left = value;
        } else if (sum.right == null) {
            if (canBeMerged(sum.left, value))
                sum.left = mergeEvents(sum.left, value);
            else
                sum.right = value;
        } else {
            if (canBeMerged(sum.right, value))
                sum.right = mergeEvents(sum.right, value);
            else {
                if (sum.right.isReportableEvent(genomeLocParser))
                    sum.intervals.add(sum.right.getLoc(genomeLocParser));
                sum.right = value;
            }
        }

        return sum;
    }

    @Override
    protected void onTraversalDone() {
        if (sum.left != null && sum.left.isReportableEvent(genomeLocParser))
            sum.intervals.add(sum.left.getLoc(genomeLocParser));
        if (sum.right != null && sum.right.isReportableEvent(genomeLocParser))
            sum.intervals.add(sum.right.getLoc(genomeLocParser));

        for (GenomeLoc loc : sum.intervals) {
            result.add(loc);
        }
    }

    public List<GenomeLoc> getTargetIntervals() {
        return result;
    }

    private boolean canBeMerged(Event left, Event right) {
        return left.loc.getContigIndex() == right.loc.getContigIndex() && left.furthestStopPos >= right.loc.getStart();
    }

    static private Event mergeEvents(Event left, Event right) {
        left.merge(right);
        return left;
    }
}
