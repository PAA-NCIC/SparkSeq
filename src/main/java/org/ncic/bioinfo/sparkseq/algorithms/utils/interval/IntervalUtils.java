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
package org.ncic.bioinfo.sparkseq.algorithms.utils.interval;

import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLoc;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLocParser;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLocSortedSet;
import org.ncic.bioinfo.sparkseq.algorithms.utils.Utils;
import org.ncic.bioinfo.sparkseq.exceptions.UserException;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
 * Author: wbc
 */
public class IntervalUtils {

    /**
     * Returns true if the interval string is the "unmapped" interval
     * @param interval Interval to check
     * @return true if the interval string is the "unmapped" interval
     */
    public static boolean isUnmapped(String interval) {
        return (interval != null && interval.trim().toLowerCase().equals("unmapped"));
    }

    /**
     * merge two interval lists, using an interval set rule
     * @param setOne a list of genomeLocs, in order (cannot be NULL)
     * @param setTwo a list of genomeLocs, also in order (cannot be NULL)
     * @param rule the rule to use for merging, i.e. union, intersection, etc
     * @return a list, correctly merged using the specified rule
     */
    public static List<GenomeLoc> mergeListsBySetOperator(List<GenomeLoc> setOne, List<GenomeLoc> setTwo, IntervalSetRule rule) {
        // shortcut, if either set is zero, return the other set
        if (setOne == null || setOne.size() == 0 || setTwo == null || setTwo.size() == 0)
            return Collections.unmodifiableList((setOne == null || setOne.size() == 0) ? setTwo : setOne);

        // our master list, since we can't guarantee removal time in a generic list
        LinkedList<GenomeLoc> retList = new LinkedList<GenomeLoc>();

        // if we're set to UNION, just add them all
        if (rule == null || rule == IntervalSetRule.UNION) {
            retList.addAll(setOne);
            retList.addAll(setTwo);
            return Collections.unmodifiableList(retList);
        }

        // else we're INTERSECTION, create two indexes into the lists
        int iOne = 0;
        int iTwo = 0;

        // merge the second into the first using the rule
        while (iTwo < setTwo.size() && iOne < setOne.size())
            // if the first list is ahead, drop items off the second until we overlap
            if (setTwo.get(iTwo).isBefore(setOne.get(iOne)))
                iTwo++;
                // if the second is ahead, drop intervals off the first until we overlap
            else if (setOne.get(iOne).isBefore(setTwo.get(iTwo)))
                iOne++;
                // we overlap, intersect the two intervals and add the result.  Then remove the interval that ends first.
            else {
                retList.add(setOne.get(iOne).intersect(setTwo.get(iTwo)));
                if (setOne.get(iOne).getStop() < setTwo.get(iTwo).getStop()) iOne++;
                else iTwo++;
            }

        //if we have an empty list, throw an exception.  If they specified intersection and there are no items, this is bad.
        if (retList.size() == 0)
            throw new UserException.BadInput("The INTERSECTION of your -L options produced no intervals.");

        // we don't need to add the rest of remaining locations, since we know they don't overlap. return what we have
        return Collections.unmodifiableList(retList);
    }

    /**
     * Sorts and merges an interval list.  Multiple techniques are available for merging: ALL, which combines
     * all overlapping and abutting intervals into an interval that spans the union of all covered bases, and
     * OVERLAPPING_ONLY, which unions overlapping intervals but keeps abutting intervals separate.
     *
     * @param parser Genome loc parser for the intervals.
     * @param intervals A collection of intervals to merge.
     * @param mergingRule A descriptor for the type of merging to perform.
     * @return A sorted, merged version of the intervals passed in.
     */
    public static GenomeLocSortedSet sortAndMergeIntervals(GenomeLocParser parser, List<GenomeLoc> intervals, IntervalMergingRule mergingRule) {
        // Make a copy of the (potentially unmodifiable) list to be sorted
        intervals = new ArrayList<GenomeLoc>(intervals);
        // sort raw interval list
        Collections.sort(intervals);
        // now merge raw interval list
        intervals = mergeIntervalLocations(intervals, mergingRule);

        return GenomeLocSortedSet.createSetFromList(parser,intervals);
    }

    /**
     * computes whether the test interval list is equivalent to master.  To be equivalent, test must
     * contain GenomeLocs covering every base in master, exactly once.  Note that this algorithm
     * assumes that master genomelocs are all discontiguous (i.e., we don't have locs like 1-3 and 4-6 but
     * rather just 1-6).  In order to use this algorithm with contiguous genomelocs first merge them.  The algorithm
     * doesn't assume that test has discontinuous genomelocs.
     *
     * Returns a null string if there are no differences, otherwise returns a string describing the difference
     * (useful for UnitTests).  Assumes both lists are sorted
     *
     * @param masterArg sorted master genome locs
     * @param testArg sorted test genome locs
     * @return null string if there are no difference, otherwise a string describing the difference
     */
    public static String equateIntervals(List<GenomeLoc> masterArg, List<GenomeLoc> testArg) {
        LinkedList<GenomeLoc> master = new LinkedList<GenomeLoc>(masterArg);
        LinkedList<GenomeLoc> test = new LinkedList<GenomeLoc>(testArg);

        while ( ! master.isEmpty() ) { // there's still unchecked bases in master
            final GenomeLoc masterHead = master.pop();
            final GenomeLoc testHead = test.pop();

            if ( testHead.overlapsP(masterHead) ) {
                // remove the parts of test that overlap master, and push the remaining
                // parts onto master for further comparison.
                for ( final GenomeLoc masterPart : Utils.reverse(masterHead.subtract(testHead)) ) {
                    master.push(masterPart);
                }
            } else {
                // testHead is incompatible with masterHead, so we must have extra bases in testHead
                // that aren't in master
                return "Incompatible locs detected masterHead=" + masterHead + ", testHead=" + testHead;
            }
        }

        if ( test.isEmpty() ) // everything is equal
            return null; // no differences
        else
            return "Remaining elements found in test: first=" + test.peek();
    }

    /**
     * Converts a GenomeLoc to a picard interval.
     * @param loc The GenomeLoc.
     * @param locIndex The loc index for use in the file.
     * @return The picard interval.
     */
    private static htsjdk.samtools.util.Interval toInterval(GenomeLoc loc, int locIndex) {
        return new htsjdk.samtools.util.Interval(loc.getContig(), loc.getStart(), loc.getStop(), false, "interval_" + locIndex);
    }

    /**
     * merge a list of genome locs that may be overlapping, returning the list of unique genomic locations
     *
     * @param raw the unchecked genome loc list
     * @param rule the merging rule we're using
     *
     * @return the list of merged locations
     */
    public static List<GenomeLoc> mergeIntervalLocations(final List<GenomeLoc> raw, IntervalMergingRule rule) {
        if (raw.size() <= 1)
            return Collections.unmodifiableList(raw);
        else {
            ArrayList<GenomeLoc> merged = new ArrayList<GenomeLoc>();
            Iterator<GenomeLoc> it = raw.iterator();
            GenomeLoc prev = it.next();
            while (it.hasNext()) {
                GenomeLoc curr = it.next();
                if (prev.overlapsP(curr)) {
                    prev = prev.merge(curr);
                } else if (prev.contiguousP(curr) && (rule == null || rule == IntervalMergingRule.ALL)) {
                    prev = prev.merge(curr);
                } else {
                    merged.add(prev);
                    prev = curr;
                }
            }
            merged.add(prev);
            return Collections.unmodifiableList(merged);
        }
    }

    public static long intervalSize(final List<GenomeLoc> locs) {
        long size = 0;
        for ( final GenomeLoc loc : locs )
            size += loc.size();
        return size;
    }
}

