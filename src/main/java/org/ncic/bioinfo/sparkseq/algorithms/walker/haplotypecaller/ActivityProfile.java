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

import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLoc;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLocParser;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLocSortedSet;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

/**
 * Author: wbc
 */
public class ActivityProfile {
    protected final List<ActivityProfileState> stateList;
    protected final GenomeLocParser parser;
    protected final GenomeLocSortedSet restrictToIntervals;

    protected final int maxProbPropagationDistance;
    protected final double activeProbThreshold;

    protected GenomeLoc regionStartLoc = null;
    protected GenomeLoc regionStopLoc = null;

    /**
     * A cached value of the regionStartLoc contig length, to make calls to
     * getCurrentContigLength efficient
     */
    protected int contigLength = -1;

    /**
     * Create a empty ActivityProfile, restricting output to profiles overlapping intervals, if not null
     *
     * @param parser                     the parser we can use to create genome locs, cannot be null
     * @param maxProbPropagationDistance region probability propagation distance beyond it's maximum size
     * @param activeProbThreshold        threshold for the probability of a profile state being active
     * @param intervals                  only include states that are within these intervals, if not null
     */
    public ActivityProfile(final GenomeLocParser parser, final int maxProbPropagationDistance, final double activeProbThreshold, final GenomeLocSortedSet intervals) {
        if (parser == null) throw new IllegalArgumentException("parser cannot be null");

        this.parser = parser;
        this.stateList = new ArrayList<ActivityProfileState>();
        this.restrictToIntervals = intervals;
        this.maxProbPropagationDistance = maxProbPropagationDistance;
        this.activeProbThreshold = activeProbThreshold;
    }

    @Override
    public String toString() {
        return "ActivityProfile{" +
                "start=" + regionStartLoc +
                ", stop=" + regionStopLoc +
                '}';
    }

    /**
     * How far away can probability mass be moved around in this profile?
     * <p>
     * This distance puts an upper limit on how far, in bp, we will ever propagate probability max around
     * when adding a new ActivityProfileState.  For example, if the value of this function is
     * 10, and you are looking at a state at bp 5, and we know that no states beyond 5 + 10 will have
     * their probability propagated back to that state.
     *
     * @return a positive integer distance in bp
     */
    public int getMaxProbPropagationDistance() {
        return maxProbPropagationDistance;
    }

    /**
     * How many profile results are in this profile?
     *
     * @return the number of profile results
     */
    public int size() {
        return stateList.size();
    }

    /**
     * Is this profile empty?
     *
     * @return true if the profile is empty
     */
    public boolean isEmpty() {
        return stateList.isEmpty();
    }

    /**
     * Get the span of this activity profile, which is from the start of the first state to the stop of the last
     *
     * @return a potentially null GenomeLoc.  Will be null if this profile is empty
     */
    public GenomeLoc getSpan() {
        return isEmpty() ? null : regionStartLoc.endpointSpan(regionStopLoc);
    }

    public int getContigIndex() {
        return regionStartLoc.getContigIndex();
    }

    public int getStop() {
        return regionStopLoc.getStop();
    }

    /**
     * Get the list of active profile results in this object
     *
     * @return a non-null, ordered list of active profile results
     */
    protected List<ActivityProfileState> getStateList() {
        return stateList;
    }

    /**
     * Get the probabilities of the states as a single linear array of doubles
     *
     * @return a non-null array
     */
    protected double[] getProbabilitiesAsArray() {
        final double[] probs = new double[getStateList().size()];
        int i = 0;
        for (final ActivityProfileState state : getStateList())
            probs[i++] = state.isActiveProb;
        return probs;
    }

    /**
     * Helper function that gets the genome loc for a site offset from relativeLoc, protecting ourselves from
     * falling off the edge of the contig.
     *
     * @param relativeLoc the location offset is relative to
     * @param offset      the offset from relativeLoc where we'd like to create a GenomeLoc
     * @return a genome loc with relativeLoc.start + offset, if this is on the contig, null otherwise
     */
    protected GenomeLoc getLocForOffset(final GenomeLoc relativeLoc, final int offset) {
        final int start = relativeLoc.getStart() + offset;
        if (start < 0 || start > getCurrentContigLength()) {
            return null;
        } else {
            return parser.createGenomeLoc(regionStartLoc.getContig(), regionStartLoc.getContigIndex(), start, start);
        }
    }

    /**
     * Get the length of the current contig
     *
     * @return the length in bp
     */
    private int getCurrentContigLength() {
        return contigLength;
    }

    // --------------------------------------------------------------------------------
    //
    // routines to add states to a profile
    //
    // --------------------------------------------------------------------------------

    /**
     * Add the next ActivityProfileState to this profile.
     * <p>
     * Must be contiguous with the previously added result, or an IllegalArgumentException will be thrown
     *
     * @param state a well-formed ActivityProfileState result to incorporate into this profile
     */
    public void add(final ActivityProfileState state) {
        final GenomeLoc loc = state.getLoc();

        if (regionStartLoc == null) {
            regionStartLoc = loc;
            regionStopLoc = loc;
            contigLength = parser.getContigInfo(regionStartLoc.getContig()).getSequenceLength();
        } else {
            if (regionStopLoc.getStart() != loc.getStart() - 1)
                throw new IllegalArgumentException("Bad add call to ActivityProfile: loc " + loc + " not immediately after last loc " + regionStopLoc);
            regionStopLoc = loc;
        }

        final Collection<ActivityProfileState> processedStates = processState(state);
        for (final ActivityProfileState processedState : processedStates) {
            incorporateSingleState(processedState);
        }
    }

    /**
     * Incorporate a single activity profile state into the current list of states
     * <p>
     * If state's position occurs immediately after the last position in this profile, then
     * the state is appended to the state list.  If it's within the existing states list,
     * the prob of stateToAdd is added to its corresponding state in the list.  If the
     * position would be before the start of this profile, stateToAdd is simply ignored.
     *
     * @param stateToAdd the state we want to add to the states list
     */
    private void incorporateSingleState(final ActivityProfileState stateToAdd) {
        final int position = stateToAdd.getOffset(regionStartLoc);

        if (position > size())
            // should we allow this?  probably not
            throw new IllegalArgumentException("Must add state contiguous to existing states: adding " + stateToAdd);

        if (position >= 0) {
            // ignore states starting before this region's start
            if (position < size()) {
                stateList.get(position).isActiveProb += stateToAdd.isActiveProb;
            } else {
                if (position != size())
                    throw new IllegalStateException("position == size but it wasn't");
                stateList.add(stateToAdd);
            }
        }
    }

    /**
     * Process justAddedState, returning a collection of derived states that actually be added to the stateList
     * <p>
     * The purpose of this function is to transform justAddedStates, if needed, into a series of atomic states
     * that we actually want to track.  For example, if state is for soft clips, we transform that single
     * state into a list of states that surround the state up to the distance of the soft clip.
     * <p>
     * Can be overridden by subclasses to transform states in any way
     * <p>
     * There's no particular contract for the output states, except that they can never refer to states
     * beyond the current end of the stateList unless the explicitly include preceding states before
     * the reference.  So for example if the current state list is [1, 2, 3] this function could return
     * [1,2,3,4,5] but not [1,2,3,5].
     *
     * @param justAddedState the state our client provided to use to add to the list
     * @return a list of derived states that should actually be added to this profile's state list
     */
    protected Collection<ActivityProfileState> processState(final ActivityProfileState justAddedState) {
        if (justAddedState.resultState.equals(ActivityProfileState.Type.HIGH_QUALITY_SOFT_CLIPS)) {
            // special code to deal with the problem that high quality soft clipped bases aren't added to pileups
            final List<ActivityProfileState> states = new LinkedList<ActivityProfileState>();
            // add no more than the max prob propagation distance num HQ clips
            final int numHQClips = Math.min(justAddedState.resultValue.intValue(), getMaxProbPropagationDistance());
            for (int jjj = -numHQClips; jjj <= numHQClips; jjj++) {
                final GenomeLoc loc = getLocForOffset(justAddedState.getLoc(), jjj);
                if (loc != null)
                    states.add(new ActivityProfileState(loc, justAddedState.isActiveProb));
            }

            return states;
        } else {
            return Collections.singletonList(justAddedState);
        }
    }

    // --------------------------------------------------------------------------------
    //
    // routines to get active regions from the profile
    //
    // --------------------------------------------------------------------------------

    /**
     * Get the next completed active regions from this profile, and remove all states supporting them from this profile
     * <p>
     * Takes the current profile and finds all of the active / inactive from the start of the profile that are
     * ready.  By ready we mean unable to have their probability modified any longer by future additions to the
     * profile.  The regions that are popped off the profile take their states with them, so the start of this
     * profile will always be after the end of the last region returned here.
     * <p>
     * The regions are returned sorted by genomic position.
     * <p>
     * This function may not return anything in the list, if no regions are ready
     * <p>
     * No returned region will be larger than maxRegionSize.
     *
     * @param activeRegionExtension the extension value to provide to the constructed regions
     * @param minRegionSize         the minimum region size, in the case where we have to cut up regions that are too large
     * @param maxRegionSize         the maximize size of the returned region
     * @param forceConversion       if true, we'll return a region whose end isn't sufficiently far from the end of the
     *                              stateList.  Used to close out the active region when we've hit some kind of end (such
     *                              as the end of the contig)
     * @return a non-null list of active regions
     */
    public List<ActiveRegion> popReadyActiveRegions(final int activeRegionExtension, final int minRegionSize, final int maxRegionSize, final boolean forceConversion) {
        if (activeRegionExtension < 0)
            throw new IllegalArgumentException("activeRegionExtension must be >= 0 but got " + activeRegionExtension);
        if (minRegionSize < 1)
            throw new IllegalArgumentException("minRegionSize must be >= 1 but got " + minRegionSize);
        if (maxRegionSize < 1)
            throw new IllegalArgumentException("maxRegionSize must be >= 1 but got " + maxRegionSize);

        final LinkedList<ActiveRegion> regions = new LinkedList<ActiveRegion>();

        curStateIdx = 0;  // 遍历stateList的index
        while (true) {
            final ActiveRegion nextRegion = popNextReadyActiveRegion(activeRegionExtension, minRegionSize, maxRegionSize, forceConversion);
            if (nextRegion == null)
                return regions;
            else {
                if (restrictToIntervals == null)
                    regions.add(nextRegion);
                else
                    regions.addAll(nextRegion.splitAndTrimToIntervals(restrictToIntervals));
            }
        }
    }

    private int curStateIdx = 0;    // 在pop active region的时候，用于遍历stateList的下标

    /**
     * Helper function for popReadyActiveRegions that pops the first ready region off the front of this profile
     * <p>
     * If a region is returned, modifies the state of this profile so that states used to make the region are
     * no longer part of the profile.  Associated information (like the region start position) of this profile
     * are also updated.
     *
     * @param activeRegionExtension the extension value to provide to the constructed regions
     * @param minRegionSize         the minimum region size, in the case where we have to cut up regions that are too large
     * @param maxRegionSize         the maximize size of the returned region
     * @param forceConversion       if true, we'll return a region whose end isn't sufficiently far from the end of the
     *                              stateList.  Used to close out the active region when we've hit some kind of end (such
     *                              as the end of the contig)
     * @return a fully formed active region, or null if none can be made
     */
    private ActiveRegion popNextReadyActiveRegion(final int activeRegionExtension,
                                                  final int minRegionSize, final int maxRegionSize, final boolean forceConversion) {
        if (curStateIdx >= stateList.size())
            return null;

        final ActivityProfileState first = stateList.get(curStateIdx);
        final boolean isActiveRegion = first.isActiveProb > activeProbThreshold;
        //得到的active region是[curStateIdx， offsetOfNextRegionEnd]
        final int offsetOfNextRegionEnd = findEndOfRegion(isActiveRegion, curStateIdx, minRegionSize, maxRegionSize, forceConversion);

        /*if (offsetOfNextRegionEnd == -1)
            // couldn't find a valid ending offset, so we return null
            return null;*/ //我们不会提前判断，也就没有-1

        int stateCount = offsetOfNextRegionEnd - curStateIdx + 1;
        List<ActivityProfileState> supportingStates = new ArrayList<>(stateCount);
        for(int i = curStateIdx; i <= offsetOfNextRegionEnd; i ++) {
            supportingStates.add(stateList.get(i));
        }

        curStateIdx = offsetOfNextRegionEnd + 1;
        /*// update the start and stop locations as necessary
        if (curStateIdx >= stateList.size()) {
            regionStartLoc = regionStopLoc = null;
        } else {
            regionStartLoc = stateList.get(0).getLoc();
        }*/
        final GenomeLoc regionLoc = parser.createGenomeLoc(first.getLoc().getContig(),
                first.getLoc().getStart(), stateList.get(offsetOfNextRegionEnd).getLoc().getStart());
        return new ActiveRegion(regionLoc, supportingStates, isActiveRegion, parser, activeRegionExtension);
    }

    /**
     * Find the end of the current region, returning the index into the element isActive element, or -1 if the region isn't done
     * <p>
     * The current region is defined from the start of the stateList, looking for elements that have the same isActiveRegion
     * flag (i.e., if isActiveRegion is true we are looking for states with isActiveProb > threshold, or alternatively
     * for states < threshold).  The maximize size of the returned region is maxRegionSize.  If forceConversion is
     * true, then we'll return the region end even if this isn't safely beyond the max prob propagation distance.
     * <p>
     * Note that if isActiveRegion is true, and we can construct a active region > maxRegionSize in bp, we
     * find the further local minimum within that max region, and cut the region there, under the constraint
     * that the resulting region must be at least minRegionSize in bp.
     *
     * @param isActiveRegion  is the region we're looking for an active region or inactive region?
     * @param startSearchIdx  开始查找的state的位置下标
     * @param minRegionSize   the minimum region size, in the case where we have to cut up regions that are too large
     * @param maxRegionSize   the maximize size of the returned region
     * @param forceConversion if true, we'll return a region whose end isn't sufficiently far from the end of the
     *                        stateList.  Used to close out the active region when we've hit some kind of end (such
     *                        as the end of the contig)
     * @return the index into stateList of the last element of this region, or -1 if it cannot be found
     */
    private int findEndOfRegion(final boolean isActiveRegion, final int startSearchIdx, final int minRegionSize, final int maxRegionSize, final boolean forceConversion) {
        /*if (!forceConversion && stateList.size() < maxRegionSize + getMaxProbPropagationDistance()) {
            // we really haven't finalized at the probability mass that might affect our decision, so keep
            // waiting until we do before we try to make any decisions
            return -1;
        }*/ // 因为我们在这里的实现是把所有state都算出来再一起pop，所以不会有“等等再决定”这种事情发生，因此这种判断不会发生 by wangbc

        // 更换了end of active region的定义，这里的active region的范围是从start search index到end of active region(左闭右开区间)
        int endOfActiveRegion = findFirstActivityBoundary(startSearchIdx, isActiveRegion, maxRegionSize);
        int regionSize = endOfActiveRegion - startSearchIdx;

        if (isActiveRegion && regionSize == maxRegionSize)
            // we've run to the end of the region, let's find a good place to cut
            endOfActiveRegion = findBestCutSite(startSearchIdx, endOfActiveRegion, minRegionSize);

        // we're one past the end, so i must be decremented
        return endOfActiveRegion - 1;
    }

    /**
     * Find the the local minimum within 0 - endOfActiveRegion where we should divide region
     * <p>
     * This algorithm finds the global minimum probability state within the region [minRegionSize, endOfActiveRegion)
     * (exclusive of endOfActiveRegion), and returns the state index of that state.
     * that it
     *
     * @param startSearchIdx    开始查找的state的位置下标
     * @param endOfActiveRegion the last state of the current active region (exclusive)
     * @param minRegionSize     the minimum of the left-most region, after cutting
     * @return the index of state after the cut site (just like endOfActiveRegion)
     */
    private int findBestCutSite(final int startSearchIdx, final int endOfActiveRegion, final int minRegionSize) {
        int minI = endOfActiveRegion - 1;
        double minP = Double.MAX_VALUE;

        // 从后往前找，找到最小的那个极小值
        int minEndIdxOfRegion = startSearchIdx + minRegionSize - 1;
        // 如果不是在stateList边上，就不要尝试边上那个点，因为不可能是极小值（为了省掉isMinimum里面那个判断）
        int startValue = (minI == stateList.size() - 1) ? minI - 1 : minI;
        for (int i = startValue; i >= minEndIdxOfRegion; i--) {
            double cur = getProb(i);
            if (cur < minP && isMinimum(i)) {
                minP = cur;
                minI = i;
            }
        }

        return minI + 1;
    }

    /**
     * Find the first index into the state list where the state is considered ! isActiveRegion
     * <p>
     * Note that each state has a probability of being active, and this function thresholds that
     * value on activeProbThreshold, coloring each state as active or inactive.  Finds the
     * largest contiguous stretch of states starting at the first state (index 0) with the same isActive
     * state as isActiveRegion.  If the entire state list has the same isActive value, then returns
     * maxRegionSize
     *
     * @param startSearchIdx 开始查找的state的位置下标
     * @param isActiveRegion are we looking for a stretch of active states, or inactive ones?
     * @param maxRegionSize  don't look for a boundary that would yield a region of size > maxRegionSize
     * @return the index of the first state in the state list with isActive value != isActiveRegion, 如果超过了maxRegionSize，会提前截断
     */
    private int findFirstActivityBoundary(final int startSearchIdx, final boolean isActiveRegion, final int maxRegionSize) {
        final int nStates = stateList.size();
        int endOfActiveRegion = startSearchIdx;

        int regionSize = 0;
        while (endOfActiveRegion < nStates && regionSize < maxRegionSize) {
            if (stateList.get(endOfActiveRegion).isActiveProb > activeProbThreshold != isActiveRegion) {
                break;
            }
            endOfActiveRegion++;
            regionSize++;
        }

        return endOfActiveRegion;
    }

    /**
     * Helper function to get the probability of the state at offset index
     *
     * @param index a valid offset into the state list
     * @return the isActiveProb of the state at index
     */
    private double getProb(final int index) {
        return stateList.get(index).isActiveProb;
    }

    /**
     * Is the probability at index in a local minimum?
     * <p>
     * Checks that the probability at index is <= both the probabilities to either side.
     * Returns false if index is at the end or the start of the state list.
     *
     * @param index the index of the state we want to test
     * @return true if prob at state is a minimum, false otherwise
     */
    private boolean isMinimum(final int index) {
        /*if (index == stateList.size() - 1)
            // we cannot be at a minimum if the current position is the last in the state list
            return false;
        else if (index < 1)
            // we cannot be at a minimum if the current position is the first or second
            return false;*/ //我们不会让这种调用发生 by wangbc
        final double indexP = getProb(index);
        return indexP <= getProb(index + 1) && indexP < getProb(index - 1);
    }
}