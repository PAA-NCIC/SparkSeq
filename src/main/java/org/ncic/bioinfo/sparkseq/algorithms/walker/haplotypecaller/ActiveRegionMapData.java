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

import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefMetaDataTracker;

/**
 * Data to use in the ActiveRegionWalker.map function produced by the NanoScheduler input iterator
 *
 * Author: wbc
 */
public final class ActiveRegionMapData {
    public ActiveRegion activeRegion;
    public RefMetaDataTracker tracker;
    public byte[] fullReferenceWithPadding;
    public byte[] refBases;

    public ActiveRegionMapData(ActiveRegion activeRegion, RefMetaDataTracker tracker,
                               byte[] fullReferenceWithPadding, byte[] refBases) {
        this.activeRegion = activeRegion;
        this.tracker = tracker;
        this.fullReferenceWithPadding = fullReferenceWithPadding;
        this.refBases = refBases;
    }

    public ActiveRegionMapData(ActiveRegion activeRegion, RefMetaDataTracker tracker,
                               RefContentProvider refContentProvider) {
        this.activeRegion = activeRegion;
        this.tracker = tracker;
        this.fullReferenceWithPadding = activeRegion.getActiveRegionReference(refContentProvider, HaplotypeCaller.REFERENCE_PADDING);
        this.refBases = activeRegion.getActiveRegionReference(refContentProvider);
    }
}