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
package org.ncic.bioinfo.sparkseq.algorithms.data.sam;

import org.ncic.bioinfo.sparkseq.exceptions.ReviewedGATKException;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLoc;
import org.ncic.bioinfo.sparkseq.algorithms.utils.HasGenomeLocation;

/**
 * Useful class for forwarding on locusContext data from this iterator
 * Author: wbc
 */
public class AlignmentContext implements HasGenomeLocation {
    protected GenomeLoc loc = null;
    protected ReadBackedPileup basePileup = null;
    protected boolean hasPileupBeenDownsampled;

    /**
     * The number of bases we've skipped over in the reference since the last map invocation.
     * Only filled in by RodTraversals right now.  By default, nothing is being skipped, so skippedBases == 0.
     */
    private long skippedBases = 0;

    public AlignmentContext(GenomeLoc loc, ReadBackedPileup basePileup) {
        this(loc, basePileup, 0, false);
    }

    public AlignmentContext(GenomeLoc loc, ReadBackedPileup basePileup, long skippedBases, boolean hasPileupBeenDownsampled) {
        if (loc == null)
            throw new ReviewedGATKException("BUG: GenomeLoc in Alignment context is null");
        if (basePileup == null)
            throw new ReviewedGATKException("BUG: ReadBackedPileup in Alignment context is null");
        if (skippedBases < 0)
            throw new ReviewedGATKException("BUG: skippedBases is -1 in Alignment context");

        this.loc = loc;
        this.basePileup = basePileup;
        this.skippedBases = skippedBases;
        this.hasPileupBeenDownsampled = hasPileupBeenDownsampled;
    }

    /**
     * Returns base pileup over the current genomic location. May return null if this context keeps only
     * extended event (indel) pileup.
     *
     * @return
     */
    public ReadBackedPileup getBasePileup() {
        return basePileup;
    }

    /**
     * Returns true if any reads have been filtered out of the pileup due to excess DoC.
     *
     * @return True if reads have been filtered out.  False otherwise.
     */
    public boolean hasPileupBeenDownsampled() {
        return hasPileupBeenDownsampled;
    }

    /**
     * How many reads cover this locus?
     *
     * @return
     */
    public int size() {
        return basePileup.getNumberOfElements();
    }

    public String getContig() {
        return getLocation().getContig();
    }

    public long getPosition() {
        return getLocation().getStart();
    }

    public GenomeLoc getLocation() {
        return loc;
    }

}