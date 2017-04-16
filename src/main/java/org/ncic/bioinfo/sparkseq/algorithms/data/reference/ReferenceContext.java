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
package org.ncic.bioinfo.sparkseq.algorithms.data.reference;

import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLoc;

/**
 * Author: wbc
 */
public class ReferenceContext {

    final int contigId;
    final private GenomeLoc locus;
    final private byte[] bases;

    /**
     * The window of reference information around the current locus.
     */
    final private GenomeLoc window;

    public ReferenceContext(GenomeLoc locus, int contigId, byte[] bases) {
        this.locus = locus;
        this.window = locus;
        this.contigId = contigId;
        this.bases = bases;
    }

    public ReferenceContext(GenomeLoc locus, GenomeLoc window, int contigId, byte[] bases) {
        this.locus = locus;
        this.window = window;
        this.contigId = contigId;
        this.bases = bases;
    }

    /**
     * Contig id of this reference
     *
     * @return contig id
     */
    public int getContigId() {
        return contigId;
    }

    /**
     * The locus currently being examined.
     *
     * @return The current locus.
     */
    public GenomeLoc getLocus() {
        return locus;
    }

    public GenomeLoc getWindow() {
        return window;
    }

    /**
     * Get the base at the given locus.
     *
     * @return The base at the given locus from the reference.
     */
    public byte getBase() {
        return bases[0];
    }

    /**
     * All the bases in the window currently being examined.
     *
     * @return All bases available.  If the window is of size [0,0], the array will
     * contain only the base at the given locus.
     */
    public byte[] getBases() {
        return bases;
    }

    public byte[] getForwardBases() {
        final byte[] bases = getBases();
        final int mid = locus.getStart() - window.getStart();
        // todo -- warning of performance problem, especially if this is called over and over
        return new String(bases).substring(mid).getBytes();
    }

}

