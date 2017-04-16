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
package org.ncic.bioinfo.sparkseq.algorithms.walker.indelrealigner;

import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLoc;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLocParser;
import org.ncic.bioinfo.sparkseq.algorithms.utils.HasGenomeLocation;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.ReferenceContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.GATKSAMRecord;

import java.util.ArrayList;
import java.util.List;

/**
 * Author: wbc
 */
class ReadBin implements HasGenomeLocation {

    private final ArrayList<GATKSAMRecord> reads = new ArrayList<GATKSAMRecord>();
    private byte[] reference = null;
    private GenomeLoc loc = null;
    private final GenomeLocParser parser;
    private final int referencePadding;

    public ReadBin(final GenomeLocParser parser, final int referencePadding) {
        this.parser = parser;
        this.referencePadding = referencePadding;
    }

    // Return false if we can't process this read bin because the reads are not correctly overlapping.
    // This can happen if e.g. there's a large known indel with no overlapping reads.
    public void add(GATKSAMRecord read) {

        final int readStart = read.getSoftStart();
        final int readStop = read.getSoftEnd();
        if ( loc == null )
            loc = parser.createGenomeLoc(read.getReferenceName(), readStart, Math.max(readStop, readStart)); // in case it's all an insertion
        else if ( readStop > loc.getStop() )
            loc = parser.createGenomeLoc(loc.getContig(), loc.getStart(), readStop);

        reads.add(read);
    }

    public List<GATKSAMRecord> getReads() {
        return reads;
    }

    public byte[] getReference(RefContentProvider refContentProvider) {
        // set up the reference if we haven't done so yet
        if ( reference == null ) {
            ReferenceContext context = refContentProvider.getReferenceContext(loc, referencePadding);
            reference = context.getBases();
            loc = context.getLocus();
        }

        return reference;
    }

    public GenomeLoc getLocation() {
        return loc;
    }

    public int size() {
        return reads.size();
    }

    public void clear() {
        reads.clear();
        reference = null;
        loc = null;
    }

}