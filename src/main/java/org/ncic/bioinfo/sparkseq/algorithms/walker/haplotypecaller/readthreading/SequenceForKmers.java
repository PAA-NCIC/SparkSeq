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
package org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.readthreading;

/**
 * Keeps track of the information needed to add a sequence to the read threading assembly graph
 * <p>
 * Author: wbc
 */
final class SequenceForKmers {
    final String name;
    final byte[] sequence;
    final int start, stop;
    final int count;
    final boolean isRef;

    /**
     * Create a new sequence for creating kmers
     */
    SequenceForKmers(final String name, byte[] sequence, int start, int stop, int count, boolean ref) {
        if (start < 0) throw new IllegalArgumentException("Invalid start " + start);
        if (stop < start) throw new IllegalArgumentException("Invalid stop " + stop);
        if (sequence == null) throw new IllegalArgumentException("Sequence is null ");
        if (count < 1) throw new IllegalArgumentException("Invalid count " + count);

        this.name = name;
        this.sequence = sequence;
        this.start = start;
        this.stop = stop;
        this.count = count;
        this.isRef = ref;
    }
}
