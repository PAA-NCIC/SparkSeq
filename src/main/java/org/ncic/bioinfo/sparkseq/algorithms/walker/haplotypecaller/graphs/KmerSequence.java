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
package org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.graphs;

import htsjdk.samtools.SAMRecord;
import org.ncic.bioinfo.sparkseq.algorithms.utils.haplotype.Haplotype;

import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;

/**
 * Author: wbc
 */
public class KmerSequence implements List<Kmer> {
    private final byte[] sequence;
    private final int start;
    private final int size;
    private final int kmerSize;
    private final int rawLength;

    /**
     * Creates a kmer sequence from a read's sequence.
     *
     * @param read the read to represent as a sequence of kmers.
     * @param kmerSize the kmer size.
     */
    public KmerSequence(final SAMRecord read, final int kmerSize) {
        this(read.getReadBases(), kmerSize);
    }

    /**
     * Creates a kmer sequence from a haplotype's sequence.
     *
     * @param hap the haplotype to represent as a sequence of kmers.
     * @param kmerSize the kmer size.
     */
    @SuppressWarnings("unused")
    public KmerSequence(final Haplotype hap, final int kmerSize) {
        this(hap.getBases(), kmerSize);
    }

    /**
     * Creates a kmer sequence out of a byte sequence.
     *
     * @param sequence the byte array to represent as a kmer sequence.
     * @param kmerSize the kmer size.
     */
    public KmerSequence(final byte[] sequence, final int kmerSize) {
        this(sequence,0,Math.max(0,sequence.length - kmerSize + 1),kmerSize, sequence.length);
    }

    /**
     * Creates a kmer sequence out of a range of a byte array
     *
     * @param sequence the input array.
     * @param start inclusive first position of the array that maps to the first position in the first kmer.
     * @param size number kmers in the output.
     * @param kmerSize kmer length in bases.
     * @param rawLength the of the range in bases.
     */
    protected KmerSequence(final byte[] sequence, final int start, final int size, final int kmerSize, final int rawLength) {
        if (sequence == null) {
            throw new IllegalArgumentException("start must be 0 or greater");
        }
        if (rawLength > sequence.length - start) {
            throw new IllegalArgumentException("the raw sequence length goes beyond the array capacity");
        }
        if (size < 0) {
            throw new IllegalArgumentException("the length cannot be negative");
        }
        if (start < 0) {
            throw new IllegalArgumentException("start must be 0 or greater");
        }
        if (size > 0 && size + kmerSize - 1 > rawLength) {
            throw new IllegalArgumentException(
                    String.format("the kmerSize (%d) + size (%d) - 1 cannot be larger than rawLength (%d)",kmerSize,size,rawLength) );
        }
        this.sequence = sequence;
        this.start = start;
        this.size = size;
        this.kmerSize = kmerSize;
        this.rawLength = rawLength;
    }

    public int kmerSize() {
        return kmerSize;
    }

    public KmerSequence subsequence(final int from, final int to) {
        if (from < 0 || from > to) {
            throw new IllegalArgumentException();
        }
        if (to > size) {
            throw new IllegalArgumentException();
        }
        return new KmerSequence(sequence,this.start + from,to - from,kmerSize,rawLength - from - (size - to));
    }


    @Override
    public int size() {
        return size;
    }

    @Override
    public boolean isEmpty() {
        return size == 0;
    }

    @Override
    public boolean contains(final Object o) {
        if (o instanceof Kmer) {
            if (o instanceof MyKmer) {
                final MyKmer k = (MyKmer) o;
                if (k.bases == sequence && k.start >= start && k.length == kmerSize && k.start < start + size) {
                    return true;
                }
            }
            final Kmer k = (Kmer) o;
            if (k.length != kmerSize) {
                return false;
            }
            for (int i = 0; i < size; i++) {
                int j;
                for (j = 0; j < kmerSize; j++) {
                    if (sequence[start + i + j] != k.bases[k.start + j]) {
                        break;
                    }
                }
                if (j == kmerSize) {
                    return true;
                }
            }
            return false;
        } else {
            return false;
        }
    }

    @Override
    public Iterator<Kmer> iterator() {
        return new Iterator<Kmer>() {

            private int offset = 0;

            @Override
            public boolean hasNext() {
                return offset < size;
            }

            @Override
            public Kmer next() {
                return new Kmer(sequence,start + offset,kmerSize);
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    @Override
    public Object[] toArray() {
        return toArray(new Kmer[size()]);
    }

    @Override
    @SuppressWarnings("unchecked")
    public <T> T[] toArray(final T[] a) {
        if (a == null) {
            throw new IllegalArgumentException();
        } else if (!a.getClass().getComponentType().isAssignableFrom(Kmer.class)) {
            throw new IllegalArgumentException();
        } else {
            T[] result;
            if (a.length < size) {
                result = (T[]) Array.newInstance(a.getClass().getComponentType(), size);
            } else {
                result = a;
            }
            for (int i = 0; i < size; i++) {
                result[i] = (T) new Kmer(sequence,start + i,kmerSize);
            }
            return result;
        }
    }

    @Override
    public boolean add(final Kmer kmer) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean remove(final Object o) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean containsAll(final Collection<?> c) {
        for (final Object o : c)
            if (!contains(o))
                return false;
        return true;
    }

    @Override
    public boolean addAll(final Collection<? extends Kmer> c) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean addAll(final int index, final Collection<? extends Kmer> c) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean removeAll(final Collection<?> c) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean retainAll(final Collection<?> c) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void clear() {
        throw new UnsupportedOperationException();
    }

    @Override
    public Kmer get(final int index) {
        if (index < 0 || index >= size) {
            throw new IllegalArgumentException();
        }
        return new Kmer(sequence,start + index,kmerSize);
    }

    @Override
    public Kmer set(final int index, final Kmer element) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void add(final int index, final Kmer element) {
        throw new UnsupportedOperationException();
    }

    @Override
    public Kmer remove(final int index) {
        throw new UnsupportedOperationException();
    }

    @Override
    public int indexOf(final Object o) {
        if (o instanceof Kmer) {
            final Kmer k = (Kmer) o;
            if (k.length != kmerSize) {
                return -1;
            }
            for (int i = 0; i < size; i++) {
                int j;
                for (j = 0; j < kmerSize; j++) {
                    if (sequence[start + i + j] != k.bases[k.start + j]) {
                        break;
                    }
                }
                if (j == kmerSize) {
                    return i;
                }
            }
            return -1;
        } else {
            return -1;
        }
    }

    @Override
    public int lastIndexOf(final Object o) {
        if (o instanceof Kmer) {
            final Kmer k = (Kmer) o;
            if (k.length != kmerSize) {
                return -1;
            }
            for (int i = size - 1; i >= 0; i--) {
                int j;
                for (j = kmerSize - 1; j >= 0; j--) {
                    if (sequence[start + i + j] != k.bases[k.start + j]) {
                        break;
                    }
                }
                if (j == 0) {
                    return i;
                }
            }
            return -1;
        } else {
            return -1;
        }
    }

    @Override
    public ListIterator<Kmer> listIterator() {
        return new MyListIterator(0);
    }

    @Override
    public ListIterator<Kmer> listIterator(final int index) {
        return new MyListIterator(index);
    }

    @Override
    public List<Kmer> subList(final int fromIndex, final int toIndex) {
        return subsequence(fromIndex,toIndex);
    }

    /**
     * Returns the byte array representation of the kmer sequence.
     * @return never {@code null}.
     */
    public byte[] getBytes() {
        if (start == 0 && rawLength == sequence.length)
            return sequence;
        else
            return Arrays.copyOfRange(sequence, start, rawLength + start);
    }

    /**
     * Internal class that implements the {@link Kmer} more efficiently
     * making reference to the sequence's own byte array.
     */
    protected class MyKmer extends Kmer {

        /**
         * Create a new instance give the offset in the byte array.
         * @param start the start base offset for the kmer.
         */
        public MyKmer(final int start) {
            super(sequence,start,kmerSize);
        }
    }

    /**
     * Iterator implementation of Kmer elements.
     */
    private class MyListIterator implements ListIterator<Kmer> {

        private int i = 0;

        /**
         * Creates a iterator at certain offset in the sequence.
         * @param idx the start position or kmer offset.
         */
        private MyListIterator(final int idx) {
            i = idx;
        }

        @Override
        public boolean hasNext() {
            return i < size;
        }

        @Override
        public Kmer next() {
            return new Kmer(sequence,start + i++,kmerSize);
        }

        @Override
        public boolean hasPrevious() {
            return i > 0;
        }

        @Override
        public Kmer previous() {
            return new Kmer(sequence,start + --i,kmerSize);
        }

        @Override
        public int nextIndex() {
            return i;
        }

        @Override
        public int previousIndex() {
            return i - 1;
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException();
        }

        @Override
        public void set(final Kmer kmer) {
            throw new UnsupportedOperationException();
        }

        @Override
        public void add(final Kmer kmer) {
            throw new UnsupportedOperationException();
        }

    }

}