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

import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.graphs.Kmer;

import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * Author: wbc
 */
public class KMerCounter {
    //private final static Logger logger = Logger.getLogger(KMerCounter.class);

    /**
     * A map of for each kmer to its num occurrences in addKmers
     */
    private final Map<Kmer, CountedKmer> countsByKMer = new HashMap<Kmer, CountedKmer>();
    private final int kmerLength;

    /**
     * Create a new kmer counter
     *
     * @param kmerLength the length of kmers we'll be counting to error correct, must be >= 1
     */
    public KMerCounter(final int kmerLength) {
        if ( kmerLength < 1 ) throw new IllegalArgumentException("kmerLength must be > 0 but got " + kmerLength);
        this.kmerLength = kmerLength;
    }

    /**
     * Get the count of kmer in this kmer counter
     * @param kmer a non-null counter to get
     * @return a positive integer
     */
    public int getKmerCount(final Kmer kmer) {
        if ( kmer == null ) throw new IllegalArgumentException("kmer cannot be null");
        final CountedKmer counted = countsByKMer.get(kmer);
        return counted == null ? 0 : counted.count;
    }

    /**
     * Get an unordered collection of the counted kmers in this counter
     * @return a non-null collection
     */
    public Collection<CountedKmer> getCountedKmers() {
        return countsByKMer.values();
    }

    /**
     * Get kmers that have minCount or greater in this counter
     * @param minCount only return kmers with count >= this value
     * @return a non-null collection of kmers
     */
    public Collection<Kmer> getKmersWithCountsAtLeast(final int minCount) {
        final List<Kmer> result = new LinkedList<Kmer>();
        for ( final CountedKmer countedKmer : getCountedKmers() ) {
            if ( countedKmer.count >= minCount )
                result.add(countedKmer.kmer);
        }
        return result;
    }

    /**
     * Remove all current counts, resetting the counter to an empty state
     */
    public void clear() {
        countsByKMer.clear();
    }

    /**
     * Add a kmer that occurred kmerCount times
     *
     * @param kmer a kmer
     * @param kmerCount the number of occurrences
     */
    public void addKmer(final Kmer kmer, final int kmerCount) {
        if ( kmer.length() != kmerLength ) throw new IllegalArgumentException("bad kmer length " + kmer + " expected size " + kmerLength);
        if ( kmerCount < 0 ) throw new IllegalArgumentException("bad kmerCount " + kmerCount);

        CountedKmer countFromMap = countsByKMer.get(kmer);
        if ( countFromMap == null ) {
            countFromMap = new CountedKmer(kmer);
            countsByKMer.put(kmer, countFromMap);
        }
        countFromMap.count += kmerCount;
    }

    @Override
    public String toString() {
        final StringBuilder b = new StringBuilder("KMerCounter{");
        b.append("counting ").append(countsByKMer.size()).append(" distinct kmers");
        b.append("\n}");
        return b.toString();
    }

    public static class CountedKmer implements Comparable<CountedKmer> {
        public final Kmer kmer;
        public int count = 0;

        private CountedKmer(final Kmer kmer) {
            this.kmer = kmer;
        }

        public Kmer getKmer() {
            return kmer;
        }

        public int getCount() {
            return count;
        }

        @Override
        public String toString() {
            return "CountedKmer{" +
                    "kmer='" + kmer + '\'' +
                    ", count=" + count +
                    '}';
        }

        @Override
        public int compareTo(CountedKmer o) {
            return o.count - count;
        }
    }

    // -------------------------------------------------------------------------------------
    // Protected methods for testing purposes only
    // -------------------------------------------------------------------------------------

    /**
     * For testing purposes only
     */
    protected void addKmer(final String rawKmer, final int kmerCount) {
        addKmer(new Kmer(rawKmer), kmerCount);
    }

    /**
     * For testing purposes
     *
     * @param kmers
     */
    protected void addKmers(final String ... kmers) {
        for ( final String kmer : kmers )
            addKmer(kmer, 1);
    }
}
