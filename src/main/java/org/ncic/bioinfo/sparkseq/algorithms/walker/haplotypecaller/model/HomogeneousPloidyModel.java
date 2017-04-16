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
package org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.model;

import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.PloidyModel;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.SampleList;

/**
 * Author: wbc
 */
public class HomogeneousPloidyModel implements PloidyModel, SampleList {

    private SampleList sampleList;

    private final int ploidy;

    /**
     * Constructs a homogeneous ploidy model given the sample list and ploidy.
     *
     * @param samples the sample list.
     * @param ploidy the common ploidy for all samples in {@code samples}.
     *
     * @throws IllegalArgumentException if {@code samples} is {@code null},
     *    or ploidy is 0 or less.
     */
    public HomogeneousPloidyModel(final SampleList samples, final int ploidy) {
        if (ploidy <= 0)
            throw new IllegalArgumentException("does not support negative ploidy");
        this.ploidy = ploidy;

        sampleList = samples;
    }

    @Override
    public int sampleCount() {
        return sampleList.sampleCount();
    }

    @Override
    public String sampleAt(final int index) {
        return sampleList.sampleAt(index);
    }

    @Override
    public int sampleIndex(final String sample) {
        return sampleList.sampleIndex(sample);
    }

    @Override
    public int samplePloidy(final int sampleIndex) {
        checkSampleIndex(sampleIndex);
        return ploidy;
    }

    private void checkSampleIndex(final int sampleIndex) {
        if (sampleIndex < 0)
            throw new IllegalArgumentException("the sample index cannot be negative: " + sampleIndex);
        if (sampleIndex >= sampleList.sampleCount())
            throw new IllegalArgumentException("the sample index is equal or larger than the sample count: " + sampleIndex + " >= " + sampleList.sampleCount());
    }

    @Override
    public boolean isHomogeneous() {
        return true;
    }

    @Override
    public int totalPloidy() {
        return ploidy * sampleList.sampleCount();
    }

}
