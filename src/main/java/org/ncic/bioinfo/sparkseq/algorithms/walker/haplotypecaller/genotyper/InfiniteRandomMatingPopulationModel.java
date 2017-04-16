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
package org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.ncic.bioinfo.sparkseq.algorithms.utils.AlleleListUtils;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.AlleleList;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.AlleleListPermutation;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.GenotypeLikelihoodCalculator;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.GenotypeLikelihoodCalculators;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.GenotypingData;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.GenotypingLikelihoods;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.PloidyModel;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.readlikelihood.ReadLikelihoods;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.model.GenotypingModel;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Author: wbc
 */
public class InfiniteRandomMatingPopulationModel implements GenotypingModel {

    private final int cachePloidyCapacity;
    private final int cacheAlleleCountCapacity;
    private ThreadLocal<GenotypeLikelihoodCalculator[][]> likelihoodCalculators;

    /**
     * Create a new infinite model instance.
     */
    public InfiniteRandomMatingPopulationModel() {
        this(10,50);
    }

    public InfiniteRandomMatingPopulationModel(final int calculatorCachePloidyCapacity, final int calculatorCacheAlleleCapacity) {
        cachePloidyCapacity = calculatorCachePloidyCapacity;
        cacheAlleleCountCapacity = calculatorCachePloidyCapacity;
        likelihoodCalculators = new ThreadLocal<GenotypeLikelihoodCalculator[][]>( ) {
            @Override
            public GenotypeLikelihoodCalculator[][] initialValue() {
                return new GenotypeLikelihoodCalculator[calculatorCachePloidyCapacity][calculatorCacheAlleleCapacity];
            }
        };
    }

    @Override
    public <A extends Allele> GenotypingLikelihoods<A> calculateLikelihoods(final AlleleList<A> genotypingAlleles, final GenotypingData<A> data) {
        if (genotypingAlleles == null)
            throw new IllegalArgumentException("the allele cannot be null");
        if (data == null)
            throw new IllegalArgumentException("the genotyping data cannot be null");

        final AlleleListPermutation<A> permutation = AlleleListUtils.permutation(data, genotypingAlleles);
        final AlleleLikelihoodMatrixMapper<A> alleleLikelihoodMatrixMapper = AlleleLikelihoodMatrixMapper.newInstance(permutation);

        final int sampleCount = data.sampleCount();

        switch (sampleCount) {
            case 0: return noSampleLikelihoods(permutation,data);
            case 1: return singleSampleLikelihoods(genotypingAlleles,data,alleleLikelihoodMatrixMapper);
            default:
                final PloidyModel ploidyModel = data.ploidyModel();
                return ploidyModel.isHomogeneous() ? multiSampleHomogeneousPloidyModelLikelihoods(genotypingAlleles, data, alleleLikelihoodMatrixMapper, sampleCount, ploidyModel)
                        : multiSampleHeterogeneousPloidyModelLikelihoods(genotypingAlleles, data, alleleLikelihoodMatrixMapper, sampleCount, ploidyModel);
        }
    }

    private <A extends Allele> GenotypingLikelihoods<A> noSampleLikelihoods(final AlleleList<A> genotypingAlleles,
                                                                            final GenotypingData<A> data) {
        @SuppressWarnings("unchecked")
        final List<GenotypeLikelihoods> likelihoods = Collections.EMPTY_LIST;
        return new GenotypingLikelihoods<>(genotypingAlleles,data.ploidyModel(), likelihoods);

    }

    private <A extends Allele> GenotypingLikelihoods<A> singleSampleLikelihoods(final AlleleList<A> genotypingAlleles,
                                                                                final GenotypingData<A> data,
                                                                                final AlleleLikelihoodMatrixMapper<A> alleleLikelihoodMatrixMapper) {
        final PloidyModel ploidyModel = data.ploidyModel();
        final int samplePloidy = ploidyModel.samplePloidy(0);
        final int alleleCount = genotypingAlleles.alleleCount();
        final GenotypeLikelihoodCalculator likelihoodsCalculator = getLikelihoodsCalculator(samplePloidy,alleleCount);
        final ReadLikelihoods.Matrix<A> sampleLikelihoods = alleleLikelihoodMatrixMapper.map(data.readLikelihoods().sampleMatrix(0));
        final List<GenotypeLikelihoods> genotypeLikelihoods = Collections.singletonList(likelihoodsCalculator.genotypeLikelihoods(sampleLikelihoods));
        return new GenotypingLikelihoods<>(genotypingAlleles,ploidyModel,genotypeLikelihoods);
    }

    private GenotypeLikelihoodCalculator getLikelihoodsCalculator(final int samplePloidy, final int alleleCount) {
        if (samplePloidy >= cacheAlleleCountCapacity)
            return GenotypeLikelihoodCalculators.getInstance(samplePloidy, alleleCount);
        else if (alleleCount >= cacheAlleleCountCapacity)
            return GenotypeLikelihoodCalculators.getInstance(samplePloidy, alleleCount);
        final GenotypeLikelihoodCalculator[][] cache = likelihoodCalculators.get();
        final GenotypeLikelihoodCalculator result = cache[samplePloidy][alleleCount];
        return result != null ? result : (cache[samplePloidy][alleleCount] = GenotypeLikelihoodCalculators.getInstance(samplePloidy, alleleCount));
    }

    private <A extends Allele> GenotypingLikelihoods<A> multiSampleHeterogeneousPloidyModelLikelihoods(final AlleleList<A> genotypingAlleles,
                                                                                                       final GenotypingData<A> data,
                                                                                                       final AlleleLikelihoodMatrixMapper<A> alleleLikelihoodMatrixMapper,
                                                                                                       final int sampleCount,
                                                                                                       final PloidyModel ploidyModel) {
        final List<GenotypeLikelihoods> genotypeLikelihoods = new ArrayList<>(sampleCount);
        final int alleleCount = genotypingAlleles.alleleCount();
        for (int i = 0; i < sampleCount; i++) {
            final int samplePloidy = ploidyModel.samplePloidy(i);
            final GenotypeLikelihoodCalculator likelihoodsCalculator = getLikelihoodsCalculator(samplePloidy,alleleCount);
            final ReadLikelihoods.Matrix<A> sampleLikelihoods = alleleLikelihoodMatrixMapper.map(data.readLikelihoods().sampleMatrix(i));
            genotypeLikelihoods.add(likelihoodsCalculator.genotypeLikelihoods(sampleLikelihoods));
        }
        return new GenotypingLikelihoods<>(genotypingAlleles,ploidyModel,genotypeLikelihoods);
    }

    private <A extends Allele> GenotypingLikelihoods<A> multiSampleHomogeneousPloidyModelLikelihoods(final AlleleList<A> genotypingAlleles,
                                                                                                     final GenotypingData<A> data,
                                                                                                     final AlleleLikelihoodMatrixMapper<A> alleleLikelihoodMatrixMapper,
                                                                                                     final int sampleCount,
                                                                                                     final PloidyModel ploidyModel) {
        final int samplePloidy = ploidyModel.samplePloidy(0);
        final List<GenotypeLikelihoods> genotypeLikelihoods = new ArrayList<>(sampleCount);
        final int alleleCount = genotypingAlleles.alleleCount();
        final GenotypeLikelihoodCalculator likelihoodsCalculator = getLikelihoodsCalculator(samplePloidy,alleleCount);
        for (int i = 0; i < sampleCount; i++) {
            final ReadLikelihoods.Matrix<A> sampleLikelihoods = alleleLikelihoodMatrixMapper.map(data.readLikelihoods().sampleMatrix(i));
            genotypeLikelihoods.add(likelihoodsCalculator.genotypeLikelihoods(sampleLikelihoods));
        }
        return new GenotypingLikelihoods<>(genotypingAlleles,ploidyModel,genotypeLikelihoods);
    }
}
