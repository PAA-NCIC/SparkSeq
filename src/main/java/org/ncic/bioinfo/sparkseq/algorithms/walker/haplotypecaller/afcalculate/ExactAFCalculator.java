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
package org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.afcalculate;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GATKVariantContextUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.MathUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.PriorityQueue;

/**
 * Author: wbc
 */
abstract class ExactAFCalculator extends AFCalculator {

    protected static final int HOM_REF_INDEX = 0;  // AA likelihoods are always first
    /**
     * Sorts {@link ExactAFCalculator.LikelihoodSum} instances where those with higher likelihood are first.
     */
    protected static final Comparator<LikelihoodSum> LIKELIHOOD_SUM_COMPARATOR = new Comparator<LikelihoodSum>() {

        @Override
        public int compare(final LikelihoodSum o1, final LikelihoodSum o2) {
            return - Double.compare(o1.sum,o2.sum);
        }
    };
    /**
     * Sorts {@link ExactAFCalculator.LikelihoodSum} instances where those with higher likelihood are first but make sure that
     * NON_REF alleles are place are last.
     */
    protected static final Comparator<LikelihoodSum> LIKELIHOOD_NON_REF_THEN_SUM_COMPARATOR = new Comparator<LikelihoodSum>() {
        @Override
        public int compare(final LikelihoodSum o1, final LikelihoodSum o2) {
            if (o1.allele == GATKVariantContextUtils.NON_REF_SYMBOLIC_ALLELE)
                return 1;
            else if (o2.allele == GATKVariantContextUtils.NON_REF_SYMBOLIC_ALLELE)
                return -1;
            else
                return o1.compareTo(o2);
        }
    };
    /**
     * Sorts {@link ExactAFCalculator.LikelihoodSum} instances where those with lower alternative allele index are first regardless of
     * the likelihood sum.
     */
    protected static final Comparator<LikelihoodSum> LIKELIHOOD_INDEX_COMPARATOR = new Comparator<LikelihoodSum>() {
        @Override
        public int compare(final LikelihoodSum o1, final LikelihoodSum o2) {
            return Integer.compare(o1.index, o2.index);
        }
    };

    protected ExactAFCalculator() {

    }

    /**
     * Wrapper class that compares two likelihoods associated with two alleles
     */
    protected static final class LikelihoodSum implements Comparable<LikelihoodSum> {
        public double sum = 0.0;
        public final Allele allele;
        public final int index;

        public LikelihoodSum(final Allele allele, final int index) { this.allele = allele; this.index = index; }

        public int compareTo(LikelihoodSum other) {
            final double diff = sum - other.sum;
            return ( diff < 0.0 ) ? 1 : (diff > 0.0 ) ? -1 : 0;
        }
    }

    /**
     * Unpack GenotypesContext into arraylist of doubel values
     * @param GLs            Input genotype context
     * @return               ArrayList of doubles corresponding to GL vectors
     */
    protected static ArrayList<double[]> getGLs(final GenotypesContext GLs, final boolean includeDummy) {
        final ArrayList<double[]> genotypeLikelihoods = new ArrayList<>(GLs.size() + 1);

        if ( includeDummy ) genotypeLikelihoods.add(new double[]{0.0,0.0,0.0}); // dummy
        for ( Genotype sample : GLs.iterateInSampleNameOrder() ) {
            if ( sample.hasLikelihoods() ) {
                final double[] gls = sample.getLikelihoods().getAsVector();

                if ( MathUtils.sum(gls) < GATKVariantContextUtils.SUM_GL_THRESH_NOCALL )
                    genotypeLikelihoods.add(gls);
            }
        }

        return genotypeLikelihoods;
    }

    @Override
    protected VariantContext reduceScope(final VariantContext vc, final int defaultPloidy, final int maximumAlternativeAlleles) {
        // don't try to genotype too many alternate alleles
        final List<Allele> inputAltAlleles = vc.getAlternateAlleles();
        final List<Allele> outputAltAlleles = reduceScopeAlleles(vc, defaultPloidy, maximumAlternativeAlleles);

        // only if output allele has reduced from the input alt allele set size we should care.
        final int altAlleleReduction = inputAltAlleles.size() - outputAltAlleles.size();

        if (altAlleleReduction == 0)
            return vc;
        else if (altAlleleReduction != 0) {
            logger.warn("this tool is currently set to genotype at most " + maximumAlternativeAlleles
                    + " alternate alleles in a given context, but the context at " + vc.getChr() + ":" + vc.getStart()
                    + " has " + (vc.getAlternateAlleles().size())
                    + " alternate alleles so only the top alleles will be used; see the --max_alternate_alleles argument");

            final List<Allele> alleles = new ArrayList<>(maximumAlternativeAlleles + 1);
            alleles.add(vc.getReference());
            alleles.addAll(reduceScopeAlleles(vc, defaultPloidy, maximumAlternativeAlleles));
            final VariantContextBuilder builder = new VariantContextBuilder(vc);
            builder.alleles(alleles);
            builder.genotypes(reduceScopeGenotypes(vc, defaultPloidy, alleles));
            if (altAlleleReduction < 0)
                throw new IllegalStateException("unexpected: reduction increased the number of alt. alleles!: " + - altAlleleReduction + " " + vc + " " + builder.make());
            return builder.make();
        } else // if (altAlleleReduction < 0)
            throw new IllegalStateException("unexpected: reduction increased the number of alt. alleles!: " + - altAlleleReduction + " " + vc);
    }

    /**
     * Returns a the new set of alleles to use.
     * @param vc target variant context.
     * @param numAllelesToChoose number of alleles to keep.
     * @return the list of alternative allele to keep.
     */
    protected List<Allele> reduceScopeAlleles(final VariantContext vc, final int defaultPloidy, final int numAllelesToChoose) {

        // Look  for the <NON_REF> allele to exclude it from the pruning if present.
        final int numOriginalAltAlleles = vc.getAlternateAlleles().size();

        final int nonRefAltAlleleIndex = GATKVariantContextUtils.indexOfAltAllele(vc,
                GATKVariantContextUtils.NON_REF_SYMBOLIC_ALLELE, false);
        final boolean nonRefAltAllelePresent = nonRefAltAlleleIndex >= 0;

        // <NON_REF> should not be considered in the downsizing, so we need to count it out when
        // considering if alt. allele downsizing is required.
        final int numProperOriginalAltAlleles = numOriginalAltAlleles - (nonRefAltAllelePresent ? 1 : 0);

        // Avoid pointless allele reduction:
        if (numAllelesToChoose >= numProperOriginalAltAlleles)
            return vc.getAlternateAlleles();

        final LikelihoodSum[] likelihoodSums = new LikelihoodSum[numOriginalAltAlleles];
        for ( int i = 0; i < numOriginalAltAlleles; i++ ) {
            final Allele allele = vc.getAlternateAllele(i);
            likelihoodSums[i] = new LikelihoodSum(allele,i);
        }

        // Calculate the allele likelihood sums.
        reduceScopeCalculateLikelihoodSums(vc, defaultPloidy, likelihoodSums);

        // sort them by probability mass and choose the best ones
        // Make sure that the <NON_REF> allele is last if present.
        Collections.sort(Arrays.asList(likelihoodSums), nonRefAltAllelePresent ? LIKELIHOOD_NON_REF_THEN_SUM_COMPARATOR : LIKELIHOOD_SUM_COMPARATOR);

        // We need to return the best likelihood alleles in the original alternative allele index order.
        // This heap will keep track of that index order.
        final PriorityQueue<LikelihoodSum> mostLikelyAllelesHeapByIndex = new PriorityQueue<>(numOriginalAltAlleles, LIKELIHOOD_INDEX_COMPARATOR);

        for ( int i = 0; i < numAllelesToChoose; i++ )
            mostLikelyAllelesHeapByIndex.add(likelihoodSums[i]);

        // guaranteed no to have been added at this point thanks for checking on whether reduction was
        // needed in the first place.
        if (nonRefAltAllelePresent)
            mostLikelyAllelesHeapByIndex.add(likelihoodSums[nonRefAltAlleleIndex]);

        final ArrayList<Allele> orderedBestAlleles = new ArrayList<>(numAllelesToChoose);

        while (!mostLikelyAllelesHeapByIndex.isEmpty())
            orderedBestAlleles.add(mostLikelyAllelesHeapByIndex.remove().allele);

        return orderedBestAlleles;
    }

    protected static final int PL_INDEX_OF_HOM_REF = 0;

    /**
     * Update the likelihood sums with using the variant context genotype likelihoods.
     * @param vc source variant context.
     * @param likelihoodSums where to update the likelihood sums.
     */
    protected abstract void reduceScopeCalculateLikelihoodSums(final VariantContext vc, final int defaultPloidy, final LikelihoodSum[] likelihoodSums);

    /**
     * Transforms the genotypes of the variant context according to the new subset of possible alleles.
     *
     * @param vc original variant-context.
     * @param allelesToUse possible alleles.
     * @return never {@code null}, the new set of genotype calls for the reduced scope.
     */
    protected abstract GenotypesContext reduceScopeGenotypes(final VariantContext vc, final int defaultPloidy, final List<Allele> allelesToUse);
}
