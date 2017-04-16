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
package org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.argcollection;

import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.HomoSapiensConstants;

import java.util.Collections;
import java.util.List;

/**
 * Author: wbc
 */
public class GenotypeCalculationArgumentCollection implements Cloneable {

    /**
     * Depending on the value of the --max_alternate_alleles argument, we may genotype only a fraction of the alleles being sent on for genotyping.
     * Using this argument instructs the genotyper to annotate (in the INFO field) the number of alternate alleles that were originally discovered at the site.
     */
    public boolean ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED = false;

    /**
     * The expected heterozygosity value used to compute prior probability that a locus is non-reference.
     * <p>
     * The default priors are for provided for humans:
     * <p>
     * het = 1e-3
     * <p>
     * which means that the probability of N samples being hom-ref at a site is:
     * <p>
     * 1 - sum_i_2N (het / i)
     * <p>
     * Note that heterozygosity as used here is the population genetics concept:
     * <p>
     * http://en.wikipedia.org/wiki/Zygosity#Heterozygosity_in_population_genetics
     * <p>
     * That is, a hets value of 0.01 implies that two randomly chosen chromosomes from the population of organisms
     * would differ from each other (one being A and the other B) at a rate of 1 in 100 bp.
     * <p>
     * Note that this quantity has nothing to do with the likelihood of any given sample having a heterozygous genotype,
     * which in the GATK is purely determined by the probability of the observed data P(D | AB) under the model that there
     * may be a AB het genotype.  The posterior probability of this AB genotype would use the het prior, but the GATK
     * only uses this posterior probability in determining the prob. that a site is polymorphic.  So changing the
     * het parameters only increases the chance that a site will be called non-reference across all samples, but
     * doesn't actually change the output genotype likelihoods at all, as these aren't posterior probabilities at all.
     * <p>
     * The quantity that changes whether the GATK considers the possibility of a het genotype at all is the ploidy,
     * which determines how many chromosomes each individual in the species carries.
     */
    public Double snpHeterozygosity = HomoSapiensConstants.SNP_HETEROZYGOSITY;

    /**
     * This argument informs the prior probability of having an indel at a site.
     */
    public double indelHeterozygosity = HomoSapiensConstants.INDEL_HETEROZYGOSITY;

    /**
     * The minimum phred-scaled Qscore threshold to separate high confidence from low confidence calls. Only genotypes with
     * confidence >= this threshold are emitted as called sites. A reasonable threshold is 30 for high-pass calling (this
     * is the default).
     */
    public double STANDARD_CONFIDENCE_FOR_CALLING = 30.0;

    /**
     * This argument allows you to emit low quality calls as filtered records.
     */
    public double STANDARD_CONFIDENCE_FOR_EMITTING = 30.0;

    /**
     * If there are more than this number of alternate alleles presented to the genotyper (either through discovery or GENOTYPE_GIVEN ALLELES),
     * then only this many alleles will be used.  Note that genotyping sites with many alternate alleles is both CPU and memory intensive and it
     * scales exponentially based on the number of alternate alleles.  Unless there is a good reason to change the default value, we highly recommend
     * that you not play around with this parameter.
     * <p>
     * As of GATK 2.2 the genotyper can handle a very large number of events, so the default maximum has been increased to 6.
     */
    public int MAX_ALTERNATE_ALLELES = 6;

    /**
     * By default, the prior specified with the argument --heterozygosity/-hets is used for variant discovery at a particular locus, using an infinite sites model,
     * see e.g. Waterson (1975) or Tajima (1996).
     * This model asserts that the probability of having a population of k variant sites in N chromosomes is proportional to theta/k, for 1=1:N
     * <p>
     * There are instances where using this prior might not be desireable, e.g. for population studies where prior might not be appropriate,
     * as for example when the ancestral status of the reference allele is not known.
     * By using this argument, user can manually specify priors to be used for calling as a vector for doubles, with the following restriciotns:
     * a) User must specify 2N values, where N is the number of samples.
     * b) Only diploid calls supported.
     * c) Probability values are specified in double format, in linear space.
     * d) No negative values allowed.
     * e) Values will be added and Pr(AC=0) will be 1-sum, so that they sum up to one.
     * f) If user-defined values add to more than one, an error will be produced.
     * <p>
     * If user wants completely flat priors, then user should specify the same value (=1/(2*N+1)) 2*N times,e.g.
     * -inputPrior 0.33 -inputPrior 0.33
     * for the single-sample diploid case.
     */
    public List<Double> inputPrior = Collections.emptyList();

    /**
     * Sample ploidy - equivalent to number of chromosomes per pool. In pooled experiments this should be = # of samples in pool * individual sample ploidy
     */
    public int samplePloidy = HomoSapiensConstants.DEFAULT_PLOIDY;

    /**
     * Creates a copy of this configuration.
     *
     * @return never {@code null}.
     */
    @Override
    public GenotypeCalculationArgumentCollection clone() {
        try {
            return (GenotypeCalculationArgumentCollection) super.clone();
        } catch (CloneNotSupportedException e) {
            throw new IllegalStateException("unreachable code");
        }
    }
}
