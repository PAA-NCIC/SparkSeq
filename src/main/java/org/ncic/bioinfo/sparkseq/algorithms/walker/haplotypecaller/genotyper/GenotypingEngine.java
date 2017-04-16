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
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.log4j.Logger;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GATKVariantContextUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLoc;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLocParser;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenotypingGivenAllelesUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.MathUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.QualityUtils;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefMetaDataTracker;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.ReferenceContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.AlignmentContext;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.GenotypingOutputMode;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.OutputMode;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.PerReadAlleleLikelihoodMap;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.SampleList;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.afcalculate.AFCalculationResult;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.afcalculate.AFCalculator;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.afcalculate.AFCalculatorProvider;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.VariantAnnotatorEngine;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.argcollection.StandardCallerArgumentCollection;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.model.GenotypeLikelihoodsCalculationModel;
import org.ncic.bioinfo.sparkseq.exceptions.UserException;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Author: wbc
 */
public abstract class GenotypingEngine<Config extends StandardCallerArgumentCollection> {

    public static final String NUMBER_OF_DISCOVERED_ALLELES_KEY = "NDA";

    public static final String LOW_QUAL_FILTER_NAME = "LowQual";

    protected final AFCalculatorProvider afCalculatorProvider   ;

    protected Logger logger;

    protected final Config configuration;

    protected VariantAnnotatorEngine annotationEngine;

    protected final int numberOfGenomes;

    protected final SampleList samples;


    private final AFPriorProvider log10AlleleFrequencyPriorsSNPs;

    private final AFPriorProvider log10AlleleFrequencyPriorsIndels;

    protected final GenomeLocParser genomeLocParser;

    /**
     * Construct a new genotyper engine, on a specific subset of samples.
     *
     * @param configuration engine configuration object.
     * @param samples subset of sample to work on identified by their names. If {@code null}, the full toolkit
     *                    sample set will be used instead.
     * @param genomeLocParser the genome-loc-parser
     *
     * @throws IllegalArgumentException if any of {@code samples}, {@code configuration} or {@code genomeLocParser} is {@code null}.
     */
    protected GenotypingEngine(final Config configuration, final SampleList samples,
                               final GenomeLocParser genomeLocParser, final AFCalculatorProvider afCalculatorProvider) {

        if (configuration == null)
            throw new IllegalArgumentException("the configuration cannot be null");
        if (samples == null)
            throw new IllegalArgumentException("the sample list provided cannot be null");
        if (afCalculatorProvider == null)
            throw new IllegalArgumentException("the AF calculator provider cannot be null");

        this.afCalculatorProvider = afCalculatorProvider;
        this.configuration = configuration;
        logger = Logger.getLogger(getClass());
        this.samples = samples;
        numberOfGenomes = this.samples.sampleCount() * configuration.genotypeArgs.samplePloidy;
        MathUtils.Log10Cache.ensureCacheContains(numberOfGenomes * 2);
        log10AlleleFrequencyPriorsSNPs = composeAlleleFrequencyPriorProvider(numberOfGenomes,
                configuration.genotypeArgs.snpHeterozygosity, configuration.genotypeArgs.inputPrior);
        log10AlleleFrequencyPriorsIndels = composeAlleleFrequencyPriorProvider(numberOfGenomes,
                configuration.genotypeArgs.indelHeterozygosity, configuration.genotypeArgs.inputPrior);
        this.genomeLocParser = genomeLocParser;
    }

    /**
     * Changes the logger for this genotyper engine.
     *
     * @param logger new logger.
     *
     * @throws IllegalArgumentException if {@code logger} is {@code null}.
     */
    public void setLogger(final Logger logger) {
        if (logger == null)
            throw new IllegalArgumentException("the logger cannot be null");
        this.logger = logger;
    }

    public Set<VCFInfoHeaderLine> getAppropriateVCFInfoHeaders() {
        Set<VCFInfoHeaderLine> headerInfo = new HashSet<>();
        if ( configuration.genotypeArgs.ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED )
            headerInfo.add(new VCFInfoHeaderLine(UnifiedGenotypingEngine.NUMBER_OF_DISCOVERED_ALLELES_KEY, 1, VCFHeaderLineType.Integer, "Number of alternate alleles discovered (but not necessarily genotyped) at this site"));
        return headerInfo;
    }

    /**
     * Returns a reference to the engine configuration
     *
     * @return never {@code null}.
     */
    public Config getConfiguration() {
        return configuration;
    }

    /**
     * Completes a variant context with genotype calls and associated annotations given the genotype likelihoods and
     *  the model that need to be applied.
     *
     * @param vc variant-context to complete.
     * @param model model name.
     *
     * @throws IllegalArgumentException if {@code model} or {@code vc} is {@code null}.
     *
     * @return can be {@code null} indicating that genotyping it not possible with the information provided.
     */
    public VariantCallContext calculateGenotypes(final VariantContext vc, final GenotypeLikelihoodsCalculationModel.Model model) {
        if (vc == null)
            throw new IllegalArgumentException("vc cannot be null");
        if (model == null)
            throw new IllegalArgumentException("the model cannot be null");
        return calculateGenotypes(null,null,null,null,vc,model,false,null);
    }

    /**
     * Main entry function to calculate genotypes of a given VC with corresponding GL's that is shared across genotypers (namely UG and HC).
     *
     * @param tracker                            Tracker
     * @param refContext                         Reference context
     * @param rawContext                         Raw context
     * @param stratifiedContexts                 Stratified alignment contexts
     * @param vc                                 Input VC
     * @param model                              GL calculation model
     * @param inheritAttributesFromInputVC       Output VC will contain attributes inherited from input vc
     * @return                                   VC with assigned genotypes
     */
    protected VariantCallContext calculateGenotypes(final RefMetaDataTracker tracker, final ReferenceContext refContext,
                                                    final AlignmentContext rawContext, Map<String, AlignmentContext> stratifiedContexts,
                                                    final VariantContext vc, final GenotypeLikelihoodsCalculationModel.Model model,
                                                    final boolean inheritAttributesFromInputVC,
                                                    final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {

        final boolean limitedContext = tracker == null || refContext == null || rawContext == null || stratifiedContexts == null;
        // if input VC can't be genotyped, exit with either null VCC or, in case where we need to emit all sites, an empty call
        if (hasTooManyAlternativeAlleles(vc) || vc.getNSamples() == 0)
            return emptyCallContext(tracker,refContext,rawContext);

        final int defaultPloidy = configuration.genotypeArgs.samplePloidy;
        final int maxAltAlleles = configuration.genotypeArgs.MAX_ALTERNATE_ALLELES;
        final AFCalculator afCalculator = afCalculatorProvider.getInstance(vc,defaultPloidy,maxAltAlleles);
        final AFCalculationResult AFresult = afCalculator.getLog10PNonRef(vc, defaultPloidy,maxAltAlleles, getAlleleFrequencyPriors(vc,defaultPloidy,model));

        final OutputAlleleSubset outputAlternativeAlleles = calculateOutputAlleleSubset(AFresult);

        final double PoFGT0 = Math.pow(10, AFresult.getLog10PosteriorOfAFGT0());

        // note the math.abs is necessary because -10 * 0.0 => -0.0 which isn't nice
        double log10Confidence =
                ! outputAlternativeAlleles.siteIsMonomorphic ||
                        configuration.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES || configuration.annotateAllSitesWithPLs
                        ? AFresult.getLog10PosteriorOfAFEq0() + 0.0
                        : AFresult.getLog10PosteriorOfAFGT0() + 0.0 ;


        // Add 0.0 removes -0.0 occurrences.
        final double phredScaledConfidence = (-10.0 * log10Confidence) + 0.0;

        // return a null call if we don't pass the confidence cutoff or the most likely allele frequency is zero
        if ( !passesEmitThreshold(phredScaledConfidence, outputAlternativeAlleles.siteIsMonomorphic) && !forceSiteEmission())
            // technically, at this point our confidence in a reference call isn't accurately estimated
            //  because it didn't take into account samples with no data, so let's get a better estimate
            return limitedContext ? null : estimateReferenceConfidence(vc, stratifiedContexts, getModelTheta(model), true, PoFGT0);

        // start constructing the resulting VC
        final GenomeLocParser genomeLocParser = this.genomeLocParser;
        if (genomeLocParser == null)
            throw new IllegalStateException("this UG engine was created without a valid genomeLocParser and no refContext was provided");
        final GenomeLoc loc = genomeLocParser.createGenomeLoc(vc);
        final List<Allele> outputAlleles = outputAlternativeAlleles.outputAlleles(vc.getReference());
        final VariantContextBuilder builder = new VariantContextBuilder(callSourceString(), loc.getContig(), loc.getStart(), loc.getStop(), outputAlleles);

        // Seems that when log10PError is 0.0, you must pass -0.0 to get a nice output at the other end otherwise is a "-0".
        // Truth is that this should be fixed in the "variant" dependency code but perhaps it can be amended also in the VariantContextWriter.
        //TODO Please remove this comment when this has been fixed (PT https://www.pivotaltracker.com/story/show/69492530)
        //TODO and change the code below accordingly.
        builder.log10PError(log10Confidence == 0.0 ? -0.0 : log10Confidence);
        if ( ! passesCallThreshold(phredScaledConfidence) )
            builder.filter(LOW_QUAL_FILTER_NAME);

        // create the genotypes

        final GenotypesContext genotypes = afCalculator.subsetAlleles(vc, defaultPloidy, outputAlleles, true);
        builder.genotypes(genotypes);

        // *** note that calculating strand bias involves overwriting data structures, so we do that last
        final Map<String, Object> attributes = composeCallAttributes(inheritAttributesFromInputVC, vc, rawContext, stratifiedContexts, tracker, refContext,
                outputAlternativeAlleles.alternativeAlleleMLECounts(), outputAlternativeAlleles.siteIsMonomorphic, AFresult, outputAlternativeAlleles.outputAlleles(vc.getReference()),genotypes,model,perReadAlleleLikelihoodMap);

        builder.attributes(attributes);

        VariantContext vcCall = builder.make();

        /*if ( annotationEngine != null && !limitedContext ) { // limitedContext callers need to handle annotations on their own by calling their own annotationEngine
            // Note: we want to use the *unfiltered* and *unBAQed* context for the annotations
            final ReadBackedPileup pileup = rawContext.getBasePileup();
            stratifiedContexts = AlignmentContextUtils.splitContextBySampleName(pileup);

            vcCall = annotationEngine.annotateContext(tracker, refContext, stratifiedContexts, vcCall, perReadAlleleLikelihoodMap);
        }*/

        // if we are subsetting alleles (either because there were too many or because some were not polymorphic)
        // then we may need to trim the alleles (because the original VariantContext may have had to pad at the end).
        if ( outputAlleles.size() != vc.getAlleles().size() && !limitedContext ) // limitedContext callers need to handle allele trimming on their own to keep their perReadAlleleLikelihoodMap alleles in sync
            vcCall = GATKVariantContextUtils.reverseTrimAlleles(vcCall);

        return new VariantCallContext(vcCall, confidentlyCalled(phredScaledConfidence, PoFGT0));
    }

    /**
     * What string to use as source of variant-context generated by this genotyper-engine.
     * @return never {@code null} nor empty.
     */
    protected abstract String callSourceString();

    /**
     * Holds information about the alternative allele subsetting based on supporting evidence, genotyping and
     * output modes.
     */
    private static class OutputAlleleSubset {
        private  final Allele[] alleles;
        private  final boolean siteIsMonomorphic;
        private  final int[] mleCounts;
        private  final int count;

        private OutputAlleleSubset(final int count, final Allele[] alleles, final int[] mleCounts, final boolean siteIsMonomorphic) {
            this.siteIsMonomorphic = siteIsMonomorphic;
            this.count = count;
            this.alleles = alleles;
            this.mleCounts = mleCounts;
        }

        private List<Allele> outputAlleles(final Allele referenceAllele) {
            final ArrayList<Allele> result = new ArrayList<>(count + 1);
            result.add(referenceAllele);
            for (int i = 0; i < count; i++)
                result.add(alleles[i]);
            return result;
        }

        public List<Integer> alternativeAlleleMLECounts() {
            final List<Integer> result = new ArrayList<>(count);
            for (int i = 0; i < count; i++)
                result.add(mleCounts[i]);
            return result;
        }
    }


    /**
     * Provided the exact mode computations it returns the appropiate subset of alleles that progress to genotyping.
     * @param afcr the exact model calcualtion result.
     * @return never {@code null}.
     */
    private OutputAlleleSubset calculateOutputAlleleSubset(final AFCalculationResult afcr) {
        final List<Allele> alleles = afcr.getAllelesUsedInGenotyping();

        final int alternativeAlleleCount = alleles.size() - 1;
        Allele[] outputAlleles = new Allele[alternativeAlleleCount];
        int[] mleCounts = new int[alternativeAlleleCount];
        int outputAlleleCount = 0;
        boolean siteIsMonomorphic = true;
        for (final Allele alternativeAllele : alleles) {
            if (alternativeAllele.isReference()) continue;
            final boolean isPlausible = afcr.isPolymorphicPhredScaledQual(alternativeAllele, configuration.genotypeArgs.STANDARD_CONFIDENCE_FOR_EMITTING);
            final boolean toOutput = isPlausible || forceKeepAllele(alternativeAllele);

            siteIsMonomorphic &= ! isPlausible;
            if (!toOutput) continue;
            outputAlleles[outputAlleleCount] = alternativeAllele;
            mleCounts[outputAlleleCount++] = afcr.getAlleleCountAtMLE(alternativeAllele);
        }

        return new OutputAlleleSubset(outputAlleleCount,outputAlleles,mleCounts,siteIsMonomorphic);
    }

    /**
     * Checks whether even if the allele is not well supported by the data, we should keep it for genotyping.
     *
     * @param allele target allele.
     *
     * @return {@code true} iff we need to keep this alleles even if does not seem plausible.
     */
    protected abstract boolean forceKeepAllele(final Allele allele);

    /**
     * Checks whether a variant site seems confidently called base on user threshold that the score provided
     * by the exact model.
     *
     * @param conf
     * @param PofF
     * @return {@code true} iff the variant is confidently called.
     */
    protected final boolean confidentlyCalled(final double conf, final double PofF) {
        return conf >= configuration.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING ||
                (configuration.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES
                        && QualityUtils.phredScaleErrorRate(PofF) >= configuration.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING);
    }

    /**
     * Based in the model used, returns the appropriate heterozygosity argument value.
     * @param model genotyping model.
     *
     * @return a valid heterozygosity in (0,1).
     */
    private double getModelTheta(final GenotypeLikelihoodsCalculationModel.Model model) {
        switch (model) {
            case SNP:
            case GENERALPLOIDYSNP:
                return configuration.genotypeArgs.snpHeterozygosity;
            case INDEL:
            case GENERALPLOIDYINDEL:
                return configuration.genotypeArgs.indelHeterozygosity;
            default:
                throw new IllegalArgumentException("Unexpected GenotypeCalculationModel " + model);
        }
    }


    /**
     * Checks whether the variant context has too many alternative alleles for progress to genotyping the site.
     * <p>
     *     AF calculation may get intro trouble with too many alternative alleles.
     * </p>
     *
     * @param vc the variant context to evaluate.
     *
     * @throws NullPointerException if {@code vc} is {@code null}.
     *
     * @return {@code true} iff there is too many alternative alleles based on
     * {@link GenotypeLikelihoods#MAX_ALT_ALLELES_THAT_CAN_BE_GENOTYPED}.
     */
    protected final boolean hasTooManyAlternativeAlleles(final VariantContext vc) {
        // protect against too many alternate alleles that we can't even run AF on:
        if (vc.getNAlleles() <= GenotypeLikelihoods.MAX_ALT_ALLELES_THAT_CAN_BE_GENOTYPED)
            return false;
        logger.warn("Attempting to genotype more than "+GenotypeLikelihoods.MAX_ALT_ALLELES_THAT_CAN_BE_GENOTYPED +
                " alleles. Site will be skipped at location "+vc.getChr()+":"+vc.getStart());
        return true;
    }

    /**
     * Produces an empty variant-call context to output when there is no enough data provided to call anything.
     *
     * @param tracker meta-data tracker.
     * @param ref the reference context.
     * @param rawContext the read alignment at that location.
     * @return it might be null if no enough information is provided to do even an empty call. For example when
     * we have limited-context (i.e. any of the tracker, reference or alignment is {@code null}.
     */
    protected final VariantCallContext emptyCallContext(final RefMetaDataTracker tracker, final ReferenceContext ref,
                                                        final AlignmentContext rawContext) {
        if (tracker == null || ref == null || rawContext == null)
            return null;

        if (!forceSiteEmission())
            return null;

        VariantContext vc;

        if ( configuration.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES ) {
            final VariantContext ggaVc = GenotypingGivenAllelesUtils.composeGivenAllelesVariantContextFromRod(tracker,
                    rawContext.getLocation(), false, logger, configuration.alleles);
            if (ggaVc == null)
                return null;
            vc = new VariantContextBuilder(callSourceString(), ref.getLocus().getContig(), ggaVc.getStart(),
                    ggaVc.getEnd(), ggaVc.getAlleles()).make();
        } else {
            // deal with bad/non-standard reference bases
            if ( !Allele.acceptableAlleleBases(new byte[]{ref.getBase()}) )
                return null;
            final Set<Allele> alleles = new HashSet<>(Collections.singleton(Allele.create(ref.getBase(),true)));
            vc = new VariantContextBuilder(callSourceString(), ref.getLocus().getContig(),
                    ref.getLocus().getStart(), ref.getLocus().getStart(), alleles).make();
        }

        /*
        if ( vc != null && annotationEngine != null ) {
            // Note: we want to use the *unfiltered* and *unBAQed* context for the annotations
            final ReadBackedPileup pileup = rawContext.getBasePileup();
            vc = annotationEngine.annotateContext(tracker, ref,  AlignmentContextUtils.splitContextBySampleName(pileup), vc);
        }*/

        return new VariantCallContext(vc, false);
    }

    /**
     * Indicates whether we have to emit any site no matter what.
     * <p>
     *     Note: this has been added to allow differences between UG and HC GGA modes where the latter force emmitions of all given alleles
     *     sites even if there is no enough confidence.
     * </p>
     *
     * @return {@code true} iff we force emissions.
     */
    protected abstract boolean forceSiteEmission();

    protected final VariantCallContext estimateReferenceConfidence(VariantContext vc, Map<String, AlignmentContext> contexts, double theta, boolean ignoreCoveredSamples, double initialPofRef) {
        if ( contexts == null )
            return null;

        double log10POfRef = Math.log10(initialPofRef);

        // for each sample that we haven't examined yet
        final int sampleCount = samples.sampleCount();
        for (int i = 0; i < sampleCount; i++) {
            final String sample = samples.sampleAt(i);
            final AlignmentContext context = contexts.get(sample);
            if ( ignoreCoveredSamples && context != null )
                continue;
            final int depth = context == null ? 0 : context.getBasePileup().depthOfCoverage();
            log10POfRef += estimateLog10ReferenceConfidenceForOneSample(depth, theta);
        }

        return new VariantCallContext(vc, QualityUtils.phredScaleLog10CorrectRate(log10POfRef) >= configuration.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING, false);
    }

    /**
     * Returns the log10 prior probability for all possible allele counts from 0 to N where N is the total number of
     * genomes (total-ploidy).
     *
     * @param vc the target variant-context, use to determine the total ploidy thus the possible ACs.
     * @param defaultPloidy default ploidy to be assume if we do not have the ploidy for some sample in {@code vc}.
     * @param model the calculation model (SNP,INDEL or MIXED) whose priors are to be retrieved.
     * @throws java.lang.NullPointerException if either {@code vc} or {@code model} is {@code null}
     * @return never {@code null}, an array with exactly <code>total-ploidy(vc) + 1</code> positions.
     */
    protected final double[] getAlleleFrequencyPriors( final VariantContext vc, final int defaultPloidy, final GenotypeLikelihoodsCalculationModel.Model model ) {
        final int totalPloidy = GATKVariantContextUtils.totalPloidy(vc,defaultPloidy);
        switch (model) {
            case SNP:
            case GENERALPLOIDYSNP:
                return log10AlleleFrequencyPriorsSNPs.forTotalPloidy(totalPloidy);
            case INDEL:
            case GENERALPLOIDYINDEL:
                return log10AlleleFrequencyPriorsIndels.forTotalPloidy(totalPloidy);
            default:
                throw new IllegalArgumentException("Unexpected GenotypeCalculationModel " + model);
        }
    }

    /**
     * Compute the log10 probability of a sample with sequencing depth and no alt allele is actually truly homozygous reference
     *
     * Assumes the sample is diploid
     *
     * @param depth the depth of the sample
     * @param theta the heterozygosity of this species (between 0 and 1)
     *
     * @throws IllegalArgumentException if {@code depth} is less than 0 or {@code theta} is not in the (0,1) range.
     *
     * @return a valid log10 probability of the sample being hom-ref
     */
    protected final double estimateLog10ReferenceConfidenceForOneSample(final int depth, final double theta) {
        if (theta <= 0 || theta >= 1) throw new IllegalArgumentException("theta must be greater than 0 and less than 1");
        final double log10PofNonRef = Math.log10(theta / 2.0) + getRefBinomialProbLog10(depth);
        return MathUtils.log10OneMinusX(Math.pow(10.0, log10PofNonRef));
    }

    /**
     * Calculates the hom-reference binomial log10 probability given the depth.
     *
     * @param depth the query depth.
     *
     * @throws IllegalArgumentException if {@code depth} is less than 0.
     *
     * @return a valid log10 probability between 0 and {@link Double#NEGATIVE_INFINITY}.
     */
    protected final double getRefBinomialProbLog10(final int depth) {
        if (depth < 0)
            throw new IllegalArgumentException("depth cannot be less than 0");
        return MathUtils.log10BinomialProbability(depth, 0);
    }

    /**
     * Function that fills vector with allele frequency priors. By default, infinite-sites, neutral variation prior is used,
     * where Pr(AC=i) = theta/i where theta is heterozygosity
     * @param N                                Number of chromosomes
     * @param heterozygosity                   default heterozygosity to use, if inputPriors is empty
     * @param inputPriors                      Input priors to use (in which case heterozygosity is ignored)
     *
     * @return never {@code null}.
     */
    private static AFPriorProvider composeAlleleFrequencyPriorProvider(final int N, final double heterozygosity, final List<Double> inputPriors) {

        final double[] priors = new double[N + 1];
        double sum = 0.0;
        final AFPriorProvider result;

        if (!inputPriors.isEmpty()) {
            // user-specified priors
            if (inputPriors.size() != N)
                throw new UserException.BadArgumentValue("inputPrior","Invalid length of inputPrior vector: vector length must be equal to # samples +1 ");
            return new CustomAFPriorProvider(inputPriors);
        }
        else
            return new HeterozygosityAFPriorProvider(heterozygosity);
    }

    /**
     * Function that fills vector with allele frequency priors. By default, infinite-sites, neutral variation prior is used,
     * where Pr(AC=i) = theta/i where theta is heterozygosity
     * @param N                                Number of chromosomes
     * @param priors                           (output) array to be filled with priors
     * @param heterozygosity                   default heterozygosity to use, if inputPriors is empty
     * @param inputPriors                      Input priors to use (in which case heterozygosity is ignored)
     */
    public static void computeAlleleFrequencyPriors(final int N, final double[] priors, final double heterozygosity, final List<Double> inputPriors) {


        double sum = 0.0;

        if (!inputPriors.isEmpty()) {
            // user-specified priors
            if (inputPriors.size() != N)
                throw new UserException.BadArgumentValue("inputPrior","Invalid length of inputPrior vector: vector length must be equal to # samples +1 ");

            int idx = 1;
            for (final double prior: inputPriors) {
                if (prior < 0.0)
                    throw new UserException.BadArgumentValue("Bad argument: negative values not allowed","inputPrior");
                priors[idx++] = Math.log10(prior);
                sum += prior;
            }
        }
        else {
            // for each i
            for (int i = 1; i <= N; i++) {
                final double value = heterozygosity / (double)i;
                priors[i] = Math.log10(value);
                sum += value;
            }
        }

        // protection against the case of heterozygosity too high or an excessive number of samples (which break population genetics assumptions)
        if (sum > 1.0) {
            throw new UserException.BadArgumentValue("heterozygosity","The heterozygosity value is set too high relative to the number of samples to be processed, or invalid values specified if input priors were provided - try reducing heterozygosity value or correct input priors.");
        }
        // null frequency for AF=0 is (1 - sum(all other frequencies))
        priors[0] = Math.log10(1.0 - sum);
    }


    protected final boolean passesEmitThreshold(double conf, boolean bestGuessIsRef) {
        return (configuration.outputMode == OutputMode.EMIT_ALL_CONFIDENT_SITES || !bestGuessIsRef) &&
                conf >= Math.min(configuration.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING,
                        configuration.genotypeArgs.STANDARD_CONFIDENCE_FOR_EMITTING);
    }

    protected final boolean passesCallThreshold(double conf) {
        return conf >= configuration.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING;
    }

    private static final String DOWNSAMPLED_KEY = "DS";
    private static final String MLE_ALLELE_COUNT_KEY = "MLEAC";
    private static final String MLE_ALLELE_FREQUENCY_KEY = "MLEAF";

    protected Map<String,Object> composeCallAttributes(final boolean inheritAttributesFromInputVC, final VariantContext vc,
                                                       final AlignmentContext rawContext, final Map<String, AlignmentContext> stratifiedContexts, final RefMetaDataTracker tracker, final ReferenceContext refContext, final List<Integer> alleleCountsofMLE, final boolean bestGuessIsRef,
                                                       final AFCalculationResult AFresult, final List<Allele> allAllelesToUse, final GenotypesContext genotypes,
                                                       final GenotypeLikelihoodsCalculationModel.Model model, final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {
        final HashMap<String, Object> attributes = new HashMap<>();

        final boolean limitedContext = tracker == null || refContext == null || rawContext == null || stratifiedContexts == null;

        // inherit attributes from input vc if requested
        if (inheritAttributesFromInputVC)
            attributes.putAll(vc.getAttributes());
        // if the site was down-sampled, record that fact
        if ( !limitedContext && rawContext.hasPileupBeenDownsampled() )
            attributes.put(DOWNSAMPLED_KEY, true);

        // add the MLE AC and AF annotations
        if ( alleleCountsofMLE.size() > 0 ) {
            attributes.put(MLE_ALLELE_COUNT_KEY, alleleCountsofMLE);
            final ArrayList<Double> MLEfrequencies = calculateMLEAlleleFrequencies(alleleCountsofMLE, genotypes);
            attributes.put(MLE_ALLELE_FREQUENCY_KEY, MLEfrequencies);
        }

        if ( configuration.genotypeArgs.ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED )
            attributes.put(NUMBER_OF_DISCOVERED_ALLELES_KEY, vc.getAlternateAlleles().size());


        return attributes;
    }

    private ArrayList<Double> calculateMLEAlleleFrequencies(List<Integer> alleleCountsofMLE, GenotypesContext genotypes) {
        int AN = 0;
        for (final Genotype g : genotypes)
            for (final Allele a : g.getAlleles())
                if (!a.isNoCall()) AN++;

        final ArrayList<Double> MLEfrequencies = new ArrayList<Double>(alleleCountsofMLE.size());
        // the MLEAC is allowed to be larger than the AN (e.g. in the case of all PLs being 0, the GT is ./. but the exact model may arbitrarily choose an AC>1)
        for (final int AC : alleleCountsofMLE )
            MLEfrequencies.add(Math.min(1.0, (double)AC / (double)AN));
        return MLEfrequencies;
    }

    /**
     * Changes the annotation engine for this genotyping-engine.
     *
     * @param annotationEngine the new annotation engine (can be {@code null}).
     */
    public void setAnnotationEngine(final VariantAnnotatorEngine annotationEngine) {
        this.annotationEngine = annotationEngine;
    }

}
