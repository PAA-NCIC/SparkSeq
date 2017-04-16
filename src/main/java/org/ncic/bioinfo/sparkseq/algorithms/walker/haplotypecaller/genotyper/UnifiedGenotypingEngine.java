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
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.log4j.Logger;
import org.ncic.bioinfo.sparkseq.algorithms.utils.AlignmentContextUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.BaseUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GATKVariantContextUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLocParser;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenotypingGivenAllelesUtils;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefMetaDataTracker;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.ReferenceContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.AlignmentContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.PileupElement;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.ReadBackedPileup;
import org.ncic.bioinfo.sparkseq.algorithms.walker.baserecalibrator.BAQ;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.GenotypingOutputMode;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.OutputMode;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.PerReadAlleleLikelihoodMap;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.SampleList;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.afcalculate.AFCalculationResult;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.afcalculate.AFCalculatorProvider;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.argcollection.UnifiedArgumentCollection;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.model.GeneralPloidyGenotypeLikelihoodsCalculationModel;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.model.GenotypeLikelihoodsCalculationModel;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.model.IndelGenotypeLikelihoodsCalculationModel;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.model.SNPGenotypeLikelihoodsCalculationModel;
import org.ncic.bioinfo.sparkseq.exceptions.UserException;

import java.io.PrintStream;
import java.lang.reflect.Constructor;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Author: wbc
 */
public class UnifiedGenotypingEngine extends GenotypingEngine<UnifiedArgumentCollection> {

    public static final String PL_FOR_ALL_SNP_ALLELES_KEY = "APL";

    private static final int SNP_MODEL = 0;
    private static final int INDEL_MODEL = 1;

    // the model used for calculating genotypes
    private ThreadLocal<Map<String, GenotypeLikelihoodsCalculationModel>> glcm;
    private final List<GenotypeLikelihoodsCalculationModel.Model> modelsToUse = new ArrayList<>(2);


    // the various loggers and writers
    private PrintStream verboseWriter;

    private final boolean BAQEnabledOnCMDLine;

    // ---------------------------------------------------------------------------------------------------------
    //
    // Public interface functions
    //
    // ---------------------------------------------------------------------------------------------------------


    /**
     * Creates a new unified genotyping given the UG configuration parameters and the GA engine.
     *
     * @param configuration the UG configuration.
     * @throws NullPointerException if either {@code configuration} or {@code toolkit} is {@code null}.
     */
    public UnifiedGenotypingEngine(final UnifiedArgumentCollection configuration,
                                   final AFCalculatorProvider afCalculatorProvider,
                                   final SampleList sampleList,
                                   final GenomeLocParser genomeLocParser,
                                   final BAQ.CalculationMode baqMode) {
        this(configuration, sampleList, genomeLocParser, afCalculatorProvider, baqMode);
    }


    /**
     * Creates a new unified genotyping given the UG configuration parameters, the targeted set of samples and
     * a genome location parser.
     *
     * @param configuration      the UG configuration.
     * @param samples            {@inheritDoc}
     * @param baqCalculationMode the BAQ calculation mode.
     * @throws NullPointerException     if any of {@code configuration}, {@code samples} or {@code genomeLocParser} is {@code null}.
     * @throws IllegalArgumentException if {@code baqCalculationMode} is {@code null}.
     */
    public UnifiedGenotypingEngine(final UnifiedArgumentCollection configuration,
                                   final SampleList samples, final GenomeLocParser genomeLocParser, final AFCalculatorProvider afCalculatorProvider,
                                   final BAQ.CalculationMode baqCalculationMode) {

        super(configuration, samples, genomeLocParser, afCalculatorProvider);

        if (baqCalculationMode == null)
            throw new IllegalArgumentException("the BAQ calculation mode cannot be null");

        this.BAQEnabledOnCMDLine = baqCalculationMode != BAQ.CalculationMode.OFF;

        determineGLModelsToUse();

        initializeGenotypeLikelihoodsCalculationModels();
    }

    /**
     * Changes the verbose output writer for this engine.
     *
     * @param writer the new writer; it can be {@code null}.
     */
    public void setVerboseWriter(final PrintStream writer) {
        verboseWriter = writer;
    }

    /**
     * Initialize {@link #glcm}.
     */
    private void initializeGenotypeLikelihoodsCalculationModels() {
        glcm = new ThreadLocal<Map<String, GenotypeLikelihoodsCalculationModel>>() {

            @Override
            protected Map<String, GenotypeLikelihoodsCalculationModel> initialValue() {
                return getGenotypeLikelihoodsCalculationObject(logger, UnifiedGenotypingEngine.this.configuration);
            }
        };
    }

    /**
     * Compute full calls at a given locus. Entry point for engine calls from the UnifiedGenotyper.
     *
     * @param tracker    the meta data tracker
     * @param refContext the reference base
     * @param rawContext contextual information around the locus
     * @return the VariantCallContext object
     */
    public List<VariantCallContext> calculateLikelihoodsAndGenotypes(final RefMetaDataTracker tracker,
                                                                     final ReferenceContext refContext,
                                                                     final AlignmentContext rawContext) {
        final List<VariantCallContext> results = new ArrayList<>(2);

        final List<GenotypeLikelihoodsCalculationModel.Model> models = getGLModelsToUse(tracker, rawContext);

        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap = new HashMap<>();

        final VariantCallContext defaultResult = configuration.outputMode == OutputMode.EMIT_ALL_SITES
                && configuration.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES
                ? emptyCallContext(tracker, refContext, rawContext)
                : null;

        if (models.isEmpty())
            results.add(defaultResult);
        else {
            for (final GenotypeLikelihoodsCalculationModel.Model model : models) {
                perReadAlleleLikelihoodMap.clear();
                final Map<String, AlignmentContext> stratifiedContexts = getFilteredAndStratifiedContexts(refContext, rawContext, model);
                if (stratifiedContexts == null)
                    results.add(defaultResult);
                else {
                    final VariantContext vc = calculateLikelihoods(tracker, refContext, stratifiedContexts, AlignmentContextUtils.ReadOrientation.COMPLETE, null, true, model, perReadAlleleLikelihoodMap);
                    if (vc != null)
                        results.add(calculateGenotypes(tracker, refContext, rawContext, stratifiedContexts, vc, model, true, perReadAlleleLikelihoodMap));
// todo - uncomment if we want to also emit a null ref call (with no QUAL) if there's no evidence for REF and if EMIT_ALL_SITES is set
//                    else if (UAC.OutputMode == OUTPUT_MODE.EMIT_ALL_SITES)
//                        results.add(generateEmptyContext(tracker, refContext, null, rawContext));

                }
            }
        }
        return results;
    }

    private static Map<String, GenotypeLikelihoodsCalculationModel> getGenotypeLikelihoodsCalculationObject(Logger logger, UnifiedArgumentCollection UAC) {

        final Map<String, GenotypeLikelihoodsCalculationModel> glcm = new HashMap<>();
        final List<Class<? extends GenotypeLikelihoodsCalculationModel>> glmClasses = new ArrayList<>();
        glmClasses.add(GeneralPloidyGenotypeLikelihoodsCalculationModel.class);

        glmClasses.add(IndelGenotypeLikelihoodsCalculationModel.class);
        glmClasses.add(SNPGenotypeLikelihoodsCalculationModel.class);
        for (final Class<? extends GenotypeLikelihoodsCalculationModel> glmClass : glmClasses) {
            final String key = glmClass.getSimpleName().replaceAll("GenotypeLikelihoodsCalculationModel", "").toUpperCase();
            try {
                final Object args[] = new Object[]{UAC, logger};
                final Constructor c = glmClass.getDeclaredConstructor(UnifiedArgumentCollection.class, Logger.class);
                glcm.put(key, (GenotypeLikelihoodsCalculationModel) c.newInstance(args));
            } catch (Exception e) {
                throw new UserException("The likelihoods model provided for the -glm argument (" + UAC.GLmodel + ") is not a valid option: " + e.getMessage());
            }
        }

        return glcm;
    }

    /**
     * Compute GLs at a given locus. Entry point for engine calls from UGCalcLikelihoods.
     *
     * @param tracker                    the meta data tracker
     * @param refContext                 the reference base
     * @param rawContext                 contextual information around the locus
     * @param perReadAlleleLikelihoodMap Map to store per-sample, per-read, per-allele likelihoods (only used for indels)
     * @return the VariantContext object
     */
    public VariantContext calculateLikelihoods(final RefMetaDataTracker tracker,
                                               final ReferenceContext refContext,
                                               final AlignmentContext rawContext,
                                               final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {
        final List<GenotypeLikelihoodsCalculationModel.Model> models = getGLModelsToUse(tracker, rawContext);
        if (models.isEmpty()) {
            return null;
        }

        for (final GenotypeLikelihoodsCalculationModel.Model model : models) {
            final Map<String, AlignmentContext> stratifiedContexts = getFilteredAndStratifiedContexts(refContext, rawContext, model);
            // return the first valid one we encounter
            if (stratifiedContexts != null)
                return calculateLikelihoods(tracker, refContext, stratifiedContexts, AlignmentContextUtils.ReadOrientation.COMPLETE, null, true, model, perReadAlleleLikelihoodMap);

        }

        return null;
    }

    /**
     * Compute genotypes at a given locus. Entry point for engine calls from UGCallVariants.
     *
     * @param tracker    the meta data tracker
     * @param refContext the reference base
     * @param rawContext contextual information around the locus
     * @param vc         the GL-annotated variant context
     * @return the VariantCallContext object
     */
    public VariantCallContext calculateGenotypes(final RefMetaDataTracker tracker,
                                                 final ReferenceContext refContext,
                                                 final AlignmentContext rawContext,
                                                 final VariantContext vc) {
        final List<GenotypeLikelihoodsCalculationModel.Model> models = getGLModelsToUse(tracker, rawContext);
        if (models.isEmpty()) {
            return null;
        }

        // return the first one
        final GenotypeLikelihoodsCalculationModel.Model model = models.get(0);
        final Map<String, AlignmentContext> stratifiedContexts = getFilteredAndStratifiedContexts(refContext, rawContext, model);
        return calculateGenotypes(tracker, refContext, rawContext, stratifiedContexts, vc, model, null);
    }

    /**
     * Compute genotypes at a given locus.
     *
     * @param vc the GL-annotated variant context
     * @return the VariantCallContext object
     */
    public VariantCallContext calculateGenotypes(VariantContext vc) {
        return calculateGenotypes(null, null, null, null, vc, GenotypeLikelihoodsCalculationModel.Model.valueOf("SNP"), null);
    }


    // ---------------------------------------------------------------------------------------------------------
    //
    // Private implementation helpers
    //
    // ---------------------------------------------------------------------------------------------------------

    // private method called by both UnifiedGenotyper and UGCalcLikelihoods entry points into the engine
    private VariantContext calculateLikelihoods(final RefMetaDataTracker tracker,
                                                final ReferenceContext refContext,
                                                final Map<String, AlignmentContext> stratifiedContexts,
                                                final AlignmentContextUtils.ReadOrientation type,
                                                final List<Allele> alternateAllelesToUse,
                                                final boolean useBAQedPileup,
                                                final GenotypeLikelihoodsCalculationModel.Model model,
                                                final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {

        return glcm.get().get(model.name()).getLikelihoods(tracker, refContext, stratifiedContexts, type, alternateAllelesToUse, useBAQedPileup && BAQEnabledOnCMDLine,
                genomeLocParser, perReadAlleleLikelihoodMap);
    }


    public VariantCallContext calculateGenotypes(final VariantContext vc, final GenotypeLikelihoodsCalculationModel.Model model) {
        return calculateGenotypes(null, null, null, null, vc, model, null);
    }

    public VariantCallContext calculateGenotypes(final RefMetaDataTracker tracker,
                                                 final ReferenceContext refContext,
                                                 final AlignmentContext rawContext,
                                                 final Map<String, AlignmentContext> stratifiedContexts,
                                                 final VariantContext vc,
                                                 final GenotypeLikelihoodsCalculationModel.Model model,
                                                 final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {
        return calculateGenotypes(tracker, refContext, rawContext, stratifiedContexts, vc, model, false, perReadAlleleLikelihoodMap);
    }

    @Override
    protected boolean forceKeepAllele(final Allele allele) {
        return configuration.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES || configuration.annotateAllSitesWithPLs;
    }

    @Override
    public VariantCallContext calculateGenotypes(final RefMetaDataTracker tracker, final ReferenceContext refContext,
                                                 final AlignmentContext rawContext, Map<String, AlignmentContext> stratifiedContexts,
                                                 final VariantContext vc, final GenotypeLikelihoodsCalculationModel.Model model,
                                                 final boolean inheritAttributesFromInputVC,
                                                 final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {
        boolean limitedContext = tracker == null || refContext == null || rawContext == null || stratifiedContexts == null;
        final VariantCallContext result = super.calculateGenotypes(tracker, refContext, rawContext, stratifiedContexts, vc, model, inheritAttributesFromInputVC, perReadAlleleLikelihoodMap);
        if (verboseWriter != null && !limitedContext)
            printVerboseData(refContext.getLocus().toString(), vc, model);
        return result;
    }

    @Override
    protected String callSourceString() {
        return "UG_call";
    }

    @Override
    protected boolean forceSiteEmission() {
        return configuration.outputMode == OutputMode.EMIT_ALL_SITES;
    }


    @Override
    protected Map<String, Object> composeCallAttributes(final boolean inheritAttributesFromInputVC, final VariantContext vc,
                                                        final AlignmentContext rawContext, final Map<String, AlignmentContext> stratifiedContexts, final RefMetaDataTracker tracker, final ReferenceContext refContext, final List<Integer> alleleCountsofMLE, final boolean bestGuessIsRef,
                                                        final AFCalculationResult AFresult, final List<Allele> allAllelesToUse, final GenotypesContext genotypes,
                                                        final GenotypeLikelihoodsCalculationModel.Model model, final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {
        final Map<String, Object> result = super.composeCallAttributes(inheritAttributesFromInputVC, vc, rawContext, stratifiedContexts, tracker, refContext, alleleCountsofMLE, bestGuessIsRef,
                AFresult, allAllelesToUse, genotypes, model, perReadAlleleLikelihoodMap);

        final boolean limitedContext = tracker == null || refContext == null || rawContext == null || stratifiedContexts == null;

        if (configuration.COMPUTE_SLOD && !limitedContext && !bestGuessIsRef) {
            final double strandScore = calculateSLOD(stratifiedContexts, tracker, refContext, AFresult, allAllelesToUse, model, perReadAlleleLikelihoodMap);
            if (!Double.isNaN(strandScore))
                result.put("SB", strandScore);
        }
        return result;
    }

    private double calculateSLOD(final Map<String, AlignmentContext> stratifiedContexts,
                                 final RefMetaDataTracker tracker,
                                 final ReferenceContext refContext, final AFCalculationResult AFresult,
                                 final List<Allele> allAllelesToUse,
                                 final GenotypeLikelihoodsCalculationModel.Model model,
                                 final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {
        // the overall lod
        //double overallLog10PofNull = AFresult.log10AlleleFrequencyPosteriors[0];
        final double overallLog10PofF = AFresult.getLog10LikelihoodOfAFGT0();
        //if ( DEBUG_SLOD ) System.out.println("overallLog10PofF=" + overallLog10PofF);

        // the forward lod
        final AFCalculationResult forwardAFresult = getDirectionalAfCalcResult(AlignmentContextUtils.ReadOrientation.FORWARD, stratifiedContexts, tracker, refContext, allAllelesToUse, model, perReadAlleleLikelihoodMap);
        final double forwardLog10PofNull = forwardAFresult.getLog10LikelihoodOfAFEq0();
        final double forwardLog10PofF = forwardAFresult.getLog10LikelihoodOfAFGT0();
        //if ( DEBUG_SLOD ) System.out.println("forwardLog10PofNull=" + forwardLog10PofNull + ", forwardLog10PofF=" + forwardLog10PofF);

        // the reverse lod
        final AFCalculationResult reverseAFresult = getDirectionalAfCalcResult(AlignmentContextUtils.ReadOrientation.REVERSE, stratifiedContexts, tracker, refContext, allAllelesToUse, model, perReadAlleleLikelihoodMap);
        final double reverseLog10PofNull = reverseAFresult.getLog10LikelihoodOfAFEq0();
        final double reverseLog10PofF = reverseAFresult.getLog10LikelihoodOfAFGT0();
        //if ( DEBUG_SLOD ) System.out.println("reverseLog10PofNull=" + reverseLog10PofNull + ", reverseLog10PofF=" + reverseLog10PofF);

        final double forwardLod = forwardLog10PofF + reverseLog10PofNull - overallLog10PofF;
        final double reverseLod = reverseLog10PofF + forwardLog10PofNull - overallLog10PofF;
        //if ( DEBUG_SLOD ) System.out.println("forward lod=" + forwardLod + ", reverse lod=" + reverseLod);

        // strand score is max bias between forward and reverse strands
        double strandScore = Math.max(forwardLod, reverseLod);
        // rescale by a factor of 10
        strandScore *= 10.0;
        //logger.debug(String.format("SLOD=%f", strandScore));
        return strandScore;
    }

    private AFCalculationResult getDirectionalAfCalcResult(final AlignmentContextUtils.ReadOrientation orientation,
                                                           final Map<String, AlignmentContext> stratifiedContexts,
                                                           final RefMetaDataTracker tracker,
                                                           final ReferenceContext refContext, List<Allele> allAllelesToUse,
                                                           final GenotypeLikelihoodsCalculationModel.Model model,
                                                           final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {
        final VariantContext vc = calculateLikelihoods(tracker, refContext, stratifiedContexts, orientation,
                allAllelesToUse, false, model, perReadAlleleLikelihoodMap);
        final int defaultPloidy = configuration.genotypeArgs.samplePloidy;
        final int maxAltAlleles = configuration.genotypeArgs.MAX_ALTERNATE_ALLELES;
        return afCalculatorProvider.getInstance(vc, defaultPloidy, maxAltAlleles).getLog10PNonRef(vc, defaultPloidy, maxAltAlleles, getAlleleFrequencyPriors(vc, defaultPloidy, model));
    }

    private Map<String, AlignmentContext> getFilteredAndStratifiedContexts(final ReferenceContext refContext,
                                                                           final AlignmentContext rawContext,
                                                                           final GenotypeLikelihoodsCalculationModel.Model model) {

        if (!BaseUtils.isRegularBase(refContext.getBase()))
            return null;


        switch (model) {
            case INDEL:
            case GENERALPLOIDYINDEL:

                final ReadBackedPileup pileup = rawContext.getBasePileup().getMappingFilteredPileup(configuration.MIN_BASE_QUALTY_SCORE);
                // don't call when there is no coverage
                if (pileup.getNumberOfElements() == 0 && configuration.outputMode != OutputMode.EMIT_ALL_SITES)
                    return null;

                // stratify the AlignmentContext and cut by sample
                return AlignmentContextUtils.splitContextBySampleName(pileup);
            case SNP:
            case GENERALPLOIDYSNP:

                if (!(configuration.outputMode == OutputMode.EMIT_ALL_SITES && configuration.genotypingOutputMode != GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES)) {
                    int numDeletions = 0;
                    for (final PileupElement p : rawContext.getBasePileup()) {
                        if (p.isDeletion())
                            numDeletions++;
                    }
                    if (((double) numDeletions) / ((double) rawContext.getBasePileup().depthOfCoverage()) > configuration.MAX_DELETION_FRACTION) {
                        return null;
                    }
                }
                // stratify the AlignmentContext and cut by sample
                return AlignmentContextUtils.splitContextBySampleName(rawContext.getBasePileup());
            default:
                throw new IllegalStateException("unexpected model: " + model);
        }
    }

    protected void printVerboseData(final String pos, final VariantContext vc, final GenotypeLikelihoodsCalculationModel.Model model) {
        Allele refAllele = null, altAllele = null;
        for (Allele allele : vc.getAlleles()) {
            if (allele.isReference())
                refAllele = allele;
            else
                altAllele = allele;
        }

        for (int i = 0; i <= numberOfGenomes; i++) {
            StringBuilder AFline = new StringBuilder("AFINFO\t");
            AFline.append(pos);
            AFline.append('\t');
            AFline.append(refAllele);
            AFline.append('\t');
            if (altAllele != null)
                AFline.append(altAllele);
            else
                AFline.append("N/A");
            AFline.append('\t');
            AFline.append(i).append('/').append(numberOfGenomes).append('\t');
            AFline.append(String.format("%.2f\t", ((float) i) / numberOfGenomes));
            AFline.append(String.format("%.8f\t", getAlleleFrequencyPriors(vc, configuration.genotypeArgs.samplePloidy, model)[i]));
            verboseWriter.println(AFline.toString());
        }

        verboseWriter.println("Qscore = " + vc.getLog10PError());
        verboseWriter.println();
    }

    /**
     * Determine the models to be use for genotype likelihood calculation.
     * <p>
     * <p>
     * Whether to use the general ones or the concrete diploid ones need to depend on what the user has said
     * in its parameter glm. Iff he selected GENERALPLOIDYINDEL or GENERALPLOIDYSNP is the general set, otherwise
     * </p>
     * <p>
     * the standard set (SNP and INDEL).
     * <p>
     * Even if the user did not select to use both models, GGA force the inclusion of both: snp and indels.
     * </p>
     * <p>
     * Also, we must fail
     * </p>
     * <p>
     * The models are added to the field {@link #modelsToUse}.
     */
    private void determineGLModelsToUse() {

        modelsToUse.clear();

        boolean useGeneral = false;
        boolean useSNP = false;
        boolean useINDEL = false;

        switch (configuration.GLmodel) {
            case BOTH:
                useSNP = useINDEL = true;
                break;
            case SNP:
                useSNP = true;
                break;
            case INDEL:
                useINDEL = true;
                break;
            case GENERALPLOIDYINDEL:
                useINDEL = useGeneral = true;
                break;
            case GENERALPLOIDYSNP:
                useSNP = useGeneral = true;
                break;
            default: //Paranoia
                throw new IllegalStateException("unexpected genotype likelihood model " + configuration.GLmodel);
        }

        // Force the use of both models in GGA:
        if (configuration.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES)
            useSNP = useINDEL = true;

        // The non-general models only support Diploid so need to go to general if not the default_ploidy == 2.
        if (configuration.genotypeArgs.samplePloidy != GATKVariantContextUtils.DEFAULT_PLOIDY)
            useGeneral = true;

        // If annotateAllSitesWithPLs requested , SNP model must be used.
        if (!useSNP && configuration.annotateAllSitesWithPLs)
            throw new UserException.BadArgumentValue("glm", "Invalid genotype likelihood model specification: " +
                    "only diploid SNP model can be used in conjunction with option allSitePLs");

        // Finally add the relevant model
        if (useSNP)
            modelsToUse.add(useGeneral ? GenotypeLikelihoodsCalculationModel.Model.GENERALPLOIDYSNP :
                    GenotypeLikelihoodsCalculationModel.Model.SNP);
        if (useINDEL)
            modelsToUse.add(useGeneral ? GenotypeLikelihoodsCalculationModel.Model.GENERALPLOIDYINDEL :
                    GenotypeLikelihoodsCalculationModel.Model.INDEL);
    }

    // decide whether we are currently processing SNPs, indels, neither, or both
    private List<GenotypeLikelihoodsCalculationModel.Model> getGLModelsToUse(final RefMetaDataTracker tracker,
                                                                             final AlignmentContext rawContext) {
        if (configuration.genotypingOutputMode != GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES)
            return modelsToUse;

        if (modelsToUse.size() != 2)
            throw new IllegalStateException("GGA mode assumes that we have initialized both the SNP and indel models but found " + modelsToUse);

        // if we're genotyping given alleles then we need to choose the model corresponding to the variant type requested
        final VariantContext vcInput = GenotypingGivenAllelesUtils.composeGivenAllelesVariantContextFromRod(tracker, rawContext.getLocation(), false, logger, configuration.alleles);

        if (vcInput == null) {
            return Collections.emptyList(); // no work to be done
        } else if (vcInput.isSNP()) {
            return Collections.singletonList(modelsToUse.get(SNP_MODEL));
        } else if (vcInput.isIndel() || vcInput.isMixed()) {
            return Collections.singletonList(modelsToUse.get(INDEL_MODEL));
        } else {
            return Collections.emptyList(); // No support for other types yet
        }
    }

}
