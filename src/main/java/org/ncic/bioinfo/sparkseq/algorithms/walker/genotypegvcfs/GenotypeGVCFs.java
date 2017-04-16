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
package org.ncic.bioinfo.sparkseq.algorithms.walker.genotypegvcfs;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefMetaDataTracker;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.ReferenceContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.SamContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.GATKFeature;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RODContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RODNames;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RodBinding;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.Tags;
import org.ncic.bioinfo.sparkseq.algorithms.engine.RODWalker;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GATKVariantContextUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLoc;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLocParser;
import org.ncic.bioinfo.sparkseq.algorithms.walker.baserecalibrator.BAQ;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.HaplotypeCaller;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.IndexedSampleList;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.SampleList;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.afcalculate.GeneralPloidyFailOverAFCalculatorProvider;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.InbreedingCoeff;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.QualByDepth;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.StrandOddsRatio;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.VariantAnnotatorEngine;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.interfaces.AnnotatorCompatible;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.argcollection.DbsnpArgumentCollection;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.argcollection.GenotypeCalculationArgumentCollection;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.argcollection.UnifiedArgumentCollection;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.genotyper.GenotypingEngine;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.genotyper.UnifiedGenotypingEngine;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Author: wbc
 */
public class GenotypeGVCFs extends RODWalker implements AnnotatorCompatible {

    final private List<RodBinding<VariantContext>> variants = new ArrayList<>();

    public boolean INCLUDE_NON_VARIANTS = false;

    public GenotypeCalculationArgumentCollection genotypeArgs = new GenotypeCalculationArgumentCollection();

    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    public RodBinding<VariantContext> getDbsnpRodBinding() {
        return dbsnp.dbsnp;
    }

    protected List<String> infoFieldAnnotations = Arrays.asList(new String[]{"ChromosomeCounts", "FisherStrand",
            "GenotypeSummaries", "InbreedingCoeff", "QualByDepth", "StrandOddsRatio"});

    protected List<String> genotypeAnnotations = new ArrayList<>(Arrays.asList(new String[]{}));

    // the genotyping engine
    private UnifiedGenotypingEngine genotypingEngine;
    // the annotation engine
    private VariantAnnotatorEngine annotationEngine;

    public List<RodBinding<VariantContext>> getCompRodBindings() { return Collections.emptyList(); }
    public RodBinding<VariantContext> getSnpEffRodBinding() { return null; }
    public List<RodBinding<VariantContext>> getResourceRodBindings() { return Collections.emptyList(); }
    public boolean alwaysAppendDbsnpId() { return false; }

    private List<VariantContext> resultVcfRecords = new ArrayList<>();

    public GenotypeGVCFs(GenomeLocParser genomeLocParser,
                         RefContentProvider refContentProvider,
                         SamContentProvider samContentProvider,
                         List<VariantContext> gvcfList) {
        super(genomeLocParser, refContentProvider, samContentProvider,
                Lists.newArrayList(new RODContentProvider(RODNames.VARIANT, null, gvcfList, genomeLocParser)));
    }

    @Override
    protected void initialize() {
        SampleList samples = getSamplesListFromContentProvider();

        variants.add(new RodBinding<>(VariantContext.class, RODNames.VARIANT, "gvcf", "VCF", new Tags()));

        // create the genotyping engine
        genotypingEngine = new UnifiedGenotypingEngine(createUAC(), samples, genomeLocParser,
                GeneralPloidyFailOverAFCalculatorProvider.createThreadSafeProvider(genotypeArgs, logger),
                BAQ.CalculationMode.OFF);
        // create the annotation engine
        annotationEngine = new VariantAnnotatorEngine(infoFieldAnnotations, genotypeAnnotations,
                Collections.<String>emptyList(), this, genomeLocParser);

    }

    private SampleList getSamplesListFromContentProvider() {
        List<SAMReadGroupRecord> readGroupInfos =
                samContentProvider.getSamFileHeader().getReadGroups();
        List<String> samples = new ArrayList<>();
        for (SAMReadGroupRecord readGroup : readGroupInfos) {
            samples.add(readGroup.getSample());
        }
        SampleList sampleList = new IndexedSampleList(samples);
        return sampleList;
    }

    /**
     * Creates a UnifiedArgumentCollection with appropriate values filled in from the arguments in this walker
     *
     * @return a complete UnifiedArgumentCollection
     */
    private UnifiedArgumentCollection createUAC() {
        UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        uac.genotypeArgs = genotypeArgs.clone();
        return uac;
    }

    private List<VariantContext> getVCInTrackerInLocus(String rodName, final RefMetaDataTracker tracker) {
        List<Feature> features = tracker.getValues(rodName);
        List<VariantContext> vcs = new ArrayList<>(features.size());
        for (Feature feature : features) {
            vcs.add((VariantContext) ((GATKFeature.TribbleGATKFeature) feature).getUnderlyingObject());
        }
        return vcs;
    }

    @Override
    protected void map(final ReferenceContext ref,
                       final RefMetaDataTracker tracker) {
        if ( tracker == null ) // RodWalkers can make funky map calls
            return;

        final GenomeLoc loc = ref.getLocus();
        final VariantContext combinedVC = ReferenceConfidenceVariantContextMerger.merge(getVCInTrackerInLocus("variant", tracker), loc, INCLUDE_NON_VARIANTS ? ref.getBase() : null, true);
        if ( combinedVC == null )
            return;
        VariantContext finalRes = regenotypeVC(tracker, ref, combinedVC);
        if(finalRes != null) {
            resultVcfRecords.add(finalRes);
        }
    }

    /**
     * Re-genotype (and re-annotate) a combined genomic VC
     *
     * @param tracker        the ref tracker
     * @param ref            the ref context
     * @param originalVC     the combined genomic VC
     * @return a new VariantContext or null if the site turned monomorphic and we don't want such sites
     */
    protected VariantContext regenotypeVC(final RefMetaDataTracker tracker, final ReferenceContext ref, final VariantContext originalVC) {
        if ( originalVC == null ) throw new IllegalArgumentException("originalVC cannot be null");

        VariantContext result = originalVC;

        // only re-genotype polymorphic sites
        if ( result.isVariant() ) {
            VariantContext regenotypedVC = genotypingEngine.calculateGenotypes(result);
            if ( regenotypedVC == null) {
                if (!INCLUDE_NON_VARIANTS)
                    return null;
            }
            else {
                regenotypedVC = GATKVariantContextUtils.reverseTrimAlleles(regenotypedVC);
                result = addGenotypingAnnotations(originalVC.getAttributes(), regenotypedVC);
            }
        }

        // if it turned monomorphic then we either need to ignore or fix such sites
        boolean createRefGTs = false;
        if ( result.isMonomorphicInSamples() ) {
            if ( !INCLUDE_NON_VARIANTS )
                return null;
            createRefGTs = true;
        }

        // Re-annotate and fix/remove some of the original annotations.
        // Note that the order of these actions matters and is different for polymorphic and monomorphic sites.
        // For polymorphic sites we need to make sure e.g. the SB tag is sent to the annotation engine and then removed later.
        // For monomorphic sites we need to make sure e.g. the hom ref genotypes are created and only then are passed to the annotation engine.
        // We could theoretically make 2 passes to re-create the genotypes, but that gets extremely expensive with large sample sizes.
        if ( createRefGTs ) {
            result = new VariantContextBuilder(result).genotypes(cleanupGenotypeAnnotations(result, true)).make();
            result = annotationEngine.annotateContext(tracker, ref, null, result);
        } else {
            result = annotationEngine.annotateContext(tracker, ref, null, result);
            result = new VariantContextBuilder(result).genotypes(cleanupGenotypeAnnotations(result, false)).make();
        }

        return result;
    }


    public static final String MLE_ALLELE_COUNT_KEY = "MLEAC";
    public static final String MLE_ALLELE_FREQUENCY_KEY = "MLEAF";

    /**
     * Add genotyping-based annotations to the new VC
     *
     * @param originalAttributes the non-null annotations from the original VC
     * @param newVC the new non-null VC
     * @return a non-null VC
     */
    private VariantContext addGenotypingAnnotations(final Map<String, Object> originalAttributes, final VariantContext newVC) {
        // we want to carry forward the attributes from the original VC but make sure to add the MLE-based annotations
        final Map<String, Object> attrs = new HashMap<>(originalAttributes);
        attrs.put(MLE_ALLELE_COUNT_KEY, newVC.getAttribute(MLE_ALLELE_COUNT_KEY));
        attrs.put(MLE_ALLELE_FREQUENCY_KEY, newVC.getAttribute(MLE_ALLELE_FREQUENCY_KEY));
        if (newVC.hasAttribute(GenotypingEngine.NUMBER_OF_DISCOVERED_ALLELES_KEY))
            attrs.put(GenotypingEngine.NUMBER_OF_DISCOVERED_ALLELES_KEY, newVC.getAttribute(GenotypingEngine.NUMBER_OF_DISCOVERED_ALLELES_KEY));

        return new VariantContextBuilder(newVC).attributes(attrs).make();
    }

    /**
     * Cleans up genotype-level annotations that need to be updated.
     * 1. move MIN_DP to DP if present
     * 2. propagate DP to AD if not present
     * 3. remove SB if present
     * 4. change the PGT value from "0|1" to "1|1" for homozygous variant genotypes
     *
     * @param VC            the VariantContext with the Genotypes to fix
     * @param createRefGTs  if true we will also create proper hom ref genotypes since we assume the site is monomorphic
     * @return a new set of Genotypes
     */
    private List<Genotype> cleanupGenotypeAnnotations(final VariantContext VC, final boolean createRefGTs) {
        final GenotypesContext oldGTs = VC.getGenotypes();
        final List<Genotype> recoveredGs = new ArrayList<>(oldGTs.size());
        for ( final Genotype oldGT : oldGTs ) {
            final Map<String, Object> attrs = new HashMap<>(oldGT.getExtendedAttributes());

            final GenotypeBuilder builder = new GenotypeBuilder(oldGT);
            int depth = oldGT.hasDP() ? oldGT.getDP() : 0;

            // move the MIN_DP to DP
            if ( oldGT.hasExtendedAttribute("MIN_DP") ) {
                depth = Integer.parseInt((String)oldGT.getAnyAttribute("MIN_DP"));
                builder.DP(depth);
                attrs.remove("MIN_DP");
            }

            // remove SB
            attrs.remove("SB");

            // update PGT for hom vars
            if ( oldGT.isHomVar() && oldGT.hasExtendedAttribute(HaplotypeCaller.HAPLOTYPE_CALLER_PHASING_GT_KEY) ) {
                attrs.put(HaplotypeCaller.HAPLOTYPE_CALLER_PHASING_GT_KEY, "1|1");
            }

            // create AD if it's not there
            if ( !oldGT.hasAD() && VC.isVariant() ) {
                final int[] AD = new int[VC.getNAlleles()];
                AD[0] = depth;
                builder.AD(AD);
            }

            if ( createRefGTs ) {
                final int ploidy = oldGT.getPloidy();
                final List<Allele> refAlleles = Collections.nCopies(ploidy,VC.getReference());

                //keep 0 depth samples as no-call
                if (depth > 0) {
                    builder.alleles(refAlleles);
                }

                // also, the PLs are technically no longer usable
                builder.noPL();
            }

            recoveredGs.add(builder.noAttributes().attributes(attrs).make());
        }
        return recoveredGs;
    }


    @Override
    protected void onTraversalDone() {

    }

    public List<VariantContext> getResultVcfRecords() {
        return resultVcfRecords;
    }
}
