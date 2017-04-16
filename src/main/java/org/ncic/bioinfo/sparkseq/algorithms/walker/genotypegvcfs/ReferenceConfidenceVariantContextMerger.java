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

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.CommonInfo;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.ncic.bioinfo.sparkseq.algorithms.data.basic.Pair;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GATKVariantContextUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLoc;
import org.ncic.bioinfo.sparkseq.algorithms.utils.MathUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.Utils;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.GenotypeLikelihoodCalculators;
import org.ncic.bioinfo.sparkseq.exceptions.UserException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

/**
 * Author: wbc
 */
public class ReferenceConfidenceVariantContextMerger {

    private static Comparable combineAnnotationValues( final List<Comparable> array ) {
        return MathUtils.median(array); // right now we take the median but other options could be explored
    }

    /**
     * Merges VariantContexts from gVCFs into a single hybrid.
     * Assumes that none of the input records are filtered.
     *
     * @param VCs     collection of unsorted genomic VCs
     * @param loc     the current location
     * @param refBase the reference allele to use if all contexts in the VC are spanning (i.e. don't start at the location in loc); if null, we'll return null in this case
     * @param removeNonRefSymbolicAllele if true, remove the <NON_REF> allele from the merged VC
     * @return new VariantContext representing the merge of all VCs or null if it not relevant
     */
    public static VariantContext merge(final List<VariantContext> VCs, final GenomeLoc loc, final Byte refBase, final boolean removeNonRefSymbolicAllele) {
        // this can happen if e.g. you are using a dbSNP file that spans a region with no gVCFs
        if ( VCs == null || VCs.size() == 0 )
            return null;

        // establish the baseline info (sometimes from the first VC)
        final VariantContext first = VCs.get(0);
        final String name = first.getSource();

        // ref allele
        final Allele refAllele = determineReferenceAlleleGivenReferenceBase(VCs, loc, refBase);
        if ( refAllele == null )
            return null;

        // FinalAlleleSet contains the alleles of the new resulting VC
        // Using linked set in order to guarantee a stable order
        final LinkedHashSet<Allele> finalAlleleSet = new LinkedHashSet<>(10);
        // Reference goes first
        finalAlleleSet.add(refAllele);

        final Map<String, Object> attributes = new LinkedHashMap<>();
        final Set<String> rsIDs = new LinkedHashSet<>(1); // most of the time there's one id
        int depth = 0;
        final Map<String, List<Comparable>> annotationMap = new LinkedHashMap<>();
        final GenotypesContext genotypes = GenotypesContext.create();

        final int variantContextCount = VCs.size();
        // In this list we hold the mapping of each variant context alleles.
        final List<Pair<VariantContext,List<Allele>>> vcAndNewAllelePairs = new ArrayList<>(variantContextCount);
        // cycle through and add info from the other VCs
        for ( final VariantContext vc : VCs ) {

            // if this context doesn't start at the current location then it must be a spanning event (deletion or ref block)
            final boolean isSpanningEvent = loc.getStart() != vc.getStart();

            vcAndNewAllelePairs.add(new Pair<>(vc,isSpanningEvent ? replaceWithNoCalls(vc.getAlleles())
                    : remapAlleles(vc.getAlleles(), refAllele, finalAlleleSet)));
        }

        // Add <NON_REF> to the end if at all required in in the output.
        if (!removeNonRefSymbolicAllele) finalAlleleSet.add(GATKVariantContextUtils.NON_REF_SYMBOLIC_ALLELE);

        final List<Allele> allelesList = new ArrayList<>(finalAlleleSet);

        for ( final Pair<VariantContext,List<Allele>> pair : vcAndNewAllelePairs ) {
            final VariantContext vc = pair.getFirst();
            final List<Allele> remappedAlleles = pair.getSecond();

            mergeRefConfidenceGenotypes(genotypes, vc, remappedAlleles, allelesList);

            // special case DP (add it up) for all events
            if ( vc.hasAttribute(VCFConstants.DEPTH_KEY) ) {
                depth += vc.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0);
            } else { // handle the gVCF case from the HaplotypeCaller
                for( final Genotype gt : vc.getGenotypes() ) {
                    depth += (gt.hasExtendedAttribute("MIN_DP") ? (Integer)gt.getAnyAttribute("MIN_DP") : (gt.hasDP() ? gt.getDP() : 0));
                }
            }

            if ( loc.getStart() != vc.getStart() )
                continue;

            // special case ID (just preserve it)
            if ( vc.hasID() ) rsIDs.add(vc.getID());

            // add attributes
            addReferenceConfidenceAttributes(vc.getAttributes(), annotationMap);
        }

        // when combining annotations use the median value from all input VCs which had annotations provided
        for ( final Map.Entry<String, List<Comparable>> p : annotationMap.entrySet() ) {
            if ( ! p.getValue().isEmpty() ) {
                attributes.put(p.getKey(), combineAnnotationValues(p.getValue()));
            }
        }

        if ( depth > 0 )
            attributes.put(VCFConstants.DEPTH_KEY, String.valueOf(depth));

        // remove stale AC and AF based attributes
        removeStaleAttributesAfterMerge(attributes);

        final String ID = rsIDs.isEmpty() ? VCFConstants.EMPTY_ID_FIELD : Utils.join(",", rsIDs);

        final VariantContextBuilder builder = new VariantContextBuilder().source(name).id(ID).alleles(allelesList)
                .chr(loc.getContig()).start(loc.getStart()).computeEndFromAlleles(allelesList, loc.getStart(), loc.getStart())
                .genotypes(genotypes).unfiltered().attributes(new TreeMap<>(attributes)).log10PError(CommonInfo.NO_LOG10_PERROR);  // we will need to re-genotype later

        return builder.make();
    }

    /**
     * Determines the ref allele given the provided reference base at this position
     *
     * @param VCs     collection of unsorted genomic VCs
     * @param loc     the current location
     * @param refBase the reference allele to use if all contexts in the VC are spanning
     * @return new Allele or null if no reference allele/base is available
     */
    private static Allele determineReferenceAlleleGivenReferenceBase(final List<VariantContext> VCs, final GenomeLoc loc, final Byte refBase) {
        final Allele refAllele = GATKVariantContextUtils.determineReferenceAllele(VCs, loc);
        if ( refAllele == null )
            return ( refBase == null ? null : Allele.create(refBase, true) );
        return refAllele;
    }

    public static final String MLE_ALLELE_COUNT_KEY = "MLEAC";
    public static final String MLE_ALLELE_FREQUENCY_KEY = "MLEAF";

    /**
     * Remove the stale attributes from the merged set
     *
     * @param attributes the attribute map
     */
    private static void removeStaleAttributesAfterMerge(final Map<String, Object> attributes) {
        attributes.remove(VCFConstants.ALLELE_COUNT_KEY);
        attributes.remove(VCFConstants.ALLELE_FREQUENCY_KEY);
        attributes.remove(VCFConstants.ALLELE_NUMBER_KEY);
        attributes.remove(MLE_ALLELE_COUNT_KEY);
        attributes.remove(MLE_ALLELE_FREQUENCY_KEY);
        attributes.remove(VCFConstants.END_KEY);
    }

    /**
     * Adds attributes to the global map from the new context in a sophisticated manner
     *
     * @param myAttributes               attributes to add from
     * @param annotationMap              map of annotations for combining later
     */
    private static void addReferenceConfidenceAttributes(final Map<String, Object> myAttributes,
                                                         final Map<String, List<Comparable>> annotationMap) {
        for ( final Map.Entry<String, Object> p : myAttributes.entrySet() ) {
            final String key = p.getKey();
            final Object value = p.getValue();

            // add the annotation values to a list for combining later
            List<Comparable> values = annotationMap.get(key);
            if( values == null ) {
                values = new ArrayList<>();
                annotationMap.put(key, values);
            }
            try {
                final String stringValue = value.toString();
                // Branch to avoid unintentional, implicit type conversions that occur with the ? operator.
                if (stringValue.contains("."))
                    values.add(Double.parseDouble(stringValue));
                else
                    values.add(Integer.parseInt(stringValue));
            } catch (final NumberFormatException e) {
                // nothing to do
            }
        }
    }

    /**
     * This method does a couple of things:
     * <ul><li>
     *     remaps the vc alleles considering the differences between the final reference allele and its own reference,</li>
     * <li>
     *     collects alternative alleles present in variant context and add them to the {@code finalAlleles} set.
     * </li></ul>
     *
     * @param vcAlleles the variant context allele list.
     * @param refAllele final reference allele.
     * @param finalAlleles where to add the final set of non-ref called alleles.
     * @return never {@code null}
     */
    //TODO as part of a larger refactoring effort {@link #remapAlleles} can be merged with {@link GATKVariantContextUtils#remapAlleles}.
    private static List<Allele> remapAlleles(final List<Allele> vcAlleles, final Allele refAllele, final LinkedHashSet<Allele> finalAlleles) {
        final Allele vcRef = vcAlleles.get(0);
        if (!vcRef.isReference()) throw new IllegalStateException("the first allele of the vc allele list must be reference");
        final byte[] refBases = refAllele.getBases();
        final int extraBaseCount = refBases.length - vcRef.getBases().length;
        if (extraBaseCount < 0) throw new IllegalStateException("the wrong reference was selected");
        final List<Allele> result = new ArrayList<>(vcAlleles.size());

        for (final Allele a : vcAlleles) {
            if (a.isReference()) {
                result.add(refAllele);
            } else if (a.isSymbolic()) {
                result.add(a);
                // we always skip <NON_REF> when adding to finalAlleles this is done outside if applies.
                if (!a.equals(GATKVariantContextUtils.NON_REF_SYMBOLIC_ALLELE))
                    finalAlleles.add(a);
            } else if (a.isCalled()) {
                final Allele newAllele;
                if (extraBaseCount > 0) {
                    final byte[] oldBases = a.getBases();
                    final byte[] newBases = Arrays.copyOf(oldBases,oldBases.length + extraBaseCount);
                    System.arraycopy(refBases,refBases.length - extraBaseCount,newBases,oldBases.length,extraBaseCount);
                    newAllele = Allele.create(newBases,false);
                } else
                    newAllele = a;
                result.add(newAllele);
                finalAlleles.add(newAllele);
            } else { // NO_CALL and strange miscellanea
                result.add(a);
            }
        }
        return result;
    }

    /**
     * Replaces any alleles in the list with NO CALLS, except for the generic ALT allele
     *
     * @param alleles list of alleles to replace
     * @return non-null list of alleles
     */
    private static List<Allele> replaceWithNoCalls(final List<Allele> alleles) {
        if ( alleles == null ) throw new IllegalArgumentException("list of alleles cannot be null");

        final List<Allele> result = new ArrayList<>(alleles.size());
        for ( final Allele allele : alleles )
            result.add(allele.equals(GATKVariantContextUtils.NON_REF_SYMBOLIC_ALLELE) ? allele : Allele.NO_CALL);
        return result;
    }

    /**
     * Merge into the context a new genotype represented by the given VariantContext for the provided list of target alleles.
     * This method assumes that none of the alleles in the VC overlaps with any of the alleles in the set.
     *
     * @param mergedGenotypes   the genotypes context to add to
     * @param VC                the Variant Context for the sample
     * @param remappedAlleles   the list of remapped alleles for the sample
     * @param targetAlleles     the list of target alleles
     */
    private static void mergeRefConfidenceGenotypes(final GenotypesContext mergedGenotypes,
                                                    final VariantContext VC,
                                                    final List<Allele> remappedAlleles,
                                                    final List<Allele> targetAlleles) {
        final int maximumPloidy = VC.getMaxPloidy(GATKVariantContextUtils.DEFAULT_PLOIDY);
        // the map is different depending on the ploidy, so in order to keep this method flexible (mixed ploidies)
        // we need to get a map done (lazily inside the loop) for each ploidy, up to the maximum possible.
        final int[][] genotypeIndexMapsByPloidy = new int[maximumPloidy + 1][];
        final int maximumAlleleCount = Math.max(remappedAlleles.size(),targetAlleles.size());
        final int[] indexesOfRelevantAlleles = getIndexesOfRelevantAlleles(remappedAlleles, targetAlleles, VC.getStart());

        for ( final Genotype g : VC.getGenotypes() ) {
            final String name = g.getSampleName();
            if ( mergedGenotypes.containsSample(name) )
                continue;
            final int ploidy = g.getPloidy();
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(g).alleles(GATKVariantContextUtils.noCallAlleles(g.getPloidy()));
            if (g.hasPL()) {
                // lazy initialization of the genotype index map by ploidy.
                final int[] genotypeIndexMapByPloidy = genotypeIndexMapsByPloidy[ploidy] == null
                        ? GenotypeLikelihoodCalculators.getInstance(ploidy, maximumAlleleCount).genotypeIndexMap(indexesOfRelevantAlleles)
                        : genotypeIndexMapsByPloidy[ploidy];
                final int[] PLs = generatePL(g, genotypeIndexMapByPloidy);
                final int[] AD = g.hasAD() ? generateAD(g.getAD(), indexesOfRelevantAlleles) : null;
                genotypeBuilder.PL(PLs).AD(AD).noGQ();
            }
            mergedGenotypes.add(genotypeBuilder.make());
        }
    }

    /**
     * Composes a new likelihood array given the original genotype and the genotype index map.
     *
     * @param g the original genotype.
     * @param genotypeIndexMapByPloidy genotype index map. The ith element indicates what genotype in {@code g} corresponds
     *                                 to the ith genotype in the return likelihoods array.
     *
     * @throws NullPointerException if {@code g} or {@code genotypeIndexMapByPloidy} is {@code null}, or if {@code g}
     *    does not contain likelihoods.
     * @throws IndexOutOfBoundsException if {@code genotypeIndexMapByPloidy} contain non valid
     *  genotype indices given the likelihood array in {@code g}.
     *
     * @return never {@code null} but an array of exactly {@code genotypeIndexMapByPloidy.length} positions.
     */
    private static int[] generatePL(final Genotype g, final int[] genotypeIndexMapByPloidy) {
        final int[] PLs = new int[genotypeIndexMapByPloidy.length];
        final int[] oldPLs = g.getPL();
        for (int i = 0; i < PLs.length; i++)
            PLs[i] = oldPLs[genotypeIndexMapByPloidy[i]];
        return PLs;
    }

    /**
     * Determines the allele mapping from myAlleles to the targetAlleles, substituting the generic "<ALT>" as appropriate.
     * If the myAlleles set does not contain "<ALT>" as an allele, it throws an exception.
     *
     * @param remappedAlleles   the list of alleles to evaluate
     * @param targetAlleles     the target list of alleles
     * @param position          position to use for error messages
     * @return non-null array of ints representing indexes
     */
    protected static int[] getIndexesOfRelevantAlleles(final List<Allele> remappedAlleles, final List<Allele> targetAlleles, final int position) {

        if ( remappedAlleles == null || remappedAlleles.size() == 0 ) throw new IllegalArgumentException("The list of input alleles must not be null or empty");
        if ( targetAlleles == null || targetAlleles.size() == 0 ) throw new IllegalArgumentException("The list of target alleles must not be null or empty");

        if ( !remappedAlleles.contains(GATKVariantContextUtils.NON_REF_SYMBOLIC_ALLELE) )
            throw new UserException("The list of input alleles must contain " + GATKVariantContextUtils.NON_REF_SYMBOLIC_ALLELE + " as an allele but that is not the case at position " + position + "; please use the Haplotype Caller with gVCF output to generate appropriate records");
        final int indexOfGenericAlt = remappedAlleles.indexOf(GATKVariantContextUtils.NON_REF_SYMBOLIC_ALLELE);

        final int[] indexMapping = new int[targetAlleles.size()];

        // the reference alleles always match up (even if they don't appear to)
        indexMapping[0] = 0;

        // create the index mapping, using the <ALT> allele whenever such a mapping doesn't exist
        for ( int i = 1; i < targetAlleles.size(); i++ ) {
            final int indexOfRemappedAllele = remappedAlleles.indexOf(targetAlleles.get(i));
            indexMapping[i] = indexOfRemappedAllele == -1 ? indexOfGenericAlt : indexOfRemappedAllele;
        }

        return indexMapping;
    }

    /**
     * Generates a new AD array by adding zeros for missing alleles given the set of indexes of the Genotype's current
     * alleles from the original AD.
     *
     * @param originalAD    the original AD to extend
     * @param indexesOfRelevantAlleles the indexes of the original alleles corresponding to the new alleles
     * @return non-null array of new AD values
     */
    protected static int[] generateAD(final int[] originalAD, final int[] indexesOfRelevantAlleles) {
        if ( originalAD == null || indexesOfRelevantAlleles == null ) throw new IllegalArgumentException("The list of input AD values and alleles must not be null");

        final int numADs = indexesOfRelevantAlleles.length;
        final int[] newAD = new int[numADs];

        for ( int i = 0; i < numADs; i++ ) {
            final int oldIndex = indexesOfRelevantAlleles[i];
            if ( oldIndex >= originalAD.length )
                newAD[i] = 0;
            else
                newAD[i] = originalAD[oldIndex];
        }

        return newAD;
    }
}