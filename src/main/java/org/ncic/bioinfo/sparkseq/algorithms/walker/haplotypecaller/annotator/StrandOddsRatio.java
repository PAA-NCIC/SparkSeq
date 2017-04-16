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
package org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator;

import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.AlignmentContext;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.PerReadAlleleLikelihoodMap;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.interfaces.ActiveRegionBasedAnnotation;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.interfaces.StandardAnnotation;

import java.util.*;

/**
 * Strand bias estimated by the Symmetric Odds Ratio test
 *
 * <p>Strand bias is a type of sequencing bias in which one DNA strand is favored over the other, which can result in incorrect evaluation of the amount of evidence observed for one allele vs. the other. The StrandOddsRatio annotation is one of several methods that aims to evaluate whether there is strand bias in the data. It is an updated form of the Fisher Strand Test that is better at taking into account large amounts of data in high coverage situations. It is used to determine if there is strand bias between forward and reverse strands for the reference or alternate allele.</p>
 *
 * <h3>Statistical notes</h3>
 * <p> Odds Ratios in the 2x2 contingency table below are
 *
 * $$ R = \frac{X[0][0] * X[1][1]}{X[0][1] * X[1][0]} $$
 *
 * and its inverse:
 *
 * <table>
 *      <tr><td>&nbsp;</td><td>+ strand </td><td>- strand</td></tr>
 *      <tr><td>REF;</td><td>X[0][0]</td><td>X[0][1]</td></tr>
 *      <tr><td>ALT;</td><td>X[1][0]</td><td>X[1][1]</td></tr>
 * </table>
 *
 * The sum R + 1/R is used to detect a difference in strand bias for REF and for ALT (the sum makes it symmetric). A high value is indicative of large difference where one entry is very small compared to the others. A scale factor of refRatio/altRatio where
 * $$ refRatio = \frac{max(X[0][0], X[0][1])}{min(X[0][0], X[0][1} $$
 * and
 * $$ altRatio = \frac{max(X[1][0], X[1][1])}{min(X[1][0], X[1][1]} $$
 * ensures that the annotation value is large only.
 * </p>
 * <p>See the <a href="http://www.broadinstitute.org/gatk/guide/article?id=4732">method document on statistical tests</a> for a more detailed explanation of this statistical test.</p>
 *
 * <h3>Related annotations</h3>
 * <ul>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_StrandBiasBySample.php">StrandBiasBySample</a></b> outputs counts of read depth per allele for each strand orientation.</li>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_FisherStrand.php">FisherStrand</a></b> uses Fisher's Exact Test to evaluate strand bias.</li>
 * </ul>
 *
 */
public class StrandOddsRatio extends StrandBiasTest implements StandardAnnotation, ActiveRegionBasedAnnotation {
    private final static double AUGMENTATION_CONSTANT = 1.0;
    private static final int MIN_COUNT = 0;

    private static final String SOR = "SOR";

    @Override
    protected Map<String, Object> calculateAnnotationFromGTfield(GenotypesContext genotypes){
        final int[][] tableFromPerSampleAnnotations = getTableFromSamples( genotypes, MIN_COUNT );
        if ( tableFromPerSampleAnnotations != null ) {
            final double ratio = calculateSOR(tableFromPerSampleAnnotations);
            return annotationForOneTable(ratio);
        }
        return null;
    }

    @Override
    protected Map<String, Object> calculateAnnotationFromStratifiedContexts(Map<String, AlignmentContext> stratifiedContexts,
                                                                                     final VariantContext vc){
        final int[][] tableNoFiltering = getSNPContingencyTable(stratifiedContexts, vc.getReference(), vc.getAltAlleleWithHighestAlleleCount(), -1, MIN_COUNT);
        final double ratio = calculateSOR(tableNoFiltering);
        return annotationForOneTable(ratio);
    }

    @Override
    protected Map<String, Object> calculateAnnotationFromLikelihoodMap(Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap,
                                                                                final VariantContext vc){
        // either SNP with no alignment context, or indels: per-read likelihood map needed
        final int[][] table = getContingencyTable(stratifiedPerReadAlleleLikelihoodMap, vc, MIN_COUNT);
        final double ratio = calculateSOR(table);
        return annotationForOneTable(ratio);
    }

    /**
     * Computes the SOR value of a table after augmentation. Based on the symmetric odds ratio but modified to take on
     * low values when the reference +/- read count ratio is skewed but the alt count ratio is not.  Natural log is taken
     * to keep values within roughly the same range as other annotations.
     *
     * Augmentation avoids quotient by zero.
     *
     * @param originalTable The table before augmentation
     * @return the SOR annotation value
     */
    final protected double calculateSOR(final int[][] originalTable) {
        final double[][] augmentedTable = augmentContingencyTable(originalTable);

        double ratio = 0;

        ratio += (augmentedTable[0][0] / augmentedTable[0][1]) * (augmentedTable[1][1] / augmentedTable[1][0]);
        ratio += (augmentedTable[0][1] / augmentedTable[0][0]) * (augmentedTable[1][0] / augmentedTable[1][1]);

        final double refRatio = (Math.min(augmentedTable[0][0], augmentedTable[0][1])/Math.max(augmentedTable[0][0], augmentedTable[0][1]));
        final double altRatio = (Math.min(augmentedTable[1][0], augmentedTable[1][1])/Math.max(augmentedTable[1][0], augmentedTable[1][1]));

        ratio = ratio*refRatio/altRatio;

        return Math.log(ratio);
    }


    /**
     * Adds the small value AUGMENTATION_CONSTANT to all the entries of the table.
     *
     * @param table the table to augment
     * @return the augmented table
     */
    private static double[][] augmentContingencyTable(final int[][] table) {
        double[][] augmentedTable = new double[2][2];
        for ( int i = 0; i < 2; i++ ) {
            for ( int j = 0; j < 2; j++ )
                augmentedTable[i][j] = table[i][j] + AUGMENTATION_CONSTANT;
        }

        return augmentedTable;
    }

    /**
     * Returns an annotation result given a ratio
     *
     * @param ratio the symmetric odds ratio of the contingency table
     * @return a hash map from SOR
     */
    protected Map<String, Object> annotationForOneTable(final double ratio) {
        final Object value = String.format("%.3f", ratio);
        return Collections.singletonMap(SOR, value);
    }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Collections.singletonList(new VCFInfoHeaderLine(SOR, 1, VCFHeaderLineType.Float, "Symmetric Odds Ratio of 2x2 contingency table to detect strand bias"));
    }

    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(SOR);
    }
}
