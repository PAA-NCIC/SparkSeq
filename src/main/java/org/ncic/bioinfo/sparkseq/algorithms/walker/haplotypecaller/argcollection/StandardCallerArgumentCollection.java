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

import htsjdk.variant.variantcontext.VariantContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.basic.DefaultHashMap;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RodBinding;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.GenotypingOutputMode;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.OutputMode;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.afcalculate.AFCalculatorImplementation;

import java.io.File;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.lang.reflect.Modifier;
import java.util.Collections;
import java.util.Map;

/**
 * Author: wbc
 */
public class StandardCallerArgumentCollection implements Cloneable {

    public GenotypeCalculationArgumentCollection genotypeArgs = new GenotypeCalculationArgumentCollection();

    public GenotypingOutputMode genotypingOutputMode = GenotypingOutputMode.DISCOVERY;

    /**
     * When the UnifiedGenotyper is put into GENOTYPE_GIVEN_ALLELES mode it will genotype the samples using only the alleles provide in this rod binding
     */
    public RodBinding<VariantContext> alleles;

    /**
     * If this fraction is greater is than zero, the caller will aggressively attempt to remove contamination through biased down-sampling of reads.
     * Basically, it will ignore the contamination fraction of reads for each alternate allele.  So if the pileup contains N total bases, then we
     * will try to remove (N * contamination fraction) bases for each alternate allele.
     */
    public double CONTAMINATION_FRACTION = DEFAULT_CONTAMINATION_FRACTION;
    public static final double DEFAULT_CONTAMINATION_FRACTION = 0.0;

    /**
     * This argument specifies a file with two columns "sample" and "contamination" specifying the contamination level for those samples.
     * Samples that do not appear in this file will be processed with CONTAMINATION_FRACTION.
     **/
    public File CONTAMINATION_FRACTION_FILE = null;

    /**
     * Indicates whether there is some sample contamination present.
     */
    private boolean sampleContaminationWasLoaded = false;

    /**
     * @return an _Immutable_ copy of the Sample-Contamination Map, defaulting to CONTAMINATION_FRACTION so that if the sample isn't in the map map(sample)==CONTAMINATION_FRACTION
     */
    public Map<String, Double> getSampleContamination() {
        //make sure that the default value is set up right
        sampleContamination.setDefaultValue(CONTAMINATION_FRACTION);
        if (!Double.isNaN(CONTAMINATION_FRACTION) && CONTAMINATION_FRACTION > 0.0)
            sampleContaminationWasLoaded = true;
        return Collections.unmodifiableMap(sampleContamination);
    }

    public void setSampleContamination(DefaultHashMap<String, Double> sampleContamination) {
        this.sampleContamination.clear();
        this.sampleContaminationWasLoaded = !Double.isNaN(CONTAMINATION_FRACTION) && CONTAMINATION_FRACTION > 0.0;
        if (!sampleContaminationWasLoaded)
            for (final Double d : sampleContamination.values())
                if (!Double.isNaN(d) && d > 0.0) {
                    sampleContaminationWasLoaded = true;
                    break;
                }
        this.sampleContamination.putAll(sampleContamination);
        this.sampleContamination.setDefaultValue(CONTAMINATION_FRACTION);
    }

    /**
     * Returns true if there is some sample contamination present, false otherwise.
     *
     * @return {@code true} iff there is some sample contamination
     */
    public boolean isSampleContaminationPresent() {
        return (!Double.isNaN(CONTAMINATION_FRACTION) && CONTAMINATION_FRACTION > 0.0) || sampleContaminationWasLoaded;
    }

    //Needs to be here because it uses CONTAMINATION_FRACTION
    private DefaultHashMap<String, Double> sampleContamination = new DefaultHashMap<String, Double>(CONTAMINATION_FRACTION);

    /**
     * Controls the model used to calculate the probability that a site is variant plus the various sample genotypes in the data at a given locus.
     */
    public AFCalculatorImplementation requestedAlleleFrequencyCalculationModel;

    public File exactCallsLog = null;

    public OutputMode outputMode = OutputMode.EMIT_VARIANTS_ONLY;

    /**
     * Advanced, experimental argument: if SNP likelihood model is specified, and if EMIT_ALL_SITES output mode is set, when we set this argument then we will also emit PLs at all sites.
     * This will give a measure of reference confidence and a measure of which alt alleles are more plausible (if any).
     * WARNINGS:
     * - This feature will inflate VCF file size considerably.
     * - All SNP ALT alleles will be emitted with corresponding 10 PL values.
     * - An error will be emitted if EMIT_ALL_SITES is not set, or if anything other than diploid SNP model is used
     */
     public boolean annotateAllSitesWithPLs = false;

    /**
     * Creates a Standard caller argument collection with default values.
     */
    public StandardCallerArgumentCollection() {
    }

    /**
     * "Casts" a caller argument collection into another type.
     * <p>
     * <p>Common fields values are copied across</p>
     *
     * @param clazz the class of the result.
     * @param <T>   result argument collection class.
     * @return never {@code null}.
     */
    public <T extends StandardCallerArgumentCollection> T cloneTo(final Class<T> clazz) {
        // short cut: just use regular clone if it happens to be the same class.
        if (clazz == getClass())
            return (T) clone();
        try {
            final T result = clazz.newInstance();
            for (final Field field : getClass().getFields()) {
                // just copy common fields.
                if (!field.getDeclaringClass().isAssignableFrom(clazz))
                    continue;
                final int fieldModifiers = field.getModifiers();
                if ((fieldModifiers & UNCOPYABLE_MODIFIER_MASK) != 0) continue;
                //Use the clone() method if appropriate
                if (Cloneable.class.isAssignableFrom(field.getType())) {
                    Method clone = field.getType().getMethod("clone");
                    field.set(result, clone.invoke(field.get(this)));
                } else
                    field.set(result, field.get(this));
            }
            return result;
        } catch (final Exception ex) {
            throw new IllegalStateException(ex);
        }
    }

    /**
     * Creates a copy of this configuration.
     *
     * @return never {@code null}.
     */
    @Override
    public StandardCallerArgumentCollection clone() {
        try {
            StandardCallerArgumentCollection cloned = (StandardCallerArgumentCollection) super.clone();
            cloned.genotypeArgs = genotypeArgs.clone();
            return cloned;
        } catch (CloneNotSupportedException e) {
            throw new IllegalStateException("unreachable code");
        }
    }

    /**
     * Holds a modifiers mask that identifies those fields that cannot be copied between
     * StandardCallerArgumentCollections.
     */
    private final int UNCOPYABLE_MODIFIER_MASK = Modifier.PRIVATE | Modifier.STATIC | Modifier.FINAL;
}
