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
package org.ncic.bioinfo.sparkseq.algorithms.utils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.log4j.Logger;
import org.ncic.bioinfo.sparkseq.algorithms.utils.haplotype.Haplotype;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefMetaDataTracker;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RodBinding;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

/**
 * Author: wbc
 */
public class GenotypingGivenAllelesUtils {

    /**
     * Composes the given allele variant-context providing information about the rods and reference location.
     * @param tracker the meta data tracker.
     * @param loc the query location.
     * @param snpsOnly whether we only should consider SNP variation.
     * @param logger where to output warnings.
     * @param allelesBinding the target variation context binding containing the given alleles.
     * @return never {@code null}
     */
    public static VariantContext composeGivenAllelesVariantContextFromRod(final RefMetaDataTracker tracker,
                                                                          final GenomeLoc loc,
                                                                          final boolean snpsOnly,
                                                                          final Logger logger,
                                                                          final RodBinding<VariantContext> allelesBinding) {
        if ( tracker == null ) throw new IllegalArgumentException("the tracker cannot be null");
        if ( loc == null ) throw new IllegalArgumentException("the location cannot be null");
        if ( allelesBinding == null ) throw new IllegalArgumentException("the alleles binding cannot be null");

        VariantContext vc = null;

        // search for usable record
        for ( final VariantContext rodVc : tracker.getValues(allelesBinding, loc) ) {
            if ( rodVc != null && ! rodVc.isFiltered() && (! snpsOnly || rodVc.isSNP() )) {
                if ( vc == null )
                    vc = rodVc;
                else
                if (logger != null)
                    logger.warn("Multiple valid VCF records detected in the alleles input file at site "
                            + loc + ", only considering the first record");
            }
        }

        return vc;
    }

    /**
     * Create the list of artificial GGA-mode haplotypes by injecting each of the provided alternate alleles into the reference haplotype
     *
     * @param refHaplotype the reference haplotype
     * @param givenHaplotypes the list of alternate alleles in VariantContexts
     * @param activeRegionWindow the window containing the reference haplotype
     *
     * @return a non-null list of haplotypes
     */
    public static List<Haplotype> composeGivenHaplotypes(final Haplotype refHaplotype, final List<VariantContext> givenHaplotypes, final GenomeLoc activeRegionWindow) {
        if (refHaplotype == null) throw new IllegalArgumentException("the reference haplotype cannot be null");
        if (givenHaplotypes == null) throw new IllegalArgumentException("given haplotypes cannot be null");
        if (activeRegionWindow == null) throw new IllegalArgumentException("active region window cannot be null");
        if (activeRegionWindow.size() != refHaplotype.length()) throw new IllegalArgumentException("inconsistent reference haplotype and active region window");

        final Set<Haplotype> returnHaplotypes = new LinkedHashSet<>();
        final int activeRegionStart = refHaplotype.getAlignmentStartHapwrtRef();

        for( final VariantContext compVC : givenHaplotypes ) {
            if (!GATKVariantContextUtils.overlapsRegion(compVC, activeRegionWindow))
                throw new IllegalArgumentException(" some variant provided does not overlap with active region window");
            for( final Allele compAltAllele : compVC.getAlternateAlleles() ) {
                final Haplotype insertedRefHaplotype = refHaplotype.insertAllele(compVC.getReference(), compAltAllele, activeRegionStart + compVC.getStart() - activeRegionWindow.getStart(), compVC.getStart());
                if( insertedRefHaplotype != null ) { // can be null if the requested allele can't be inserted into the haplotype
                    returnHaplotypes.add(insertedRefHaplotype);
                }
            }
        }

        return new ArrayList<>(returnHaplotypes);
    }
}
