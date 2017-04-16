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
package org.ncic.bioinfo.sparkseq.algorithms.utils.vcfWriter;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GATKVariantContextUtils;
import org.ncic.bioinfo.sparkseq.exceptions.PipelineException;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

/**
 * Author: wbc
 */
public class GVCFWriter implements VariantContextWriter {
    //
    // static VCF field names
    //
    protected final static String MIN_DP_FORMAT_FIELD = "MIN_DP";

    final private List<HomRefBlock> GQPartitions;

    /** fields updated on the fly during GVCFWriter operation */
    int nextAvailableStart = -1;
    String contigOfNextAvailableStart = null;
    private String sampleName = null;
    private HomRefBlock currentBlock = null;
    private final int defaultPloidy;

    private List<VariantContext> resultGVCFWriter = new ArrayList<>();

    /**
     * Is the proposed GQ partitions well-formed?
     *
     * @param GQPartitions proposed GQ partitions
     * @return a non-null string if something is wrong (string explains issue)
     */
    protected static List<HomRefBlock> parsePartitions(final List<Integer> GQPartitions, final int defaultPloidy) {
        if ( GQPartitions == null ) throw new IllegalArgumentException("GQpartitions cannot be null");
        if ( GQPartitions.isEmpty() ) throw new IllegalArgumentException("GQpartitions cannot be empty");

        final List<HomRefBlock> result = new LinkedList<>();
        int lastThreshold = 0;
        for ( final Integer value : GQPartitions ) {
            if ( value == null ) throw new IllegalArgumentException("GQPartitions contains a null integer");
            if ( value < lastThreshold ) throw new IllegalArgumentException("GQPartitions is out of order.  Last is " + lastThreshold + " but next is " + value);
            if ( value == lastThreshold ) throw new IllegalArgumentException("GQPartitions is equal elements: Last is " + lastThreshold + " but next is " + value);
            result.add(new HomRefBlock(lastThreshold, value,defaultPloidy));
            lastThreshold = value;
        }
        result.add(new HomRefBlock(lastThreshold, Integer.MAX_VALUE,defaultPloidy));

        return result;
    }

    /**
     * Create a new GVCF writer
     *
     * Should be a non-empty list of boundaries.  For example, suppose this variable is
     *
     * [A, B, C]
     *
     * We would partition our hom-ref sites into the following bands:
     *
     * X < A
     * A <= X < B
     * B <= X < C
     * X >= C
     *
     * @param GQPartitions a well-formed list of GQ partitions
     * @param defaultPloidy the assumed ploidy for input variant context without one.
     */
    public GVCFWriter(final List<Integer> GQPartitions, final int defaultPloidy) {
        this.GQPartitions = parsePartitions(GQPartitions,defaultPloidy);
        this.defaultPloidy = defaultPloidy;
    }

    /**
     * Write the VCF header
     *
     * Adds standard GVCF fields to the header
     *
     * @param header a non-null header
     */
    @Override
    public void writeHeader(VCFHeader header) {
        throw new PipelineException("UnImplemented");
    }

    /**
     * Close this GVCF writer.  Finalizes any pending hom-ref blocks and emits those to the underlyingWriter as well
     */
    @Override
    public void close() {
        close(true);
    }

    /**
     * Horrible work around because there's no clean way to get our VCFWriter closed by the GATK
     *
     * If closeUnderlyingWriter is true, then we'll close the underlying writer, otherwise we'll leave it open
     * so the GATK closes it later
     *
     * @param closeUnderlyingWriter should we leave the underlying writer open or closed?
     */
    public void close(final boolean closeUnderlyingWriter) {
        emitCurrentBlock();
    }

    /**
     * Add hom-ref site from vc to this gVCF hom-ref state tracking, emitting any pending states if appropriate
     *
     * @param vc a non-null VariantContext
     * @param g a non-null genotype from VariantContext
     * @return a VariantContext to be emitted, or null if non is appropriate
     */
    protected VariantContext addHomRefSite(final VariantContext vc, final Genotype g) {

        if ( nextAvailableStart != -1 ) {
            // don't create blocks while the hom-ref site falls before nextAvailableStart (for deletions)
            if ( vc.getStart() <= nextAvailableStart && vc.getChr().equals(contigOfNextAvailableStart) )
                return null;
            // otherwise, reset to non-relevant
            nextAvailableStart = -1;
            contigOfNextAvailableStart = null;
        }

        final VariantContext result;
        if (genotypeCanBeMergedInCurrentBlock(g)) {
            currentBlock.add(vc.getStart(), g);
            result = null;
        } else {
            result = blockToVCF(currentBlock);
            currentBlock = createNewBlock(vc, g);
        }
        return result;
    }

    private boolean genotypeCanBeMergedInCurrentBlock(final Genotype g) {
        return currentBlock != null && currentBlock.withinBounds(g.getGQ()) && currentBlock.getPloidy() == g.getPloidy()
                && (currentBlock.getMinPLs() == null || !g.hasPL() || (currentBlock.getMinPLs().length == g.getPL().length));
    }

    /**
     * Flush the current hom-ref block, if necessary, to the underlying writer, and reset the currentBlock to null
     */
    private void emitCurrentBlock() {
        if ( currentBlock != null ) {
            // there's actually some work to do
            resultGVCFWriter.add(blockToVCF(currentBlock));
            currentBlock = null;
        }
    }

    /**
     * Convert a HomRefBlock into a VariantContext
     *
     * @param block the block to convert
     * @return a VariantContext representing the gVCF encoding for this block.
     * It will return {@code null} if input {@code block} is {@code null}, indicating that there
     * is no variant-context to be output into the VCF.
     */
    private VariantContext blockToVCF(final HomRefBlock block) {
        if ( block == null ) return null;

        final VariantContextBuilder vcb = new VariantContextBuilder(block.getStartingVC());
        vcb.attributes(new HashMap<String, Object>(2)); // clear the attributes
        vcb.stop(block.getStop());
        vcb.attribute(VCFConstants.END_KEY, block.getStop());

        // create the single Genotype with GQ and DP annotations
        final GenotypeBuilder gb = new GenotypeBuilder(sampleName, GATKVariantContextUtils.homozygousAlleleList(block.getRef(),block.getPloidy()));
        gb.noAD().noPL().noAttributes(); // clear all attributes
        gb.GQ(block.getMedianGQ());
        gb.DP(block.getMedianDP());
        gb.attribute(MIN_DP_FORMAT_FIELD, block.getMinDP());
        gb.PL(block.getMinPLs());

        // This annotation is no longer standard
        //gb.attribute(MIN_GQ_FORMAT_FIELD, block.getMinGQ());

        return vcb.genotypes(gb.make()).make();
    }

    /**
     * Helper function to create a new HomRefBlock from a variant context and current genotype
     *
     * @param vc the VariantContext at the site where want to start the band
     * @param g the genotype of the sample from vc that should be used to initialize the block
     * @return a newly allocated and initialized block containing g already
     */
    private HomRefBlock createNewBlock(final VariantContext vc, final Genotype g) {
        // figure out the GQ limits to use based on the GQ of g
        HomRefBlock partition = null;
        for ( final HomRefBlock maybePartition : GQPartitions ) {
            if ( maybePartition.withinBounds(g.getGQ()) ) {
                partition = maybePartition;
                break;
            }
        }

        if ( partition == null )
            throw new IllegalStateException("GQ " + g + " from " + vc + " didn't fit into any partition");

        // create the block, add g to it, and return it for use
        final HomRefBlock block = new HomRefBlock(vc, partition.getGQLowerBound(), partition.getGQUpperBound(), defaultPloidy);
        block.add(vc.getStart(), g);
        return block;
    }

    /**
     * Add a VariantContext to this writer for emission
     *
     * Requires that the VC have exactly one genotype
     *
     * @param vc a non-null VariantContext
     */
    @Override
    public void add(VariantContext vc) {
        if ( vc == null ) throw new IllegalArgumentException("vc cannot be null");

        if ( sampleName == null )
            sampleName = vc.getGenotype(0).getSampleName();

        if ( ! vc.hasGenotypes() ) {
            throw new IllegalArgumentException("GVCF assumes that the VariantContext has genotypes");
        } else if ( vc.getGenotypes().size() != 1 ) {
            throw new IllegalArgumentException("GVCF assumes that the VariantContext has exactly one genotype but saw " + vc.getGenotypes().size());
        } else {
            if ( currentBlock != null && ! currentBlock.isContiguous(vc) ) {
                // we've made a non-contiguous step (across interval, onto another chr), so finalize
                emitCurrentBlock();
            }

            final Genotype g = vc.getGenotype(0);
            if ( g.isHomRef() && vc.hasAlternateAllele(GATKVariantContextUtils.NON_REF_SYMBOLIC_ALLELE) && vc.isBiallelic() ) {
                // create bands
                final VariantContext maybeCompletedBand = addHomRefSite(vc, g);
                if ( maybeCompletedBand != null ) resultGVCFWriter.add(maybeCompletedBand);
            } else {
                // g is variant, so flush the bands and emit vc
                emitCurrentBlock();
                nextAvailableStart = vc.getEnd();
                contigOfNextAvailableStart = vc.getChr();
                resultGVCFWriter.add(vc);
            }
        }
    }

    public List<VariantContext> getResultGVCFWriter() {
        return resultGVCFWriter;
    }
}
