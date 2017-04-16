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

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.ncic.bioinfo.sparkseq.algorithms.utils.AlignmentUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.ReadUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.pairhmm.PairHMMIndelErrorModel;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.GATKSAMRecord;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.PileupElement;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.interfaces.StandardAnnotation;

import java.util.*;

/**
 * Rank Sum Test for relative positioning of REF vs. ALT alleles within reads
 *
 * <p>This variant-level annotation tests whether there is evidence of bias in the position of alleles within the reads that support them, between the reference and alternate alleles. Seeing an allele only near the ends of reads is indicative of error, because that is where sequencers tend to make the most errors. However, some variants located near the edges of sequenced regions will necessarily be covered by the ends of reads, so we can't just set an absolute "minimum distance from end of read" threshold. That is why we use a rank sum test to evaluate whether there is a difference in how well the reference allele and the alternate allele are supported.</p>
 *
 * <p>The ideal result is a value close to zero, which indicates there is little to no difference in where the alleles are found relative to the ends of reads. A negative value indicates that the alternate allele is found at the ends of reads more often than the reference allele. Conversely, a positive value indicates that the reference allele is found at the ends of reads more often than the alternate allele. </p>
 *
 * <p>This annotation can be used to evaluate confidence in a variant call and is a recommended covariate for variant recalibration (VQSR). Finding a statistically significant difference in relative position either way suggests that the sequencing process may have been biased or affected by an artifact. In practice, we only filter out low negative values when evaluating variant quality because the idea is to filter out variants for which the quality of the data supporting the alternate allele is comparatively low. The reverse case, where it is the quality of data supporting the reference allele that is lower (resulting in positive ranksum scores), is not really informative for filtering variants.</p>
 *
 * <h3>Statistical notes</h3>
 * <p>The value output for this annotation is the u-based z-approximation from the Mann-Whitney-Wilcoxon Rank Sum Test for site position within reads (position within reads supporting REF vs. position within reads supporting ALT). See the <a href="http://www.broadinstitute.org/gatk/guide/article?id=4732">method document on statistical tests</a> for a more detailed explanation of the ranksum test.</p>
 *
 * <h3>Caveat</h3>
 * <p>The read position rank sum test can not be calculated for sites without a mixture of reads showing both the reference and alternate alleles.</p>
 *
 */
public class ReadPosRankSumTest extends RankSumTest implements StandardAnnotation {

    @Override
    public List<String> getKeyNames() { return Arrays.asList("ReadPosRankSum"); }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFInfoHeaderLine("ReadPosRankSum", 1, VCFHeaderLineType.Float, "Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias"));
    }

    @Override
    protected Double getElementForRead(final GATKSAMRecord read, final int refLoc) {
        final int offset = ReadUtils.getReadCoordinateForReferenceCoordinate( read.getSoftStart(), read.getCigar(), refLoc, ReadUtils.ClippingTail.RIGHT_TAIL, true );
        if ( offset == ReadUtils.CLIPPING_GOAL_NOT_REACHED )
            return null;

        int readPos = AlignmentUtils.calcAlignmentByteArrayOffset( read.getCigar(), offset, false, 0, 0 );
        final int numAlignedBases = AlignmentUtils.getNumAlignedBasesCountingSoftClips( read );
        if (readPos > numAlignedBases / 2)
            readPos = numAlignedBases - (readPos + 1);
        return (double)readPos;
    }

    @Override
    protected Double getElementForPileupElement(final PileupElement p) {
        final int offset = AlignmentUtils.calcAlignmentByteArrayOffset(p.getRead().getCigar(), p, 0, 0);
        return (double)getFinalReadPosition(p.getRead(), offset);
    }

    @Override
    protected boolean isUsableBase(final PileupElement p) {
        return super.isUsableBase(p) && p.getRead().getCigar() != null;
    }

    @Override
    protected boolean isUsableRead(final GATKSAMRecord read, final int refLoc) {
        return super.isUsableRead(read, refLoc) && read.getSoftStart() + read.getCigar().getReadLength() > refLoc;
    }

    private int getFinalReadPosition(final GATKSAMRecord read, final int initialReadPosition) {
        final int numAlignedBases = getNumAlignedBases(read);

        int readPos = initialReadPosition;
        if (initialReadPosition > numAlignedBases / 2) {
            readPos = numAlignedBases - (initialReadPosition + 1);
        }
        return readPos;

    }

    private int getNumClippedBasesAtStart(final SAMRecord read) {
        // compute total number of clipped bases (soft or hard clipped)
        // check for hard clips (never consider these bases):
        final Cigar c = read.getCigar();
        final CigarElement first = c.getCigarElement(0);

        int numStartClippedBases = 0;
        if (first.getOperator() == CigarOperator.H) {
            numStartClippedBases = first.getLength();
        }
        final byte[] unclippedReadBases = read.getReadBases();
        final byte[] unclippedReadQuals = read.getBaseQualities();

        // Do a stricter base clipping than provided by CIGAR string, since this one may be too conservative,
        // and may leave a string of Q2 bases still hanging off the reads.
        for (int i = numStartClippedBases; i < unclippedReadBases.length; i++) {
            if (unclippedReadQuals[i] < PairHMMIndelErrorModel.BASE_QUAL_THRESHOLD)
                numStartClippedBases++;
            else
                break;

        }

        return numStartClippedBases;
    }

    private int getNumAlignedBases(final GATKSAMRecord read) {
        return read.getReadLength() - getNumClippedBasesAtStart(read) - getNumClippedBasesAtEnd(read);
    }

    private int getNumClippedBasesAtEnd(final GATKSAMRecord read) {
        // compute total number of clipped bases (soft or hard clipped)
        // check for hard clips (never consider these bases):
        final Cigar c = read.getCigar();
        CigarElement last = c.getCigarElement(c.numCigarElements() - 1);

        int numEndClippedBases = 0;
        if (last.getOperator() == CigarOperator.H) {
            numEndClippedBases = last.getLength();
        }
        final byte[] unclippedReadBases = read.getReadBases();
        final byte[] unclippedReadQuals = read.getBaseQualities();

        // Do a stricter base clipping than provided by CIGAR string, since this one may be too conservative,
        // and may leave a string of Q2 bases still hanging off the reads.
        for (int i = unclippedReadBases.length - numEndClippedBases - 1; i >= 0; i--) {
            if (unclippedReadQuals[i] < PairHMMIndelErrorModel.BASE_QUAL_THRESHOLD)
                numEndClippedBases++;
            else
                break;
        }

        return numEndClippedBases;
    }
}
