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
package org.ncic.bioinfo.sparkseq.algorithms.walker.indelrealigner;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import org.ncic.bioinfo.sparkseq.algorithms.engine.ReadWalker;
import org.ncic.bioinfo.sparkseq.algorithms.utils.AlignmentUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.BaseUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLoc;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLocParser;
import org.ncic.bioinfo.sparkseq.algorithms.utils.RandomGenerator;
import org.ncic.bioinfo.sparkseq.algorithms.utils.ReadUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.Utils;
import org.ncic.bioinfo.sparkseq.algorithms.data.basic.Pair;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefMetaDataTracker;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.ReferenceContext;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.GATKSAMRecord;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.SamContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.Filter;
import org.ncic.bioinfo.sparkseq.algorithms.utils.smithwaterman.Parameters;
import org.ncic.bioinfo.sparkseq.algorithms.utils.smithwaterman.SWPairwiseAlignment;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RODContentProvider;
import org.ncic.bioinfo.sparkseq.exceptions.UserException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

/**
 * Author: wbc
 */
public class IndelRealigner extends ReadWalker {

    protected int MAX_ISIZE_FOR_MOVEMENT = 3000;
    protected int MAX_POS_MOVE_ALLOWED = 200;
    protected int MAX_RECORDS_IN_MEMORY = 150000;
    protected int MAX_READS_FOR_CONSENSUSES = 120;
    protected int MAX_CONSENSUSES = 30;
    protected int MAX_READS = 20000;

    protected double MISMATCH_THRESHOLD = 0.15;
    protected double LOD_THRESHOLD = 5.0;

    protected boolean NO_ORIGINAL_ALIGNMENT_TAGS = false;
    public static final String ORIGINAL_CIGAR_TAG = "OC";
    public static final String ORIGINAL_POSITION_TAG = "OP";

    private static final int MAX_QUAL = 99;
    private static final double MISMATCH_COLUMN_CLEANED_FRACTION = 0.75;

    private final static Parameters swParameters = new Parameters(30, -10, -10, -2);

    private static final int REFERENCE_PADDING = 30;

    protected ConstrainedMateFixingManager manager = null;
    public ConsensusDeterminationModel consensusModel = ConsensusDeterminationModel.USE_READS;

    // the intervals input by the user
    private List<GenomeLoc> intervalList;
    private Iterator<GenomeLoc> intervals = null;

    // the current interval in the list
    private GenomeLoc currentInterval = null;
    private boolean sawReadInCurrentInterval = false;

    // the reads and known indels that fall into the current interval
    private ReadBin readsToClean;
    private final ArrayList<GATKSAMRecord> readsNotToClean = new ArrayList<GATKSAMRecord>();
    private final ArrayList<VariantContext> knownIndelsToTry = new ArrayList<VariantContext>();
    private final HashSet<Object> indelRodsSeen = new HashSet<Object>();
    private final HashSet<GATKSAMRecord> readsActuallyCleaned = new HashSet<GATKSAMRecord>();

    public IndelRealigner(GenomeLocParser genomeLocParser,
                          RefContentProvider refContentProvider,
                          SamContentProvider samContentProvider,
                          List<RODContentProvider> rodContentProviderList,
                          List<GenomeLoc> intervalList) {
        super(genomeLocParser, refContentProvider, samContentProvider, rodContentProviderList);
        this.intervalList = intervalList;
    }

    @Override
    public void initialize() {

        readsToClean = new ReadBin(genomeLocParser, REFERENCE_PADDING);

        intervals = intervalList.iterator();
        currentInterval = intervals.hasNext() ? intervals.next() : null;
        manager = new ConstrainedMateFixingManager(genomeLocParser, MAX_ISIZE_FOR_MOVEMENT, MAX_POS_MOVE_ALLOWED, MAX_RECORDS_IN_MEMORY);
    }

    @Override
    protected List<Filter> getFilter() {
        List<Filter> filters = new ArrayList<>();
        return filters;
    }

    @Override
    protected void map(final ReferenceContext ref,
                       final GATKSAMRecord read,
                       final RefMetaDataTracker metaDataTracker) {
        if (currentInterval == null) {
            emit(read);
            return;
        }

        // edge case: when the last target interval abuts the end of the genome, we'll get one of the
        //   unmapped reads while the currentInterval still isn't null.  We need to trigger the cleaning
        //   at this point without trying to create a GenomeLoc.
        if (read.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
            cleanAndCallMap(ref, read, metaDataTracker, null);
            return;
        }

        GenomeLoc readLoc = genomeLocParser.createGenomeLoc(read);
        // hack to get around unmapped reads having screwy locations
        if (readLoc.getStop() == 0)
            readLoc = genomeLocParser.createGenomeLoc(readLoc.getContig(), readLoc.getStart(), readLoc.getStart());

        if (readLoc.isBefore(currentInterval)) {
            if (!sawReadInCurrentInterval)
                emit(read);
            else
                readsNotToClean.add(read);
        } else if (readLoc.overlapsP(currentInterval)) {
            sawReadInCurrentInterval = true;

            if (doNotTryToClean(read)) {
                readsNotToClean.add(read);
            } else {
                readsToClean.add(read);

                // add the rods to the list of known variants
                populateKnownIndels(metaDataTracker);
            }

            if (readsToClean.size() + readsNotToClean.size() >= MAX_READS) {
                abortCleanForCurrentInterval();
            }
        } else {  // the read is past the current interval
            cleanAndCallMap(ref, read, metaDataTracker, readLoc);
        }

        return;
    }

    private void emit(final GATKSAMRecord read) {

        // check to see whether the read was modified by looking at the temporary tag
        boolean wasModified = readsActuallyCleaned.contains(read);

        try {
            manager.addRead(read, wasModified);
        } catch (RuntimeIOException e) {
            throw new UserException.ErrorWritingBamFile(e.getMessage());
        }
    }

    private void emitReadLists() {
        // pre-merge lists to sort them in preparation for constrained SAMFileWriter
        readsNotToClean.addAll(readsToClean.getReads());
        ReadUtils.sortReadsByCoordinate(readsNotToClean);
        manager.addReads(readsNotToClean, readsActuallyCleaned);
        readsToClean.clear();
        readsNotToClean.clear();
        readsActuallyCleaned.clear();
    }

    private boolean doNotTryToClean(GATKSAMRecord read) {
        return read.getReadUnmappedFlag() ||
                read.getNotPrimaryAlignmentFlag() ||
                read.getReadFailsVendorQualityCheckFlag() ||
                read.getMappingQuality() == 0 ||
                read.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START ||
                ConstrainedMateFixingManager.iSizeTooBigToMove(read, MAX_ISIZE_FOR_MOVEMENT) ||
                ReadUtils.is454Read(read) ||
                ReadUtils.isIonRead(read);
        // TODO -- it would be nice if we could use indels from 454/Ion reads as alternate consenses
    }

    private void populateKnownIndels(RefMetaDataTracker metaDataTracker) {
        for (final VariantContext vc : metaDataTracker.getValues(VariantContext.class)) {
            if (indelRodsSeen.contains(vc))
                continue;
            indelRodsSeen.add(vc);
            knownIndelsToTry.add(vc);
        }
    }

    private void abortCleanForCurrentInterval() {
        emitReadLists();
        currentInterval = intervals.hasNext() ? intervals.next() : null;
        sawReadInCurrentInterval = false;
    }

    private void cleanAndCallMap(ReferenceContext ref, GATKSAMRecord read, RefMetaDataTracker metaDataTracker, GenomeLoc readLoc) {
        if (readsToClean.size() > 0) {
            GenomeLoc earliestPossibleMove = genomeLocParser.createGenomeLoc(readsToClean.getReads().get(0));
            if (manager.canMoveReads(earliestPossibleMove))
                clean(readsToClean);
        }
        knownIndelsToTry.clear();
        indelRodsSeen.clear();

        emitReadLists();
        do {
            currentInterval = intervals.hasNext() ? intervals.next() : null;

        } while (currentInterval != null && (readLoc == null || currentInterval.isBefore(readLoc)));

        sawReadInCurrentInterval = false;

        // call back into map now that the state has been updated
        map(ref, read, metaDataTracker);
    }

    private static int mismatchQualitySumIgnoreCigar(final AlignedRead aRead, final byte[] refSeq, int refIndex, int quitAboveThisValue) {
        final byte[] readSeq = aRead.getReadBases();
        final byte[] quals = aRead.getBaseQualities();
        int sum = 0;
        for (int readIndex = 0; readIndex < readSeq.length; refIndex++, readIndex++) {
            if (refIndex >= refSeq.length) {
                sum += MAX_QUAL;
                // optimization: once we pass the threshold, stop calculating
                if (sum > quitAboveThisValue)
                    return sum;
            } else {
                byte refChr = refSeq[refIndex];
                byte readChr = readSeq[readIndex];
                if (!BaseUtils.isRegularBase(readChr) || !BaseUtils.isRegularBase(refChr))
                    continue; // do not count Ns/Xs/etc ?
                if (readChr != refChr) {
                    sum += (int) quals[readIndex];
                    // optimization: once we pass the threshold, stop calculating
                    if (sum > quitAboveThisValue)
                        return sum;
                }
            }
        }
        return sum;
    }

    private void clean(ReadBin readsToClean) {

        final List<GATKSAMRecord> reads = readsToClean.getReads();
        if (reads.size() == 0)
            return;

        byte[] reference = readsToClean.getReference(refContentProvider);
        int leftmostIndex = readsToClean.getLocation().getStart();

        final ArrayList<GATKSAMRecord> refReads = new ArrayList<GATKSAMRecord>();                 // reads that perfectly match ref
        final ArrayList<AlignedRead> altReads = new ArrayList<AlignedRead>();               // reads that don't perfectly match
        final LinkedList<AlignedRead> altAlignmentsToTest = new LinkedList<AlignedRead>();  // should we try to make an alt consensus from the read?
        final Set<Consensus> altConsenses = new LinkedHashSet<Consensus>();               // list of alt consenses

        // if there are any known indels for this region, get them and create alternate consenses
        generateAlternateConsensesFromKnownIndels(altConsenses, leftmostIndex, reference);

        // decide which reads potentially need to be cleaned;
        // if there are reads with a single indel in them, add that indel to the list of alternate consenses
        long totalRawMismatchSum = determineReadsThatNeedCleaning(reads, refReads, altReads, altAlignmentsToTest, altConsenses, leftmostIndex, reference);

        // use 'Smith-Waterman' to create alternate consenses from reads that mismatch the reference, using totalRawMismatchSum as the random seed
        if (consensusModel == ConsensusDeterminationModel.USE_SW)
            generateAlternateConsensesFromReads(altAlignmentsToTest, altConsenses, reference, leftmostIndex);

        Consensus bestConsensus = null;

        for (Consensus consensus : altConsenses) {

            for (int j = 0; j < altReads.size(); j++) {
                AlignedRead toTest = altReads.get(j);
                Pair<Integer, Integer> altAlignment = findBestOffset(consensus.str, toTest, leftmostIndex);

                // the mismatch score is the min of its alignment vs. the reference and vs. the alternate
                int myScore = altAlignment.second;

                if (myScore > toTest.getAlignerMismatchScore() || myScore >= toTest.getMismatchScoreToReference())
                    myScore = toTest.getMismatchScoreToReference();
                    // keep track of reads that align better to the alternate consensus.
                    // By pushing alignments with equal scores to the alternate, it means we'll over-call (het -> hom non ref) but are less likely to under-call (het -> ref, het non ref -> het)
                else
                    consensus.readIndexes.add(new Pair<Integer, Integer>(j, altAlignment.first));

                //logger.debug(consensus.cigar +  " vs. " + toTest.getRead().getReadName() + "-" + toTest.getRead().getReadString() + " => " + myScore + " vs. " + toTest.getMismatchScoreToReference());
                if (!toTest.getRead().getDuplicateReadFlag())
                    consensus.mismatchSum += myScore;

                // optimization: once the mismatch sum is higher than the best consensus, quit since this one can't win
                //  THIS MUST BE DISABLED IF WE DECIDE TO ALLOW MORE THAN ONE ALTERNATE CONSENSUS!
                if (bestConsensus != null && consensus.mismatchSum > bestConsensus.mismatchSum)
                    break;
            }

            if (bestConsensus == null || bestConsensus.mismatchSum > consensus.mismatchSum) {
                // we do not need this alt consensus, release memory right away!!
                if (bestConsensus != null)
                    bestConsensus.readIndexes.clear();
                bestConsensus = consensus;
            } else {
                // we do not need this alt consensus, release memory right away!!
                consensus.readIndexes.clear();
            }
        }

        // if:
        // 1) the best alternate consensus has a smaller sum of quality score mismatches than the aligned version of the reads,
        // 2) beats the LOD threshold for the sum of quality score mismatches of the raw version of the reads,
        // 3) didn't just move around the mismatching columns (i.e. it actually reduces entropy),
        // then clean!
        final double improvement = (bestConsensus == null ? -1 : ((double) (totalRawMismatchSum - bestConsensus.mismatchSum)) / 10.0);
        if (improvement >= LOD_THRESHOLD) {

            bestConsensus.cigar = AlignmentUtils.leftAlignIndel(bestConsensus.cigar, reference, bestConsensus.str, bestConsensus.positionOnReference, bestConsensus.positionOnReference, true);

            // start cleaning the appropriate reads
            for (Pair<Integer, Integer> indexPair : bestConsensus.readIndexes) {
                AlignedRead aRead = altReads.get(indexPair.first);
                if (!updateRead(bestConsensus.cigar, bestConsensus.positionOnReference, indexPair.second, aRead, leftmostIndex))
                    return;
            }
            if (consensusModel != ConsensusDeterminationModel.KNOWNS_ONLY && !alternateReducesEntropy(altReads, reference, leftmostIndex)) {

            } else {
                // finish cleaning the appropriate reads
                for (Pair<Integer, Integer> indexPair : bestConsensus.readIndexes) {
                    final AlignedRead aRead = altReads.get(indexPair.first);
                    if (aRead.finalizeUpdate()) {
                        // We need to update the mapping quality score of the cleaned reads;
                        // however we don't have enough info to use the proper MAQ scoring system.
                        // For now, we will just arbitrarily add 10 to the mapping quality. [EB, 6/7/2010].
                        // TODO -- we need a better solution here
                        GATKSAMRecord read = aRead.getRead();
                        if (read.getMappingQuality() != 255) // 255 == Unknown, so don't modify it
                            read.setMappingQuality(Math.min(aRead.getRead().getMappingQuality() + 10, 254));

                        // before we fix the attribute tags we first need to make sure we have enough of the reference sequence
                        int neededBasesToLeft = leftmostIndex - read.getAlignmentStart();
                        int neededBasesToRight = read.getAlignmentEnd() - leftmostIndex - reference.length + 1;
                        int neededBases = Math.max(neededBasesToLeft, neededBasesToRight);
                        if (neededBases > 0) {
                            GenomeLoc newLoc = new GenomeLoc(currentInterval.getContig(), currentInterval.getContigIndex(), leftmostIndex, leftmostIndex + reference.length);
                            ReferenceContext referenceContext = refContentProvider.getReferenceContext(newLoc, neededBases);
                            reference = referenceContext.getBases();
                            leftmostIndex = referenceContext.getLocus().getStart();
                        }

                        // now, fix the attribute tags
                        // TODO -- get rid of this try block when Picard does the right thing for reads aligned off the end of the reference
                        try {
                            if (read.getAttribute(SAMTag.NM.name()) != null)
                                read.setAttribute(SAMTag.NM.name(), SequenceUtil.calculateSamNmTag(read, reference, leftmostIndex - 1));
                            if (read.getAttribute(SAMTag.UQ.name()) != null)
                                read.setAttribute(SAMTag.UQ.name(), SequenceUtil.sumQualitiesOfMismatches(read, reference, leftmostIndex - 1));
                        } catch (Exception e) {
                            // ignore it
                        }
                        // TODO -- this is only temporary until Tim adds code to recalculate this value
                        if (read.getAttribute(SAMTag.MD.name()) != null)
                            read.setAttribute(SAMTag.MD.name(), null);

                        // mark that it was actually cleaned
                        readsActuallyCleaned.add(read);
                    }
                }
            }
        }
    }

    private void generateAlternateConsensesFromKnownIndels(final Set<Consensus> altConsensesToPopulate, final int leftmostIndex, final byte[] reference) {
        for (VariantContext knownIndel : knownIndelsToTry) {
            if (knownIndel == null || !knownIndel.isIndel() || knownIndel.isComplexIndel())
                continue;
            final byte[] indelStr;
            if (knownIndel.isSimpleInsertion()) {
                final byte[] fullAllele = knownIndel.getAlternateAllele(0).getBases();
                indelStr = Arrays.copyOfRange(fullAllele, 1, fullAllele.length); // remove ref padding
            } else {
                indelStr = Utils.dupBytes((byte) '-', knownIndel.getReference().length() - 1);
            }
            int start = knownIndel.getStart() - leftmostIndex + 1;
            Consensus c = createAlternateConsensus(start, reference, indelStr, knownIndel);
            if (c != null)
                altConsensesToPopulate.add(c);
        }
    }

    private long determineReadsThatNeedCleaning(final List<GATKSAMRecord> reads,
                                                final ArrayList<GATKSAMRecord> refReadsToPopulate,
                                                final ArrayList<AlignedRead> altReadsToPopulate,
                                                final LinkedList<AlignedRead> altAlignmentsToTest,
                                                final Set<Consensus> altConsenses,
                                                final int leftmostIndex,
                                                final byte[] reference) {

        long totalRawMismatchSum = 0L;
        for (final GATKSAMRecord read : reads) {

            // we can not deal with screwy records
            if (read.getCigar().numCigarElements() == 0) {
                refReadsToPopulate.add(read);
                continue;
            }

            final AlignedRead aRead = new AlignedRead(read);

            // first, move existing indels (for 1 indel reads only) to leftmost position within identical sequence
            int numBlocks = AlignmentUtils.getNumAlignmentBlocks(read);
            if (numBlocks == 2) {
                Cigar newCigar = AlignmentUtils.leftAlignIndel(unclipCigar(read.getCigar()), reference, read.getReadBases(), read.getAlignmentStart() - leftmostIndex, 0, true);
                aRead.setCigar(newCigar, false);
            }

            final int startOnRef = read.getAlignmentStart() - leftmostIndex;
            final int rawMismatchScore = mismatchQualitySumIgnoreCigar(aRead, reference, startOnRef, Integer.MAX_VALUE);

            // if this doesn't match perfectly to the reference, let's try to clean it
            if (rawMismatchScore > 0) {
                altReadsToPopulate.add(aRead);
                //logger.debug("Adding " + read.getReadName() + " with raw mismatch score " + rawMismatchScore + " to non-ref reads");

                if (!read.getDuplicateReadFlag())
                    totalRawMismatchSum += rawMismatchScore;
                aRead.setMismatchScoreToReference(rawMismatchScore);
                aRead.setAlignerMismatchScore(AlignmentUtils.mismatchingQualities(aRead.getRead(), reference, startOnRef));

                // if it has an indel, let's see if that's the best consensus
                if (consensusModel != ConsensusDeterminationModel.KNOWNS_ONLY && numBlocks == 2) {
                    Consensus c = createAlternateConsensus(startOnRef, aRead.getCigar(), reference, aRead.getReadBases());
                    if (c != null)
                        altConsenses.add(c);
                } else {
                    altAlignmentsToTest.add(aRead);
                }
            }
            // otherwise, we can emit it as is
            else {
                //logger.debug("Adding " + read.getReadName() + " with raw mismatch score " + rawMismatchScore + " to ref reads");
                refReadsToPopulate.add(read);
            }
        }

        return totalRawMismatchSum;
    }

    private void generateAlternateConsensesFromReads(final LinkedList<AlignedRead> altAlignmentsToTest,
                                                     final Set<Consensus> altConsensesToPopulate,
                                                     final byte[] reference,
                                                     final int leftmostIndex) {

        // if we are under the limit, use all reads to generate alternate consenses
        if (altAlignmentsToTest.size() <= MAX_READS_FOR_CONSENSUSES) {
            for (AlignedRead aRead : altAlignmentsToTest) {
                createAndAddAlternateConsensus(aRead.getReadBases(), altConsensesToPopulate, reference);
            }
        }
        // otherwise, choose reads for alternate consenses randomly
        else {
            int readsSeen = 0;
            while (readsSeen++ < MAX_READS_FOR_CONSENSUSES && altConsensesToPopulate.size() <= MAX_CONSENSUSES) {
                int index = RandomGenerator.getRandomGenerator().nextInt(altAlignmentsToTest.size());
                AlignedRead aRead = altAlignmentsToTest.remove(index);

                createAndAddAlternateConsensus(aRead.getReadBases(), altConsensesToPopulate, reference);
            }
        }
    }

    private void createAndAddAlternateConsensus(final byte[] read, final Set<Consensus> altConsensesToPopulate, final byte[] reference) {

        // do a pairwise alignment against the reference
        SWPairwiseAlignment swConsensus = new SWPairwiseAlignment(reference, read, swParameters);
        Consensus c = createAlternateConsensus(swConsensus.getAlignmentStart2wrt1(), swConsensus.getCigar(), reference, read);
        if (c != null)
            altConsensesToPopulate.add(c);
    }

    // create a Consensus from cigar/read strings which originate somewhere on the reference
    private Consensus createAlternateConsensus(final int indexOnRef, final Cigar c, final byte[] reference, final byte[] readStr) {
        if (indexOnRef < 0)
            return null;

        // if there are no indels, we do not need this consensus, can abort early:
        if (c.numCigarElements() == 1 && c.getCigarElement(0).getOperator() == CigarOperator.M)
            return null;

        // create the new consensus
        ArrayList<CigarElement> elements = new ArrayList<CigarElement>(c.numCigarElements() - 1);
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < indexOnRef; i++)
            sb.append((char) reference[i]);

        int indelCount = 0;
        int altIdx = 0;
        int refIdx = indexOnRef;
        boolean ok_flag = true;
        for (int i = 0; i < c.numCigarElements(); i++) {
            CigarElement ce = c.getCigarElement(i);
            int elementLength = ce.getLength();
            switch (ce.getOperator()) {
                case D:
                    refIdx += elementLength;
                    indelCount++;
                    elements.add(ce);
                    break;
                case M:
                case EQ:
                case X:
                    altIdx += elementLength;
                case N:
                    if (reference.length < refIdx + elementLength)
                        ok_flag = false;
                    else {
                        for (int j = 0; j < elementLength; j++)
                            sb.append((char) reference[refIdx + j]);
                    }
                    refIdx += elementLength;
                    elements.add(new CigarElement(elementLength, CigarOperator.M));
                    break;
                case I:
                    for (int j = 0; j < elementLength; j++) {
                        if (!BaseUtils.isRegularBase(readStr[altIdx + j])) {
                            // Insertions with N's in them cause real problems sometimes; it's better to drop them altogether
                            ok_flag = false;
                            break;
                        }
                        sb.append((char) readStr[altIdx + j]);
                    }
                    altIdx += elementLength;
                    indelCount++;
                    elements.add(ce);
                    break;
                case S:
                default:
                    break;
            }
        }
        // make sure that there is at most only a single indel and it aligns appropriately!
        if (!ok_flag || indelCount != 1 || reference.length < refIdx)
            return null;

        for (int i = refIdx; i < reference.length; i++)
            sb.append((char) reference[i]);
        byte[] altConsensus = StringUtil.stringToBytes(sb.toString()); // alternative consensus sequence we just built from the current read

        return new Consensus(altConsensus, new Cigar(elements), indexOnRef);
    }

    private Consensus createAlternateConsensus(final int indexOnRef, final byte[] reference, final byte[] indelStr, final VariantContext indel) {
        if (indexOnRef < 0 || indexOnRef >= reference.length)
            return null;

        // create the new consensus
        StringBuilder sb = new StringBuilder();
        Cigar cigar = new Cigar();
        int refIdx;

        for (refIdx = 0; refIdx < indexOnRef; refIdx++)
            sb.append((char) reference[refIdx]);
        if (indexOnRef > 0)
            cigar.add(new CigarElement(indexOnRef, CigarOperator.M));

        if (indel.isSimpleDeletion()) {
            refIdx += indelStr.length;
            cigar.add(new CigarElement(indelStr.length, CigarOperator.D));
        } else if (indel.isSimpleInsertion()) {
            for (byte b : indelStr)
                sb.append((char) b);
            cigar.add(new CigarElement(indelStr.length, CigarOperator.I));
        } else {
            throw new IllegalStateException("Creating an alternate consensus from a complex indel is not allows");
        }

        if (reference.length - refIdx > 0)
            cigar.add(new CigarElement(reference.length - refIdx, CigarOperator.M));
        for (; refIdx < reference.length; refIdx++)
            sb.append((char) reference[refIdx]);
        byte[] altConsensus = StringUtil.stringToBytes(sb.toString()); // alternative consensus sequence we just built from the current read

        return new Consensus(altConsensus, cigar, 0);
    }

    private Pair<Integer, Integer> findBestOffset(final byte[] ref, final AlignedRead read, final int leftmostIndex) {

        // optimization: try the most likely alignment first (to get a low score to beat)
        int originalAlignment = read.getOriginalAlignmentStart() - leftmostIndex;
        int bestScore = mismatchQualitySumIgnoreCigar(read, ref, originalAlignment, Integer.MAX_VALUE);
        int bestIndex = originalAlignment;

        // optimization: we can't get better than 0, so we can quit now
        if (bestScore == 0)
            return new Pair<Integer, Integer>(bestIndex, 0);

        // optimization: the correct alignment shouldn't be too far from the original one (or else the read wouldn't have aligned in the first place)
        for (int i = 0; i < originalAlignment; i++) {
            int score = mismatchQualitySumIgnoreCigar(read, ref, i, bestScore);
            if (score < bestScore) {
                bestScore = score;
                bestIndex = i;
            }
            // optimization: we can't get better than 0, so we can quit now
            if (bestScore == 0)
                return new Pair<Integer, Integer>(bestIndex, 0);
        }

        final int maxPossibleStart = ref.length - read.getReadLength();
        for (int i = originalAlignment + 1; i <= maxPossibleStart; i++) {
            int score = mismatchQualitySumIgnoreCigar(read, ref, i, bestScore);
            if (score < bestScore) {
                bestScore = score;
                bestIndex = i;
            }
            // optimization: we can't get better than 0, so we can quit now
            if (bestScore == 0)
                return new Pair<Integer, Integer>(bestIndex, 0);
        }

        return new Pair<Integer, Integer>(bestIndex, bestScore);
    }


    private boolean updateRead(final Cigar altCigar, final int altPosOnRef, final int myPosOnAlt, final AlignedRead aRead, final int leftmostIndex) {
        Cigar readCigar = new Cigar();

        // special case: there is no indel
        if (altCigar.getCigarElements().size() == 1) {
            aRead.setAlignmentStart(leftmostIndex + myPosOnAlt);
            readCigar.add(new CigarElement(aRead.getReadLength(), CigarOperator.M));
            aRead.setCigar(readCigar);
            return true;
        }

        CigarElement altCE1 = altCigar.getCigarElement(0);
        CigarElement altCE2 = altCigar.getCigarElement(1);

        int leadingMatchingBlockLength = 0; // length of the leading M element or 0 if the leading element is I

        CigarElement indelCE;
        if (altCE1.getOperator() == CigarOperator.I) {
            indelCE = altCE1;
            if (altCE2.getOperator() != CigarOperator.M) {
                return false;
            }
        } else {
            if (altCE1.getOperator() != CigarOperator.M) {
                return false;
            }
            if (altCE2.getOperator() == CigarOperator.I || altCE2.getOperator() == CigarOperator.D) {
                indelCE = altCE2;
            } else {
                return false;
            }
            leadingMatchingBlockLength = altCE1.getLength();
        }

        // the easiest thing to do is to take each case separately
        int endOfFirstBlock = altPosOnRef + leadingMatchingBlockLength;
        boolean sawAlignmentStart = false;

        // for reads starting before the indel
        if (myPosOnAlt < endOfFirstBlock) {
            aRead.setAlignmentStart(leftmostIndex + myPosOnAlt);
            sawAlignmentStart = true;

            // for reads ending before the indel
            if (myPosOnAlt + aRead.getReadLength() <= endOfFirstBlock) {
                //readCigar.add(new CigarElement(aRead.getReadLength(), CigarOperator.M));
                //aRead.setCigar(readCigar);
                aRead.setCigar(null); // reset to original alignment
                return true;
            }
            readCigar.add(new CigarElement(endOfFirstBlock - myPosOnAlt, CigarOperator.M));
        }

        // forward along the indel
        //int indelOffsetOnRef = 0, indelOffsetOnRead = 0;
        if (indelCE.getOperator() == CigarOperator.I) {
            // for reads that end in an insertion
            if (myPosOnAlt + aRead.getReadLength() < endOfFirstBlock + indelCE.getLength()) {
                int partialInsertionLength = myPosOnAlt + aRead.getReadLength() - endOfFirstBlock;
                // if we also started inside the insertion, then we need to modify the length
                if (!sawAlignmentStart)
                    partialInsertionLength = aRead.getReadLength();
                readCigar.add(new CigarElement(partialInsertionLength, CigarOperator.I));
                aRead.setCigar(readCigar);
                return true;
            }

            // for reads that start in an insertion
            if (!sawAlignmentStart && myPosOnAlt < endOfFirstBlock + indelCE.getLength()) {
                aRead.setAlignmentStart(leftmostIndex + endOfFirstBlock);
                readCigar.add(new CigarElement(indelCE.getLength() - (myPosOnAlt - endOfFirstBlock), CigarOperator.I));
                //indelOffsetOnRead = myPosOnAlt - endOfFirstBlock;
                sawAlignmentStart = true;
            } else if (sawAlignmentStart) {
                readCigar.add(indelCE);
                //indelOffsetOnRead = indelCE.getLength();
            }
        } else if (indelCE.getOperator() == CigarOperator.D) {
            if (sawAlignmentStart)
                readCigar.add(indelCE);
            //indelOffsetOnRef = indelCE.getLength();
        }

        // for reads that start after the indel
        if (!sawAlignmentStart) {
            //aRead.setAlignmentStart(leftmostIndex + myPosOnAlt + indelOffsetOnRef - indelOffsetOnRead);
            //readCigar.add(new CigarElement(aRead.getReadLength(), CigarOperator.M));
            //aRead.setCigar(readCigar);
            aRead.setCigar(null); // reset to original alignment
            return true;
        }

        int readRemaining = aRead.getReadBases().length;
        for (CigarElement ce : readCigar.getCigarElements()) {
            if (ce.getOperator() != CigarOperator.D)
                readRemaining -= ce.getLength();
        }
        if (readRemaining > 0)
            readCigar.add(new CigarElement(readRemaining, CigarOperator.M));
        aRead.setCigar(readCigar);

        return true;
    }

    private boolean alternateReducesEntropy(final List<AlignedRead> reads, final byte[] reference, final int leftmostIndex) {
        final int[] originalMismatchBases = new int[reference.length];
        final int[] cleanedMismatchBases = new int[reference.length];
        final int[] totalOriginalBases = new int[reference.length];
        final int[] totalCleanedBases = new int[reference.length];

        // set to 1 to prevent dividing by zero
        for (int i = 0; i < reference.length; i++)
            originalMismatchBases[i] = totalOriginalBases[i] = cleanedMismatchBases[i] = totalCleanedBases[i] = 0;

        for (final AlignedRead read : reads) {
            if (read.getRead().getAlignmentBlocks().size() > 1)
                continue;

            int refIdx = read.getOriginalAlignmentStart() - leftmostIndex;
            final byte[] readStr = read.getReadBases();
            final byte[] quals = read.getBaseQualities();

            for (int j = 0; j < readStr.length; j++, refIdx++) {
                if (refIdx < 0 || refIdx >= reference.length) {
                    //System.out.println( "Read: "+read.getRead().getReadName() + "; length = " + readStr.length() );
                    //System.out.println( "Ref left: "+ leftmostIndex +"; ref length=" + reference.length() + "; read alignment start: "+read.getOriginalAlignmentStart() );
                    break;
                }
                totalOriginalBases[refIdx] += quals[j];
                if (readStr[j] != reference[refIdx])
                    originalMismatchBases[refIdx] += quals[j];
            }

            // reset and now do the calculation based on the cleaning
            refIdx = read.getAlignmentStart() - leftmostIndex;
            int altIdx = 0;
            Cigar c = read.getCigar();
            for (int j = 0; j < c.numCigarElements(); j++) {
                CigarElement ce = c.getCigarElement(j);
                int elementLength = ce.getLength();
                switch (ce.getOperator()) {
                    case M:
                    case EQ:
                    case X:
                        for (int k = 0; k < elementLength; k++, refIdx++, altIdx++) {
                            if (refIdx >= reference.length)
                                break;
                            totalCleanedBases[refIdx] += quals[altIdx];
                            if (readStr[altIdx] != reference[refIdx])
                                cleanedMismatchBases[refIdx] += quals[altIdx];
                        }
                        break;
                    case I:
                        altIdx += elementLength;
                        break;
                    case D:
                        refIdx += elementLength;
                        break;
                    case S:
                    default:
                        break;
                }
            }
        }

        int originalMismatchColumns = 0, cleanedMismatchColumns = 0;
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < reference.length; i++) {
            if (cleanedMismatchBases[i] == originalMismatchBases[i])
                continue;
            boolean didMismatch = false, stillMismatches = false;
            if (originalMismatchBases[i] > totalOriginalBases[i] * MISMATCH_THRESHOLD) {
                didMismatch = true;
                originalMismatchColumns++;
                if (totalCleanedBases[i] > 0 && ((double) cleanedMismatchBases[i] / (double) totalCleanedBases[i]) > ((double) originalMismatchBases[i] / (double) totalOriginalBases[i]) * (1.0 - MISMATCH_COLUMN_CLEANED_FRACTION)) {
                    stillMismatches = true;
                    cleanedMismatchColumns++;
                }
            } else if (cleanedMismatchBases[i] > totalCleanedBases[i] * MISMATCH_THRESHOLD) {
                cleanedMismatchColumns++;
            }
        }

        final boolean reduces = (originalMismatchColumns == 0 || cleanedMismatchColumns < originalMismatchColumns);

        return reduces;
    }

    protected static Cigar unclipCigar(Cigar cigar) {
        ArrayList<CigarElement> elements = new ArrayList<CigarElement>(cigar.numCigarElements());
        for (CigarElement ce : cigar.getCigarElements()) {
            if (!isClipOperator(ce.getOperator()))
                elements.add(ce);
        }
        return new Cigar(elements);
    }

    private static boolean isClipOperator(CigarOperator op) {
        return op == CigarOperator.S || op == CigarOperator.H || op == CigarOperator.P;
    }

    protected static Cigar reclipCigar(Cigar cigar, SAMRecord read) {
        ArrayList<CigarElement> elements = new ArrayList<CigarElement>();

        int i = 0;
        int n = read.getCigar().numCigarElements();
        while (i < n && isClipOperator(read.getCigar().getCigarElement(i).getOperator()))
            elements.add(read.getCigar().getCigarElement(i++));

        elements.addAll(cigar.getCigarElements());

        i++;
        while (i < n && !isClipOperator(read.getCigar().getCigarElement(i).getOperator()))
            i++;

        while (i < n && isClipOperator(read.getCigar().getCigarElement(i).getOperator()))
            elements.add(read.getCigar().getCigarElement(i++));

        return new Cigar(elements);
    }

    private class AlignedRead {
        private final GATKSAMRecord read;
        private byte[] readBases = null;
        private byte[] baseQuals = null;
        private Cigar newCigar = null;
        private int newStart = -1;
        private int mismatchScoreToReference = 0;
        private long alignerMismatchScore = 0;

        public AlignedRead(GATKSAMRecord read) {
            this.read = read;
            mismatchScoreToReference = 0;
        }

        public GATKSAMRecord getRead() {
            return read;
        }

        public int getReadLength() {
            return readBases != null ? readBases.length : read.getReadLength();
        }

        public byte[] getReadBases() {
            if (readBases == null)
                getUnclippedBases();
            return readBases;
        }

        public byte[] getBaseQualities() {
            if (baseQuals == null)
                getUnclippedBases();
            return baseQuals;
        }

        // pull out the bases that aren't clipped out
        private void getUnclippedBases() {
            readBases = new byte[getReadLength()];
            baseQuals = new byte[getReadLength()];
            byte[] actualReadBases = read.getReadBases();
            byte[] actualBaseQuals = read.getBaseQualities();
            int fromIndex = 0, toIndex = 0;

            for (CigarElement ce : read.getCigar().getCigarElements()) {
                int elementLength = ce.getLength();
                switch (ce.getOperator()) {
                    case S:
                        fromIndex += elementLength;
                        break;
                    case M:
                    case EQ:
                    case X:
                    case I:
                        if (fromIndex + elementLength > actualReadBases.length)
                            throw new UserException.MalformedBAM(read, "the CIGAR string is inconsistent with the number of bases in the read");
                        System.arraycopy(actualReadBases, fromIndex, readBases, toIndex, elementLength);
                        System.arraycopy(actualBaseQuals, fromIndex, baseQuals, toIndex, elementLength);
                        fromIndex += elementLength;
                        toIndex += elementLength;
                    default:
                        break;
                }
            }

            // if we got clipped, trim the array
            if (fromIndex != toIndex) {
                byte[] trimmedRB = new byte[toIndex];
                byte[] trimmedBQ = new byte[toIndex];
                System.arraycopy(readBases, 0, trimmedRB, 0, toIndex);
                System.arraycopy(baseQuals, 0, trimmedBQ, 0, toIndex);
                readBases = trimmedRB;
                baseQuals = trimmedBQ;
            }
        }

        public Cigar getCigar() {
            return (newCigar != null ? newCigar : read.getCigar());
        }

        public void setCigar(Cigar cigar) {
            setCigar(cigar, true);
        }

        // tentatively sets the new Cigar, but it needs to be confirmed later
        public void setCigar(Cigar cigar, boolean fixClippedCigar) {
            if (cigar == null) {
                newCigar = null;
                return;
            }

            if (fixClippedCigar && getReadBases().length < read.getReadLength())
                cigar = reclipCigar(cigar);

            // no change?
            if (read.getCigar().equals(cigar)) {
                newCigar = null;
                return;
            }

            // no indel?
            String str = cigar.toString();
            if (!str.contains("D") && !str.contains("I")) {

            }

            newCigar = cigar;
        }

        // pull out the bases that aren't clipped out
        private Cigar reclipCigar(Cigar cigar) {
            return IndelRealigner.reclipCigar(cigar, read);
        }

        // tentatively sets the new start, but it needs to be confirmed later
        public void setAlignmentStart(int start) {
            newStart = start;
        }

        public int getAlignmentStart() {
            return (newStart != -1 ? newStart : read.getAlignmentStart());
        }

        public int getOriginalAlignmentStart() {
            return read.getAlignmentStart();
        }

        // finalizes the changes made.
        // returns true if this record actually changes, false otherwise
        public boolean finalizeUpdate() {
            // if we haven't made any changes, don't do anything
            if (newCigar == null)
                return false;
            if (newStart == -1)
                newStart = read.getAlignmentStart();
            else if (Math.abs(newStart - read.getAlignmentStart()) > MAX_POS_MOVE_ALLOWED) {
                return false;
            }

            // store the old CIGAR and start in case we need to back out
            final Cigar oldCigar = read.getCigar();
            final int oldStart = read.getAlignmentStart();

            // try updating the read with the new CIGAR and start
            read.setCigar(newCigar);
            read.setAlignmentStart(newStart);

            // back out if necessary
            if (realignmentProducesBadAlignment(read)) {
                read.setCigar(oldCigar);
                read.setAlignmentStart(oldStart);
                return false;
            }

            // annotate the record with the original cigar and start (if it changed)
            if (!NO_ORIGINAL_ALIGNMENT_TAGS) {
                read.setAttribute(ORIGINAL_CIGAR_TAG, oldCigar.toString());
                if (newStart != oldStart)
                    read.setAttribute(ORIGINAL_POSITION_TAG, oldStart);
            }

            return true;
        }

        public void setMismatchScoreToReference(int score) {
            mismatchScoreToReference = score;
        }

        public int getMismatchScoreToReference() {
            return mismatchScoreToReference;
        }

        public void setAlignerMismatchScore(long score) {
            alignerMismatchScore = score;
        }

        public long getAlignerMismatchScore() {
            return alignerMismatchScore;
        }
    }

    /**
     * Determines whether the read aligns off the end of the contig
     *
     * @param read the read to check
     * @return true if it aligns off the end
     */
    private boolean realignmentProducesBadAlignment(final GATKSAMRecord read) {
        final int contigLength = refContentProvider.getSamSequenceDictionary()
                .getSequence(currentInterval.getContig()).getSequenceLength();
        return realignmentProducesBadAlignment(read, contigLength);
    }

    /**
     * Determines whether the read aligns off the end of the contig.
     * Pulled out to make it testable.
     *
     * @param read the read to check
     * @return true if it aligns off the end
     */
    protected static boolean realignmentProducesBadAlignment(final GATKSAMRecord read, final int contigLength) {
        return read.getAlignmentEnd() > contigLength;
    }

    @Override
    protected void onTraversalDone() {
        if (readsToClean.size() > 0) {
            GenomeLoc earliestPossibleMove = genomeLocParser.createGenomeLoc(readsToClean.getReads().get(0));
            if (manager.canMoveReads(earliestPossibleMove))
                clean(readsToClean);
            emitReadLists();
        } else if (readsNotToClean.size() > 0) {
            emitReadLists();
        }

        knownIndelsToTry.clear();
        indelRodsSeen.clear();

        manager.close();
    }

    public List<SAMRecord> getResultSam() {
        return manager.getResultRecords();
    }

    private static class Consensus {
        public final byte[] str;
        public final ArrayList<Pair<Integer, Integer>> readIndexes;
        public final int positionOnReference;
        public int mismatchSum;
        public Cigar cigar;

        public Consensus(byte[] str, Cigar cigar, int positionOnReference) {
            this.str = str;
            this.cigar = cigar;
            this.positionOnReference = positionOnReference;
            mismatchSum = 0;
            readIndexes = new ArrayList<Pair<Integer, Integer>>();
        }

        @Override
        public boolean equals(Object o) {
            return (this == o || (o instanceof Consensus && Arrays.equals(this.str, (((Consensus) o).str))));
        }

        public boolean equals(Consensus c) {
            return (this == c || Arrays.equals(this.str, c.str));
        }

        @Override
        public int hashCode() {
            return Arrays.hashCode(this.str);
        }
    }
}
