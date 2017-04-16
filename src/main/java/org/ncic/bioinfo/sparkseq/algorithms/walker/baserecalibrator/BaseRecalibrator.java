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
package org.ncic.bioinfo.sparkseq.algorithms.walker.baserecalibrator;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.tribble.Feature;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.InOriginIntervalFilter;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RODNames;
import org.ncic.bioinfo.sparkseq.algorithms.engine.ReadWalker;
import org.ncic.bioinfo.sparkseq.algorithms.utils.BaseUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.EventType;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLocParser;
import org.ncic.bioinfo.sparkseq.algorithms.utils.MathUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.ReadUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.RecalUtils;
import org.ncic.bioinfo.sparkseq.algorithms.utils.clip.ReadClipper;
import org.ncic.bioinfo.sparkseq.algorithms.data.basic.Pair;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.RefMetaDataTracker;
import org.ncic.bioinfo.sparkseq.algorithms.data.reference.ReferenceContext;
import org.ncic.bioinfo.sparkseq.algorithms.utils.reports.GATKReport;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.GATKSAMRecord;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.SamContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.DuplicateReadFilter;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.FailsVendorQualityCheckFilter;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.Filter;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.MappingQualityUnavailableFilter;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.MappingQualityZeroFilter;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.NotPrimaryAlignmentFilter;
import org.ncic.bioinfo.sparkseq.algorithms.data.sam.filter.UnmappedReadFilter;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RODContentProvider;
import org.ncic.bioinfo.sparkseq.algorithms.walker.baserecalibrator.covariate.Covariate;
import org.ncic.bioinfo.sparkseq.exceptions.ReviewedGATKException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Author: wbc
 */
public class BaseRecalibrator extends ReadWalker {

    /**
     * all the command line arguments for BQSR and it's covariates
     */
    private final RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();

    public double BAQGOP = BAQ.DEFAULT_GOP;

    /**
     * an object that keeps track of the information necessary for quality score quantization
     */
    private QuantizationInfo quantizationInfo;

    /**
     * list to hold the all the covariate objects that were requested (required + standard + experimental)
     */
    private Covariate[] requestedCovariates;

    private RecalibrationEngine recalibrationEngine;

    private int minimumQToUse = 6;

    private static final String NO_DBSNP_EXCEPTION = "This calculation is critically dependent on being able to skip over known variant sites. Please provide a VCF file containing known sites of genetic variation.";

    private BAQ baq; // BAQ the reads on the fly to generate the alignment uncertainty vector
    private final static byte NO_BAQ_UNCERTAINTY = (byte) '@';

    /**
     * result
     */
    private GATKReport report = null;

    public BaseRecalibrator(GenomeLocParser genomeLocParser,
                            RefContentProvider refContentProvider,
                            SamContentProvider samContentProvider,
                            List<RODContentProvider> rodContentProviderList) {
        super(genomeLocParser, refContentProvider, samContentProvider, rodContentProviderList);
    }

    @Override
    public void initialize() {
        baq = new BAQ(BAQGOP); // setup the BAQ object with the provided gap open penalty

        Pair<ArrayList<Covariate>, ArrayList<Covariate>> covariates = RecalUtils.initializeCovariates(RAC); // initialize the required and optional covariates
        ArrayList<Covariate> requiredCovariates = covariates.getFirst();
        ArrayList<Covariate> optionalCovariates = covariates.getSecond();

        requestedCovariates = new Covariate[requiredCovariates.size() + optionalCovariates.size()];
        int covariateIndex = 0;
        for (final Covariate covariate : requiredCovariates)
            requestedCovariates[covariateIndex++] = covariate;
        for (final Covariate covariate : optionalCovariates)
            requestedCovariates[covariateIndex++] = covariate;

        for (Covariate cov : requestedCovariates) { // list all the covariates being used
            cov.initialize(RAC); // initialize any covariate member variables using the shared argument collection
        }

        initializeRecalibrationEngine();
    }

    /**
     * Initialize the recalibration engine
     */
    private void initializeRecalibrationEngine() {
        int numReadGroups = 0;
        SAMFileHeader header = samContentProvider.getSamFileHeader();
        numReadGroups += header.getReadGroups().size();

        recalibrationEngine = new RecalibrationEngine(requestedCovariates, numReadGroups, false);
    }

    private boolean isLowQualityBase( final GATKSAMRecord read, final int offset ) {
        return read.getBaseQualities()[offset] < minimumQToUse;
    }

    @Override
    protected List<Filter> getFilter() {
        List<Filter> filters = new ArrayList<>();
        filters.add(new MappingQualityZeroFilter());
        filters.add(new MappingQualityUnavailableFilter());
        filters.add(new UnmappedReadFilter());
        filters.add(new NotPrimaryAlignmentFilter());
        filters.add(new DuplicateReadFilter());
        filters.add(new FailsVendorQualityCheckFilter());
        filters.add(new InOriginIntervalFilter(this.refContentProvider));
        return filters;
    }

    @Override
    protected void map(final ReferenceContext ref,
                       final GATKSAMRecord originalRead,
                       final RefMetaDataTracker metaDataTracker) {
        final GATKSAMRecord read = ReadClipper.hardClipSoftClippedBases( ReadClipper.hardClipAdaptorSequence(originalRead) );
        if( read.isEmpty() ) { return; } // the whole read was inside the adaptor so skip it

        RecalUtils.parsePlatformForRead(read, RAC);
        if (!RecalUtils.isColorSpaceConsistent(RAC.SOLID_NOCALL_STRATEGY, read)) { // parse the solid color space and check for color no-calls
            return; // skip this read completely
        }

        final int[] isSNP = calculateIsSNP(read, ref, originalRead);
        final int[] isInsertion = calculateIsIndel(read, EventType.BASE_INSERTION);
        final int[] isDeletion = calculateIsIndel(read, EventType.BASE_DELETION);
        final int nErrors = nEvents(isSNP, isInsertion, isDeletion);

        // note for efficiency regions we don't compute the BAQ array unless we actually have
        // some error to marginalize over.  For ILMN data ~85% of reads have no error
        final byte[] baqArray = nErrors == 0 ? flatBAQArray(read) : calculateBAQArray(read);

        if( baqArray != null ) { // some reads just can't be BAQ'ed
            final ReadCovariates covariates = RecalUtils.computeCovariates(read, requestedCovariates);
            final boolean[] skip = calculateSkipArray(read, metaDataTracker); // skip known sites of variation as well as low quality and non-regular bases
            final double[] snpErrors = calculateFractionalErrorArray(isSNP, baqArray);
            final double[] insertionErrors = calculateFractionalErrorArray(isInsertion, baqArray);
            final double[] deletionErrors = calculateFractionalErrorArray(isDeletion, baqArray);

            // aggregate all of the info into our info object, and update the data
            final ReadRecalibrationInfo info = new ReadRecalibrationInfo(read, covariates, skip, snpErrors, insertionErrors, deletionErrors);
            recalibrationEngine.updateDataForRead(info);
            return;
        } else {
            return;
        }
    }


    /**
     * Compute the number of mutational events across all hasEvent vectors
     *
     * Simply the sum of entries in hasEvents
     *
     * @param hasEvents a vector a vectors of 0 (no event) and 1 (has event)
     * @return the total number of events across all hasEvent arrays
     */
    protected static int nEvents(final int[]... hasEvents) {
        int n = 0;
        for ( final int[] hasEvent : hasEvents ) {
            n += MathUtils.sum(hasEvent);
        }
        return n;
    }

    protected boolean[] calculateSkipArray( final GATKSAMRecord read, final RefMetaDataTracker metaDataTracker ) {
        final byte[] bases = read.getReadBases();
        final boolean[] skip = new boolean[bases.length];
        final boolean[] knownSites = calculateKnownSites(read, metaDataTracker.getAllValues());
        for( int iii = 0; iii < bases.length; iii++ ) {
            skip[iii] = !BaseUtils.isRegularBase(bases[iii]) || isLowQualityBase(read, iii) || knownSites[iii] || badSolidOffset(read, iii);
        }
        return skip;
    }

    protected boolean badSolidOffset( final GATKSAMRecord read, final int offset ) {
        return ReadUtils.isSOLiDRead(read) && RAC.SOLID_RECAL_MODE != RecalUtils.SOLID_RECAL_MODE.DO_NOTHING && !RecalUtils.isColorSpaceConsistent(read, offset);
    }

    protected static boolean[] calculateKnownSites( final GATKSAMRecord read, final List<Feature> features ) {
        final int readLength = read.getReadBases().length;
        final boolean[] knownSites = new boolean[readLength];
        Arrays.fill(knownSites, false);
        for( final Feature f : features ) {
            int featureStartOnRead = ReadUtils.getReadCoordinateForReferenceCoordinate(read.getSoftStart(), read.getCigar(), f.getStart(), ReadUtils.ClippingTail.LEFT_TAIL, true); // BUGBUG: should I use LEFT_TAIL here?
            if( featureStartOnRead == ReadUtils.CLIPPING_GOAL_NOT_REACHED ) {
                featureStartOnRead = 0;
            }

            int featureEndOnRead = ReadUtils.getReadCoordinateForReferenceCoordinate(read.getSoftStart(), read.getCigar(), f.getEnd(), ReadUtils.ClippingTail.LEFT_TAIL, true);
            if( featureEndOnRead == ReadUtils.CLIPPING_GOAL_NOT_REACHED ) {
                featureEndOnRead = readLength;
            }

            if( featureStartOnRead > readLength ) {
                featureStartOnRead = featureEndOnRead = readLength;
            }

            Arrays.fill(knownSites, Math.max(0, featureStartOnRead), Math.min(readLength, featureEndOnRead + 1), true);
        }
        return knownSites;
    }

    // BUGBUG: can be merged with calculateIsIndel
    protected static int[] calculateIsSNP( final GATKSAMRecord read, final ReferenceContext ref, final GATKSAMRecord originalRead ) {
        final byte[] readBases = read.getReadBases();
        final byte[] refBases = Arrays.copyOfRange(ref.getBases(), read.getAlignmentStart() - originalRead.getAlignmentStart(), ref.getBases().length + read.getAlignmentEnd() - originalRead.getAlignmentEnd());
        final int[] snp = new int[readBases.length];
        int readPos = 0;
        int refPos = 0;
        for ( final CigarElement ce : read.getCigar().getCigarElements() ) {
            final int elementLength = ce.getLength();
            switch (ce.getOperator()) {
                case M:
                case EQ:
                case X:
                    for( int iii = 0; iii < elementLength; iii++ ) {
                        snp[readPos] = ( BaseUtils.basesAreEqual(readBases[readPos], refBases[refPos]) ? 0 : 1 );
                        readPos++;
                        refPos++;
                    }
                    break;
                case D:
                case N:
                    refPos += elementLength;
                    break;
                case I:
                case S: // ReferenceContext doesn't have the soft clipped bases!
                    readPos += elementLength;
                    break;
                case H:
                case P:
                    break;
                default:
                    throw new ReviewedGATKException("Unsupported cigar operator: " + ce.getOperator());
            }
        }
        return snp;
    }

    protected static int[] calculateIsIndel( final GATKSAMRecord read, final EventType mode ) {
        final int[] indel = new int[read.getReadBases().length];
        int readPos = 0;
        for ( final CigarElement ce : read.getCigar().getCigarElements() ) {
            final int elementLength = ce.getLength();
            switch (ce.getOperator()) {
                case M:
                case EQ:
                case X:
                case S:
                {
                    readPos += elementLength;
                    break;
                }
                case D:
                {
                    final int index = ( read.getReadNegativeStrandFlag() ? readPos : readPos - 1 );
                    updateIndel(indel, index, mode, EventType.BASE_DELETION);
                    break;
                }
                case I:
                {
                    final boolean forwardStrandRead = !read.getReadNegativeStrandFlag();
                    if( forwardStrandRead ) {
                        updateIndel(indel, readPos - 1, mode, EventType.BASE_INSERTION);
                    }
                    readPos += elementLength;
                    if( !forwardStrandRead ) {
                        updateIndel(indel, readPos, mode, EventType.BASE_INSERTION);
                    }
                    break;
                }
                case N:
                case H:
                case P:
                    break;
                default:
                    throw new ReviewedGATKException("Unsupported cigar operator: " + ce.getOperator());
            }
        }
        return indel;
    }

    private static void updateIndel(final int[] indel, final int index, final EventType mode, final EventType requiredMode) {
        if ( mode == requiredMode && index >= 0 && index < indel.length )
            // protect ourselves from events at the start or end of the read (1D3M or 3M1D)
            indel[index] = 1;
    }

    protected static double[] calculateFractionalErrorArray( final int[] errorArray, final byte[] baqArray ) {
        if(errorArray.length != baqArray.length ) {
            throw new ReviewedGATKException("Array length mismatch detected. Malformed read?");
        }

        final int BLOCK_START_UNSET = -1;

        final double[] fractionalErrors = new double[baqArray.length];
        Arrays.fill(fractionalErrors, 0.0);
        boolean inBlock = false;
        int blockStartIndex = BLOCK_START_UNSET;
        int iii;
        for( iii = 0; iii < fractionalErrors.length; iii++ ) {
            if( baqArray[iii] == NO_BAQ_UNCERTAINTY ) {
                if( !inBlock ) {
                    fractionalErrors[iii] = (double) errorArray[iii];
                } else {
                    calculateAndStoreErrorsInBlock(iii, blockStartIndex, errorArray, fractionalErrors);
                    inBlock = false; // reset state variables
                    blockStartIndex = BLOCK_START_UNSET; // reset state variables
                }
            } else {
                inBlock = true;
                if( blockStartIndex == BLOCK_START_UNSET ) { blockStartIndex = iii; }
            }
        }
        if( inBlock ) {
            calculateAndStoreErrorsInBlock(iii-1, blockStartIndex, errorArray, fractionalErrors);
        }
        if( fractionalErrors.length != errorArray.length ) {
            throw new ReviewedGATKException("Output array length mismatch detected. Malformed read?");
        }
        return fractionalErrors;
    }

    private static void calculateAndStoreErrorsInBlock( final int iii,
                                                        final int blockStartIndex,
                                                        final int[] errorArray,
                                                        final double[] fractionalErrors ) {
        int totalErrors = 0;
        for( int jjj = Math.max(0,blockStartIndex-1); jjj <= iii; jjj++ ) {
            totalErrors += errorArray[jjj];
        }
        for( int jjj = Math.max(0, blockStartIndex-1); jjj <= iii; jjj++ ) {
            fractionalErrors[jjj] = ((double) totalErrors) / ((double)(iii - Math.max(0,blockStartIndex-1) + 1));
        }
    }

    /**
     * Create a BAQ style array that indicates no alignment uncertainty
     * @param read the read for which we want a BAQ array
     * @return a BAQ-style non-null byte[] counting NO_BAQ_UNCERTAINTY values
     * // TODO -- could be optimized avoiding this function entirely by using this inline if the calculation code above
     */
    protected  static byte[] flatBAQArray(final GATKSAMRecord read) {
        final byte[] baq = new byte[read.getReadLength()];
        Arrays.fill(baq, NO_BAQ_UNCERTAINTY);
        return baq;
    }

    /**
     * Compute an actual BAQ array for read, based on its quals and the reference sequence
     * @param read the read to BAQ
     * @return a non-null BAQ tag array for read
     */
    private byte[] calculateBAQArray( final GATKSAMRecord read ) {
        baq.baqRead(read, refContentProvider, BAQ.CalculationMode.RECALCULATE, BAQ.QualityMode.ADD_TAG);
        return BAQ.getBAQTag(read);
    }

    @Override
    protected void onTraversalDone() {
        recalibrationEngine.finalizeData();
        RecalibrationTables recalibrationTables = recalibrationEngine.getFinalRecalibrationTables();
        quantizationInfo = new QuantizationInfo(
                recalibrationTables , RAC.QUANTIZING_LEVELS);
        report = RecalUtils.getRecalibrationReport(RAC, quantizationInfo,
                recalibrationTables, requestedCovariates, RAC.SORT_BY_ALL_COLUMNS);
    }

    public GATKReport getReport() {
        return report;
    }
}
