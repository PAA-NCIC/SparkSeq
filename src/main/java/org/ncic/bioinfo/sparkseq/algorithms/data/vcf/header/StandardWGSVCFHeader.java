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
package org.ncic.bioinfo.sparkseq.algorithms.data.vcf.header;

import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.VCFHeaderLineIterable;

/**
 * Author: wbc
 */
public class StandardWGSVCFHeader {
    private static String[] getValues() {
        String[] values = {"##fileformat=VCFv4.1",
                "##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele at this location\">",
                "##FILTER=<ID=LowQual,Description=\"Low quality\">",
                "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">",
                "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">",
                "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum DP observed within the GVCF block\">",
                "##FORMAT=<ID=PGT,Number=1,Type=String,Description=\"Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another\">",
                "##FORMAT=<ID=PID,Number=1,Type=String,Description=\"Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group\">",
                "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">",
                "##FORMAT=<ID=SB,Number=4,Type=Integer,Description=\"Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.\">",
                "##GATKCommandLine=<ID=GenotypeGVCFs,Version=3.3-0-g37228af,Date=\"Tue Dec 06 14:09:55 GMT+08:00 2016\",Epoch=1481004595938,CommandLineOptions=\"analysis_type=GenotypeGVCFs input_file=[] showFullBamList=false read_buffer_size=null phone_home=AWS gatk_key=null tag=NA read_filter=[] intervals=null excludeIntervals=null interval_set_rule=UNION interval_merging=ALL interval_padding=0 reference_sequence=human_g1k_v37.fasta nonDeterministicRandomSeed=false disableDithering=false maxRuntime=-1 maxRuntimeUnits=MINUTES downsampling_type=BY_SAMPLE downsample_to_fraction=null downsample_to_coverage=1000 baq=OFF baqGapOpenPenalty=40.0 refactor_NDN_cigar_string=false fix_misencoded_quality_scores=false allow_potentially_misencoded_quality_scores=false useOriginalQualities=false defaultBaseQualities=-1 performanceLog=null BQSR=null quantize_quals=0 disable_indel_quals=false emit_original_quals=false preserve_qscores_less_than=6 globalQScorePrior=-1.0 validation_strictness=SILENT remove_program_records=false keep_program_records=false sample_rename_mapping_file=null unsafe=null disable_auto_index_creation_and_locking_when_reading_rods=false no_cmdline_in_header=false sites_only=false never_trim_vcf_format_field=false bcf=false bam_compression=null simplifyBAM=false disable_bam_indexing=false generate_md5=false num_threads=1 num_cpu_threads_per_data_thread=1 num_io_threads=0 monitorThreadEfficiency=false num_bam_file_handles=null read_group_black_list=null pedigree=[] pedigreeString=[] pedigreeValidationType=STRICT allow_intervals_with_unindexed_bam=false generateShadowBCF=false variant_index_type=DYNAMIC_SEEK variant_index_parameter=-1 logging_level=INFO log_to_file=null help=false version=false variant=[(RodBindingCollection [(RodBinding name=variant source=gvcf.vcf)])] out=org.broadinstitute.gatk.engine.io.stubs.VariantContextWriterStub includeNonVariantSites=false annotateNDA=false heterozygosity=0.001 indel_heterozygosity=1.25E-4 standard_min_confidence_threshold_for_calling=30.0 standard_min_confidence_threshold_for_emitting=30.0 max_alternate_alleles=6 input_prior=[] sample_ploidy=2 annotation=[InbreedingCoeff, FisherStrand, QualByDepth, ChromosomeCounts, GenotypeSummaries, StrandOddsRatio] dbsnp=(RodBinding name= source=UNBOUND) filter_reads_with_N_cigar=false filter_mismatching_base_and_quals=false filter_bases_not_stored=false\">",
                "##GATKCommandLine=<ID=HaplotypeCaller,Version=3.3-0-g37228af,Date=\"Tue Nov 29 14:26:02 GMT+08:00 2016\",Epoch=1480400762872,CommandLineOptions=\"analysis_type=HaplotypeCaller input_file=[recal_reads.bam] showFullBamList=false read_buffer_size=null phone_home=AWS gatk_key=null tag=NA read_filter=[] intervals=[intervaled.intervals] excludeIntervals=null interval_set_rule=UNION interval_merging=ALL interval_padding=0 reference_sequence=human_g1k_v37.fasta nonDeterministicRandomSeed=false disableDithering=false maxRuntime=-1 maxRuntimeUnits=MINUTES downsampling_type=BY_SAMPLE downsample_to_fraction=null downsample_to_coverage=250 baq=OFF baqGapOpenPenalty=40.0 refactor_NDN_cigar_string=false fix_misencoded_quality_scores=false allow_potentially_misencoded_quality_scores=false useOriginalQualities=false defaultBaseQualities=-1 performanceLog=null BQSR=null quantize_quals=0 disable_indel_quals=false emit_original_quals=false preserve_qscores_less_than=6 globalQScorePrior=-1.0 validation_strictness=SILENT remove_program_records=false keep_program_records=false sample_rename_mapping_file=null unsafe=null disable_auto_index_creation_and_locking_when_reading_rods=false no_cmdline_in_header=false sites_only=false never_trim_vcf_format_field=false bcf=false bam_compression=null simplifyBAM=false disable_bam_indexing=false generate_md5=false num_threads=1 num_cpu_threads_per_data_thread=1 num_io_threads=0 monitorThreadEfficiency=false num_bam_file_handles=null read_group_black_list=null pedigree=[] pedigreeString=[] pedigreeValidationType=STRICT allow_intervals_with_unindexed_bam=false generateShadowBCF=false variant_index_type=LINEAR variant_index_parameter=128000 logging_level=INFO log_to_file=null help=false version=false likelihoodCalculationEngine=PairHMM heterogeneousKmerSizeResolution=COMBO_MIN graphOutput=null bamOutput=null bamWriterType=CALLED_HAPLOTYPES disableOptimizations=false dbsnp=(RodBinding name= source=UNBOUND) dontTrimActiveRegions=false maxDiscARExtension=25 maxGGAARExtension=300 paddingAroundIndels=150 paddingAroundSNPs=20 comp=[] annotation=[ClippingRankSumTest, DepthPerSampleHC, StrandBiasBySample] excludeAnnotation=[SpanningDeletions, TandemRepeatAnnotator, ChromosomeCounts, FisherStrand, StrandOddsRatio, QualByDepth] debug=false useFilteredReadsForAnnotations=false emitRefConfidence=GVCF annotateNDA=false heterozygosity=0.001 indel_heterozygosity=1.25E-4 standard_min_confidence_threshold_for_calling=-0.0 standard_min_confidence_threshold_for_emitting=-0.0 max_alternate_alleles=6 input_prior=[] sample_ploidy=2 genotyping_mode=DISCOVERY alleles=(RodBinding name= source=UNBOUND) contamination_fraction_to_filter=0.0 contamination_fraction_per_sample_file=null p_nonref_model=null exactcallslog=null output_mode=EMIT_VARIANTS_ONLY allSitePLs=true sample_name=null kmerSize=[10, 25] dontIncreaseKmerSizesForCycles=false allowNonUniqueKmersInRef=false numPruningSamples=1 recoverDanglingHeads=false doNotRecoverDanglingBranches=false minDanglingBranchLength=4 consensus=false GVCFGQBands=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 70, 80, 90, 99] indelSizeToEliminateInRefModel=10 min_base_quality_score=10 minPruning=2 gcpHMM=10 includeUmappedReads=false useAllelesTrigger=false phredScaledGlobalReadMismappingRate=45 maxNumHaplotypesInPopulation=128 mergeVariantsViaLD=false doNotRunPhysicalPhasing=false pair_hmm_implementation=VECTOR_LOGLESS_CACHING keepRG=null justDetermineActiveRegions=false dontGenotype=false errorCorrectKmers=false debugGraphTransformations=false dontUseSoftClippedBases=false captureAssemblyFailureBAM=false allowCyclesInKmerGraphToGeneratePaths=false noFpga=false errorCorrectReads=false kmerLengthForReadErrorCorrection=25 minObservationsForKmerToBeSolid=20 pcr_indel_model=CONSERVATIVE maxReadsInRegionPerSample=1000 minReadsPerAlignmentStart=5 activityProfileOut=null activeRegionOut=null activeRegionIn=null activeRegionExtension=null forceActive=false activeRegionMaxSize=null bandPassSigma=null maxProbPropagationDistance=50 activeProbabilityThreshold=0.002 min_mapping_quality_score=20 filter_reads_with_N_cigar=false filter_mismatching_base_and_quals=false filter_bases_not_stored=false\">",
                "##GVCFBlock=minGQ=0(inclusive),maxGQ=1(exclusive)",
                "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">",
                "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">",
                "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">",
                "##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities\">",
                "##INFO=<ID=CCC,Number=1,Type=Integer,Description=\"Number of called chromosomes\">",
                "##INFO=<ID=ClippingRankSum,Number=1,Type=Float,Description=\"Z-score From Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases\">",
                "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">",
                "##INFO=<ID=DS,Number=0,Type=Flag,Description=\"Were any of the samples downsampled?\">",
                "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">",
                "##INFO=<ID=FS,Number=1,Type=Float,Description=\"Phred-scaled p-value using Fisher's exact test to detect strand bias\">",
                "##INFO=<ID=GQ_MEAN,Number=1,Type=Float,Description=\"Mean of all GQ values\">",
                "##INFO=<ID=GQ_STDDEV,Number=1,Type=Float,Description=\"Standard deviation of all GQ values\">",
                "##INFO=<ID=HWP,Number=1,Type=Float,Description=\"P value from test of Hardy Weinberg Equilibrium\">",
                "##INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description=\"Consistency of the site with at most two segregating haplotypes\">",
                "##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description=\"Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation\">",
                "##INFO=<ID=MLEAC,Number=A,Type=Integer,Description=\"Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed\">",
                "##INFO=<ID=MLEAF,Number=A,Type=Float,Description=\"Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed\">",
                "##INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMS Mapping Quality\">",
                "##INFO=<ID=MQ0,Number=1,Type=Integer,Description=\"Total Mapping Quality Zero Reads\">",
                "##INFO=<ID=MQRankSum,Number=1,Type=Float,Description=\"Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities\">",
                "##INFO=<ID=NCC,Number=1,Type=Integer,Description=\"Number of no-called samples\">",
                "##INFO=<ID=QD,Number=1,Type=Float,Description=\"Variant Confidence/Quality by Depth\">",
                "##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias\">",
                "##INFO=<ID=SOR,Number=1,Type=Float,Description=\"Symmetric Odds Ratio of 2x2 contingency table to detect strand bias\">",
                "##contig=<ID=1,length=249250621,assembly=b37>",
                "##contig=<ID=2,length=243199373,assembly=b37>",
                "##contig=<ID=3,length=198022430,assembly=b37>",
                "##contig=<ID=4,length=191154276,assembly=b37>",
                "##contig=<ID=5,length=180915260,assembly=b37>",
                "##contig=<ID=6,length=171115067,assembly=b37>",
                "##contig=<ID=7,length=159138663,assembly=b37>",
                "##contig=<ID=8,length=146364022,assembly=b37>",
                "##contig=<ID=9,length=141213431,assembly=b37>",
                "##contig=<ID=10,length=135534747,assembly=b37>",
                "##contig=<ID=11,length=135006516,assembly=b37>",
                "##contig=<ID=12,length=133851895,assembly=b37>",
                "##contig=<ID=13,length=115169878,assembly=b37>",
                "##contig=<ID=14,length=107349540,assembly=b37>",
                "##contig=<ID=15,length=102531392,assembly=b37>",
                "##contig=<ID=16,length=90354753,assembly=b37>",
                "##contig=<ID=17,length=81195210,assembly=b37>",
                "##contig=<ID=18,length=78077248,assembly=b37>",
                "##contig=<ID=19,length=59128983,assembly=b37>",
                "##contig=<ID=20,length=63025520,assembly=b37>",
                "##contig=<ID=21,length=48129895,assembly=b37>",
                "##contig=<ID=22,length=51304566,assembly=b37>",
                "##contig=<ID=X,length=155270560,assembly=b37>",
                "##contig=<ID=Y,length=59373566,assembly=b37>",
                "##contig=<ID=MT,length=16569,assembly=b37>",
                "##contig=<ID=GL000207.1,length=4262,assembly=b37>",
                "##contig=<ID=GL000226.1,length=15008,assembly=b37>",
                "##contig=<ID=GL000229.1,length=19913,assembly=b37>",
                "##contig=<ID=GL000231.1,length=27386,assembly=b37>",
                "##contig=<ID=GL000210.1,length=27682,assembly=b37>",
                "##contig=<ID=GL000239.1,length=33824,assembly=b37>",
                "##contig=<ID=GL000235.1,length=34474,assembly=b37>",
                "##contig=<ID=GL000201.1,length=36148,assembly=b37>",
                "##contig=<ID=GL000247.1,length=36422,assembly=b37>",
                "##contig=<ID=GL000245.1,length=36651,assembly=b37>",
                "##contig=<ID=GL000197.1,length=37175,assembly=b37>",
                "##contig=<ID=GL000203.1,length=37498,assembly=b37>",
                "##contig=<ID=GL000246.1,length=38154,assembly=b37>",
                "##contig=<ID=GL000249.1,length=38502,assembly=b37>",
                "##contig=<ID=GL000196.1,length=38914,assembly=b37>",
                "##contig=<ID=GL000248.1,length=39786,assembly=b37>",
                "##contig=<ID=GL000244.1,length=39929,assembly=b37>",
                "##contig=<ID=GL000238.1,length=39939,assembly=b37>",
                "##contig=<ID=GL000202.1,length=40103,assembly=b37>",
                "##contig=<ID=GL000234.1,length=40531,assembly=b37>",
                "##contig=<ID=GL000232.1,length=40652,assembly=b37>",
                "##contig=<ID=GL000206.1,length=41001,assembly=b37>",
                "##contig=<ID=GL000240.1,length=41933,assembly=b37>",
                "##contig=<ID=GL000236.1,length=41934,assembly=b37>",
                "##contig=<ID=GL000241.1,length=42152,assembly=b37>",
                "##contig=<ID=GL000243.1,length=43341,assembly=b37>",
                "##contig=<ID=GL000242.1,length=43523,assembly=b37>",
                "##contig=<ID=GL000230.1,length=43691,assembly=b37>",
                "##contig=<ID=GL000237.1,length=45867,assembly=b37>",
                "##contig=<ID=GL000233.1,length=45941,assembly=b37>",
                "##contig=<ID=GL000204.1,length=81310,assembly=b37>",
                "##contig=<ID=GL000198.1,length=90085,assembly=b37>",
                "##contig=<ID=GL000208.1,length=92689,assembly=b37>",
                "##contig=<ID=GL000191.1,length=106433,assembly=b37>",
                "##contig=<ID=GL000227.1,length=128374,assembly=b37>",
                "##contig=<ID=GL000228.1,length=129120,assembly=b37>",
                "##contig=<ID=GL000214.1,length=137718,assembly=b37>",
                "##contig=<ID=GL000221.1,length=155397,assembly=b37>",
                "##contig=<ID=GL000209.1,length=159169,assembly=b37>",
                "##contig=<ID=GL000218.1,length=161147,assembly=b37>",
                "##contig=<ID=GL000220.1,length=161802,assembly=b37>",
                "##contig=<ID=GL000213.1,length=164239,assembly=b37>",
                "##contig=<ID=GL000211.1,length=166566,assembly=b37>",
                "##contig=<ID=GL000199.1,length=169874,assembly=b37>",
                "##contig=<ID=GL000217.1,length=172149,assembly=b37>",
                "##contig=<ID=GL000216.1,length=172294,assembly=b37>",
                "##contig=<ID=GL000215.1,length=172545,assembly=b37>",
                "##contig=<ID=GL000205.1,length=174588,assembly=b37>",
                "##contig=<ID=GL000219.1,length=179198,assembly=b37>",
                "##contig=<ID=GL000224.1,length=179693,assembly=b37>",
                "##contig=<ID=GL000223.1,length=180455,assembly=b37>",
                "##contig=<ID=GL000195.1,length=182896,assembly=b37>",
                "##contig=<ID=GL000212.1,length=186858,assembly=b37>",
                "##contig=<ID=GL000222.1,length=186861,assembly=b37>",
                "##contig=<ID=GL000200.1,length=187035,assembly=b37>",
                "##contig=<ID=GL000193.1,length=189789,assembly=b37>",
                "##contig=<ID=GL000194.1,length=191469,assembly=b37>",
                "##contig=<ID=GL000225.1,length=211173,assembly=b37>",
                "##contig=<ID=GL000192.1,length=547496,assembly=b37>",
                "##reference=FAKE_REFERENCE",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1"};
        return values;
    }

    private static VCFHeader header = null;

    public static VCFHeader getHeader() {
        if(header == null) {
            VCFCodec codec = new VCFCodec();
            VCFHeaderLineIterable headerLineIterable = new VCFHeaderLineIterable(getValues());
            header = (VCFHeader) codec.readActualHeader(headerLineIterable);
        }
        return header;
    }
}
