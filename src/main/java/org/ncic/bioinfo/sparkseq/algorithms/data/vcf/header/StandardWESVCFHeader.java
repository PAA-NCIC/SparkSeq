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
public class StandardWESVCFHeader {
    protected static String[] getValues() {
        String[] values = {"##fileformat=VCFv4.1",
                "##FILTER=<ID=PASS,Description=\"Accept as a confident somatic mutation\">",
                "##FILTER=<ID=REJECT,Description=\"Rejected as a confident somatic mutation\">",
                "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">",
                "##FORMAT=<ID=BQ,Number=A,Type=Float,Description=\"Average base quality for reads supporting alleles\">",
                "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">",
                "##FORMAT=<ID=FA,Number=A,Type=Float,Description=\"Allele fraction of the alternate allele with regard to reference\">",
                "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">",
                "##FORMAT=<ID=SS,Number=1,Type=Integer,Description=\"Variant status relative to non-adjacent Normal,0=wildtype,1=germline,2=somatic,3=LOH,4=post-transcriptional modification,5=unknown\">",
                "##GATKCommandLine=<ID=MuTect,Version=3.1-0-g72492bb,Date=\"Fri Dec 23 11:00:07 GMT+08:00 2016\",Epoch=1482462007722,CommandLineOptions=\"analysis_type=MuTect input_file=[Blood_realigned.bam, Tumor_realigned.bam] showFullBamList=false read_buffer_size=null phone_home=AWS gatk_key=null tag=NA read_filter=[] intervals=[0_panel.intervals] excludeIntervals=null interval_set_rule=UNION interval_merging=ALL interval_padding=100 reference_sequence=C:\\Users\\Administrator\\IdeaProjects\\GATKSrcRead\\human_g1k_v37.fasta nonDeterministicRandomSeed=false disableDithering=false maxRuntime=-1 maxRuntimeUnits=MINUTES downsampling_type=NONE downsample_to_fraction=null downsample_to_coverage=null baq=OFF baqGapOpenPenalty=40.0 fix_misencoded_quality_scores=false allow_potentially_misencoded_quality_scores=false useOriginalQualities=false defaultBaseQualities=-1 performanceLog=null BQSR=null quantize_quals=0 disable_indel_quals=false emit_original_quals=false preserve_qscores_less_than=6 globalQScorePrior=-1.0 validation_strictness=SILENT remove_program_records=false keep_program_records=false sample_rename_mapping_file=null unsafe=null disable_auto_index_creation_and_locking_when_reading_rods=false num_threads=1 num_cpu_threads_per_data_thread=1 num_io_threads=0 monitorThreadEfficiency=false num_bam_file_handles=null read_group_black_list=null pedigree=[] pedigreeString=[] pedigreeValidationType=STRICT allow_intervals_with_unindexed_bam=false generateShadowBCF=false variant_index_type=DYNAMIC_SEEK variant_index_parameter=-1 logging_level=INFO log_to_file=null help=false version=false noop=false enable_extended_output=false enable_qscore_output=false artifact_detection_mode=false tumor_sample_name=Tumor bam_tumor_sample_name=null normal_sample_name=Blood force_output=false force_alleles=false only_passing_calls=false initial_tumor_lod=4.0 tumor_lod=6.3 fraction_contamination=0.0 minimum_mutation_cell_fraction=0.0 normal_lod=2.2 normal_artifact_lod=1.0 strand_artifact_lod=2.0 strand_artifact_power_threshold=0.9 dbsnp_normal_lod=5.5 minimum_normal_allele_fraction=0.0 tumor_f_pretest=0.005 min_qscore=5 gap_events_threshold=3 heavily_clipped_read_fraction=0.3 fraction_mapq0_threshold=0.5 pir_median_threshold=10.0 pir_mad_threshold=3.0 required_maximum_alt_allele_mapping_quality_score=20 max_alt_alleles_in_normal_count=2 max_alt_alleles_in_normal_qscore_sum=20 max_alt_allele_in_normal_fraction=0.03 power_constant_qscore=30 power_constant_af=0.30000001192092896 vcf=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub no_cmdline_in_header=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub sites_only=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub bcf=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub dbsnp=[(RodBinding name=dbsnp source=C:\\Users\\Administrator\\IdeaProjects\\CheckSum\\dbsnp_138.b37.vcf)] cosmic=[] normal_panel=[] coverage_file=null coverage_20_q20_file=null power_file=null tumor_depth_file=null normal_depth_file=null filter_reads_with_N_cigar=false filter_mismatching_base_and_quals=false filter_bases_not_stored=false\">",
                "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP Membership\">",
                "##INFO=<ID=MQ0,Number=1,Type=Integer,Description=\"Total Mapping Quality Zero Reads\">",
                "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic event\">",
                "##INFO=<ID=VT,Number=1,Type=String,Description=\"Variant type, can be SNP, INS or DEL\">",
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
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNormal\tTumor"};
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
