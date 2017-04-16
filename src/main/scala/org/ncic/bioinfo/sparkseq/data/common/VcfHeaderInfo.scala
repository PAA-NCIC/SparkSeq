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
package org.ncic.bioinfo.sparkseq.data.common

import scala.collection.mutable.ListBuffer

/**
  * Author: wbc
  */
object VcfHeaderInfo extends Serializable {

  val constHeaderLines = List(
    "##fileformat=VCFv4.1",
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
    "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">",
    "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">",
    "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">",
    "##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities\">",
    "##INFO=<ID=CCC,Number=1,Type=Integer,Description=\"Number of called chromosomes\">",
    "##INFO=<ID=ClippingRankSum,Number=1,Type=Float,Description=\"Z-score From Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases\">",
    "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP Membership\">",
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
    "##INFO=<ID=SOR,Number=1,Type=Float,Description=\"Symmetric Odds Ratio of 2x2 contingency table to detect strand bias\">"
  )

  def newHeader(refContigInfo: RefContigInfo, initExtLines: List[String]): VcfHeaderInfo = {
    val contigLines = refContigInfo.getContigIds.map(contigId => {
      val name = refContigInfo.getName(contigId)
      val length = refContigInfo.getLength(contigId)
      getContigHeaderline(name, length, refContigInfo.getRefType)
    }).toList
    var extLines = List[String]()
    if (initExtLines != null) {
      extLines = initExtLines
    }
    new VcfHeaderInfo(constHeaderLines, contigLines, extLines)
  }

  /**
    * 用于在外部读入的header lines，不区分hd和sq，全假装是hdlines。
    *
    * @param allLines
    * @return
    */
  def newHeader(allLines: List[String]): VcfHeaderInfo = {
    new VcfHeaderInfo(allLines, List(), List[String]())
  }

  private def getContigHeaderline(contigName: String, contigLength: Int, refType: String): String =
    "##contig=<ID=%s,length=%s,assembly=%s>".format(contigName, contigLength, refType)
}

class VcfHeaderInfo(hdLines: List[String],
                    sqLines: List[String],
                    extLines: List[String]) extends Serializable {

  def getHeaderLines: List[String] = {
    val buffer = ListBuffer[String]()
    buffer ++= hdLines
    buffer ++= sqLines
    buffer ++= extLines
    buffer += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    buffer.toList
  }
}
