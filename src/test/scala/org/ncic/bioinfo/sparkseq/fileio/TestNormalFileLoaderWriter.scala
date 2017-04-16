package org.ncic.bioinfo.sparkseq.fileio

import org.junit.runner.RunWith
import org.ncic.bioinfo.sparkseq.data.common.{RefContigInfo, SamHeaderInfo, VcfHeaderInfo}
import org.scalatest.FunSuite
import org.scalatest.junit.JUnitRunner

/**
  * Author: wbc
  */
@RunWith(classOf[JUnitRunner])
class TestNormalFileLoaderWriter extends FunSuite {

  test("Test read / write fastq file") {
    val filePath = getClass().getResource("/test1.fastq").getFile()
    val fastqRecords = NormalFileLoader.loadFastq(filePath)
    assert(fastqRecords.size == 2500)
    assert(fastqRecords.head.quality.equals("AAA#ACCCCCHGFHHHHHHGFHHBHFGEGEHEHAFHHFFAECBCEEFFFEG:F<E<;@;>BDCB66>>79E@EDC7<,=@#####################"))
    assert(fastqRecords.last.sequence.equals("CAGGAAATGTCGATGTCCAAGCGGGCCTCAAAGATGGACTTGCCATTGTTGATGCACTCCATAGTAGCAATTTCATCCTCCCGCTCCTTTGGGTAAAGCAG"))

    val resultFilePath = "test_result/test_result.fastq"
    NormalFileWriter.writeFastq(fastqRecords, resultFilePath)
    val resultRecords = NormalFileLoader.loadFastq(resultFilePath)
    assert(resultRecords.size == 2500)
    assert(resultRecords.head.quality.equals("AAA#ACCCCCHGFHHHHHHGFHHBHFGEGEHEHAFHHFFAECBCEEFFFEG:F<E<;@;>BDCB66>>79E@EDC7<,=@#####################"))
    assert(resultRecords.last.sequence.equals("CAGGAAATGTCGATGTCCAAGCGGGCCTCAAAGATGGACTTGCCATTGTTGATGCACTCCATAGTAGCAATTTCATCCTCCCGCTCCTTTGGGTAAAGCAG"))
  }

  test("Test read / write sam file") {
    val filePath = getClass().getResource("/test.sam").getFile()
    val refContigInfo = RefContigInfo(getClass().getResource("/human_g1k_v37.dict").getFile())
    val samRecords = NormalFileLoader.loadSam(filePath, refContigInfo)
    assert(samRecords.size == 5021)

    val resultFilePath = "test_result/test_result.sam"
    val headerInfo = SamHeaderInfo.sortedHeader(refContigInfo, List("@ for test"))
    NormalFileWriter.writeSam(headerInfo, samRecords, resultFilePath)

    val resultRecords = NormalFileLoader.loadSam(resultFilePath, refContigInfo)
    assert(resultRecords.size == 5021)
  }

  test("Test read / write vcf file") {
    val filePath = getClass().getResource("/test.vcf").getFile()
    val refContigInfo = RefContigInfo(getClass().getResource("/human_g1k_v37.dict").getFile())
    val vcfRecords = NormalFileLoader.loadVcf(filePath, refContigInfo)
    assert(vcfRecords.size == 963)

    val resultFilePath = "test_result/test_result.vcf"
    val headerInfo = VcfHeaderInfo.newHeader(refContigInfo, List("# for test"))
    NormalFileWriter.writeVcf(headerInfo, vcfRecords, resultFilePath)

    val resultRecords = NormalFileLoader.loadVcf(resultFilePath, refContigInfo)
    assert(resultRecords.size == 963)
  }

  test("Test read / write fastq file pair") {
    val filePath1 = getClass().getResource("/test1.fastq").getFile()
    val filePath2 = getClass().getResource("/test2.fastq").getFile()
    val fastqRecords = NormalFileLoader.loadFastqPair(filePath1, filePath2)
    assert(fastqRecords.size == 2500)

    val resultFilePath1 = "test_result/test_result1.fastq"
    val resultFilePath2 = "test_result/test_result2.fastq"
    NormalFileWriter.writeFastqPair(fastqRecords, resultFilePath1, resultFilePath2)
    val resultRecords = NormalFileLoader.loadFastq(resultFilePath1)
    assert(resultRecords.size == 2500)
    assert(resultRecords.head.quality.equals("AAA#ACCCCCHGFHHHHHHGFHHBHFGEGEHEHAFHHFFAECBCEEFFFEG:F<E<;@;>BDCB66>>79E@EDC7<,=@#####################"))
    assert(resultRecords.last.sequence.equals("CAGGAAATGTCGATGTCCAAGCGGGCCTCAAAGATGGACTTGCCATTGTTGATGCACTCCATAGTAGCAATTTCATCCTCCCGCTCCTTTGGGTAAAGCAG"))
  }
}
