package org.ncic.bioinfo.sparkseq.data

import org.junit.runner.RunWith
import org.ncic.bioinfo.sparkseq.data.common.{RefContigInfo, SamHeaderInfo, VcfHeaderInfo}
import org.scalatest.FunSuite
import org.scalatest.junit.JUnitRunner

/**
  * Author: wbc
  */
@RunWith(classOf[JUnitRunner])
class TestHeaderInfo extends FunSuite {

  test("Test generate sam header info") {
    val refContigInfo = RefContigInfo(getClass().getResource("/human_g1k_v37.dict").getFile())
    val headerInfo = SamHeaderInfo.sortedHeader(refContigInfo, List[String]())
    val lines = headerInfo.getHeaderLines
  }

  test("Test generate vcf header info") {
    val refContigInfo = RefContigInfo(getClass().getResource("/human_g1k_v37.dict").getFile())
    val headerInfo = VcfHeaderInfo.newHeader(refContigInfo, List[String]())
    val lines = headerInfo.getHeaderLines
    val len = lines.length
  }

}
