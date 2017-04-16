package org.ncic.bioinfo.sparkseq.data

import org.junit.runner.RunWith
import org.ncic.bioinfo.sparkseq.data.common.RefContigInfo
import org.scalatest.FunSuite
import org.scalatest.junit.JUnitRunner

/**
  * Author: wbc
  */
@RunWith(classOf[JUnitRunner])
class TestRefContigInfo extends FunSuite {

  test("Test ref contig info load") {
    val refContigInfo = RefContigInfo(getClass().getResource("/human_g1k_v37.dict").getFile())
    assert(refContigInfo.getId("GL000207.1") == 25)
    assert(refContigInfo.getLength("GL000207.1") == 4262)
  }
}
