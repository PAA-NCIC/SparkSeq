package org.ncic.bioinfo.sparkseq.consts

import org.junit.runner.RunWith
import org.ncic.bioinfo.sparkseq.const.BinTools
import org.ncic.bioinfo.sparkseq.utils.FileUtils
import org.scalatest.FunSuite
import org.scalatest.junit.JUnitRunner

/**
  * Author: wbc
  */
@RunWith(classOf[JUnitRunner])
class TestConst extends FunSuite {

  test("Test bin path") {
    assert(FileUtils.isExists(BinTools.bwaPath))
    assert(FileUtils.isExists(BinTools.confPath))
  }
}
