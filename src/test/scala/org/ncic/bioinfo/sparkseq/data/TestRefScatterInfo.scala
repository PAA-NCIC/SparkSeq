package org.ncic.bioinfo.sparkseq.data

import org.junit.runner.RunWith
import org.ncic.bioinfo.sparkseq.const.{BinTools, PipelineConst}
import org.ncic.bioinfo.sparkseq.data.common.{RefContigInfo, RefPartitionInfo}
import org.scalatest.FunSuite
import org.scalatest.junit.JUnitRunner

/**
  * Author: wbc
  */
@RunWith(classOf[JUnitRunner])
class TestRefScatterInfo extends FunSuite {

  test("Test ref contig scatter short") {
    val refContigInfo = RefContigInfo(getClass().getResource("/human_g1k_v37.dict").getFile())
    assert(refContigInfo.getLength("1") == 249250621)
    val scattered = RefPartitionInfo(refContigInfo, 1000)
    val id = refContigInfo.getId("2")
    assert(scattered.getPartitionId(id, 1) == 249251)
    assert(scattered.getPartitionId(id, 10000) == 249260)
    assert(scattered.getPartitionId(id, 10001) == 249261)
  }

  test("Test ref contig scatter long") {
    val refContigInfo = RefContigInfo(getClass().getResource("/human_g1k_v37.dict").getFile())
    assert(refContigInfo.getLength("1") == 249250621)
    val scattered = RefPartitionInfo(refContigInfo, 1249250621)
    val id = refContigInfo.getId("GL000192.1")
    assert(scattered.getPartitionId(id, 666) == 83)
  }

  test("Test get range") {
    val refContigInfo = RefContigInfo(getClass().getResource("/human_g1k_v37.dict").getFile())
    assert(refContigInfo.getLength("1") == 249250621)
    val scattered = RefPartitionInfo(refContigInfo, 10000)
    assert(scattered.getPartitionRange(0)._1 == 1)
    assert(scattered.getPartitionRange(0)._2 == 10000)
    assert(scattered.getPartitionRange(24925)._1 == 249250001)
    assert(scattered.getPartitionRange(24925)._2 == 249250621)
    assert(scattered.getPartitionRange(2492600)._2 == -1)
  }

  test("Test partition range list") {
    val refContigInfo = RefContigInfo(getClass().getResource("/human_g1k_v37.dict").getFile())
    assert(refContigInfo.getLength("1") == 249250621)
    val scattered = RefPartitionInfo(refContigInfo, 10000)
    val rangeList = scattered.getPartitionRangesInContig(0)
    assert(rangeList.length == 24926)
    assert(rangeList(0)._2 == 1)
    assert(rangeList(0)._3 == 10000)
    assert(rangeList(24925)._2 == 249250001)
    assert(rangeList(24925)._3 == 249250621)
  }

  test("Test partition range list with overlap") {
    val refContigInfo = RefContigInfo(getClass().getResource("/human_g1k_v37.dict").getFile())
    assert(refContigInfo.getLength("1") == 249250621)
    val scattered = RefPartitionInfo(refContigInfo, 10000)
    val rangeList = scattered.getPartitionRangesInContig(0, 10)
    assert(rangeList.length == 24926)
    assert(rangeList(0)._2 == 1)
    assert(rangeList(0)._3 == 10010)
    assert(rangeList(24925)._2 == 249249991)
    assert(rangeList(24925)._3 == 249250621)
  }

  test("Test partition range list test") {
    val refContigInfo = RefContigInfo(getClass().getResource("/human_g1k_v37.dict").getFile())
    assert(refContigInfo.getLength("1") == 249250621)
    val scattered = RefPartitionInfo(refContigInfo, BinTools.DEFAULT_PARTITION_LENGTH)
    val range = scattered.getPartitionRange(3172)

    print(scattered.getPartitionId(6, 14027864))
  }
}
