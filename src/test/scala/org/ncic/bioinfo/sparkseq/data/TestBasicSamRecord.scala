package org.ncic.bioinfo.sparkseq.data

import java.io.File

import org.junit.runner.RunWith
import org.ncic.bioinfo.sparkseq.data.basic.BasicSamRecord
import org.ncic.bioinfo.sparkseq.data.common.RefContigInfo
import org.scalatest.FunSuite
import org.scalatest.junit.JUnitRunner

import scala.io.Source

/**
  * Author: wbc
  */
@RunWith(classOf[JUnitRunner])
class TestBasicSamRecord extends FunSuite {

  test("Test basic sam record read and write") {

    val refContigInfo = RefContigInfo(getClass().getResource("/human_g1k_v37.dict").getFile())
    val lines = Source.fromFile(new File(getClass().getResource("/test.sam").getFile()))
      .getLines()
      .filter(line => !line.startsWith("@"))
    val records = lines.map(line => BasicSamRecord(line, refContigInfo, false)).toList

    assert(records.size == 5021)

    val record = records.head
    assert(record.readName.equals("SRR504516.1492_HWI-ST423_0087:3:1:19251:2214"))
    assert(record.flag == 99)
    assert(record.contigName.equals("1"))
    assert(record.position == 69774)
    assert(record.toString().equals("SRR504516.1492_HWI-ST423_0087:3:1:19251:2214\t99\t1\t69774\t0\t101M\t=\t69944\t271\tAGCTCTGTCCACTTTGACTGCTCACATTACAGTAGTTCTTTTGTTCTTTGGACCATTTGTCTTTATTTATGCCTGGCCATTCCACATCAAGTCATTAGATA\tHHHHHHHHHHHHHHHHHHGHHHGGFHHHGHGGHDHEGEHGGEGDGEFEEE@FE@ED?=@DBDDBBEE##################################\tNM:i:2\tMD:Z:56G26C17\tAS:i:91\tXS:i:91\tRG:Z:SRR504516\tXA:Z:15,-102462479,101M,2;19,+111362,101M,2;"))
  }
}
