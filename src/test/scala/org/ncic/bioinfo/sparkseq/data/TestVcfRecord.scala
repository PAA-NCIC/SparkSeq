package org.ncic.bioinfo.sparkseq.data

import java.io.File

import org.junit.runner.RunWith
import org.ncic.bioinfo.sparkseq.data.basic.VcfRecord
import org.ncic.bioinfo.sparkseq.data.common.RefContigInfo
import org.scalatest.FunSuite
import org.scalatest.junit.JUnitRunner

import scala.io.Source

/**
  * Author: wbc
  */
@RunWith(classOf[JUnitRunner])
class TestVcfRecord extends FunSuite {

  test("Test vcf record read and write") {
    val refContigInfo = RefContigInfo(getClass().getResource("/human_g1k_v37.dict").getFile())
    val lines = Source.fromFile(new File(getClass().getResource("/test.vcf").getFile()))
      .getLines()
      .filter(line => !line.startsWith("#"))
    val records = lines.map(line => VcfRecord(line, refContigInfo)).toList

    assert(records.size == 963)

    val record = records.head
    assert(record.contigName.equals("6"))
    assert(record.position == 25100945)
    assert(record.id.equals("."))
    assert(record.ref.equals("T"))
    assert(record.alt.equals("A"))
    assert(record.qual == 56.28)
    assert(record.filter.equals("."))
    assert(record.toString().equals("6\t25100945\t.\tT\tA\t56.28\t.\tAC=2;AF=1.00;AN=2;DP=3;FS=0.000;GQ_MEAN=9.00;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;NCC=0;QD=18.76;SOR=1.179\tGT:AD:DP:GQ:PL\t1/1:0,3:3:9:84,9,0"))
  }
}
