package org.ncic.bioinfo.sparkseq.data

import org.junit.runner.RunWith
import org.ncic.bioinfo.sparkseq.data.common.{RefContigInfo, RefPartitionInfo}
import org.ncic.bioinfo.sparkseq.partitioner.FastaPartitioner
import org.scalatest.FunSuite
import org.scalatest.junit.JUnitRunner

import scala.io.Source

/**
  * Author: wbc
  */
@RunWith(classOf[JUnitRunner])
class TestFastaPartitioner extends FunSuite {

  test("Test fasta contig reader") {
    val referencePath = getClass().getResource("/littleFasta.fasta").getFile()
    val refContigInfo = RefContigInfo(getClass().getResource("/littleFasta.dict").getFile())
    val refPartitionInfo = RefPartitionInfo(refContigInfo, 10000)

    val scatterLen = refPartitionInfo.getScatterLen
    val lines = Source.fromFile(referencePath).getLines()

    val contigList = FastaPartitioner.parseFasta(lines, refContigInfo)
    assert(contigList.length == 59)
    assert(contigList(0)._1 == refContigInfo.getId("GL000207.1"))
    assert(contigList(0)._2.length == refContigInfo.getLength("GL000207.1"))
    assert(contigList(58)._1 == refContigInfo.getId("GL000192.1"))
    assert(contigList(58)._2.length == refContigInfo.getLength("GL000192.1"))
  }

  test("Test fasta partitioner") {
    val referencePath = getClass().getResource("/littleFasta.fasta").getFile()
    val refContigInfo = RefContigInfo(getClass().getResource("/littleFasta.dict").getFile())
    val refPartitionInfo = RefPartitionInfo(refContigInfo, 100)
    val scatterLen = refPartitionInfo.getScatterLen
    val lines = Source.fromFile(referencePath).getLines()
    val fastaContigMap = FastaPartitioner.parseFasta(lines, refContigInfo).toMap

    val partitionList = FastaPartitioner.getFastaPartitions(fastaContigMap, refPartitionInfo, 10)
    assert(partitionList(0)._2.partitionId == 0)
    assert(partitionList(0)._2.contigId == 0)
    assert(partitionList(0)._2.content.length == 110 + FastaPartitioner.SAFE_OVERLAP_LEN)
    assert(partitionList(1)._2.content.length == 90 + 10 + 100 + 10 + FastaPartitioner.SAFE_OVERLAP_LEN)
  }
}
