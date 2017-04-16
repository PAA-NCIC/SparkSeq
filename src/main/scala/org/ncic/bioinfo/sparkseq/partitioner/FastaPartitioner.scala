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
package org.ncic.bioinfo.sparkseq.partitioner

import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD
import org.ncic.bioinfo.sparkseq.const.SamRecordConst
import org.ncic.bioinfo.sparkseq.data.common.{RefContigInfo, RefPartitionInfo}
import org.ncic.bioinfo.sparkseq.data.partition.FastaPartition
import org.ncic.bioinfo.sparkseq.exceptions.IllegalInputException
import org.ncic.bioinfo.sparkseq.utils.StringUtils

import scala.collection.mutable.ListBuffer
import scala.io.Source

/**
  * Author: wbc
  */
object FastaPartitioner {

  def SAFE_OVERLAP_LEN = 1000

  def partition(sc: SparkContext, refPartitionInfo: RefPartitionInfo,
                referencePath: String): RDD[(Int, FastaPartition)] = {
    partition(sc, refPartitionInfo, referencePath, 0)
  }

  def partition(sc: SparkContext, refPartitionInfo: RefPartitionInfo,
                referencePath: String, overlapLen: Int): RDD[(Int, FastaPartition)] = {
    val refContigInfo = refPartitionInfo.getRefContigInfo
    val lines = Source.fromFile(referencePath).getLines()
    val fastaContigMap = parseFasta(lines, refContigInfo).toMap

    val fastaPartitions = getFastaPartitions(fastaContigMap, refPartitionInfo, overlapLen)
    sc.makeRDD(fastaPartitions)
  }

  def getFastaPartitions(fastaContigMap: Map[Int, String],
                         refPartitionInfo: RefPartitionInfo,
                         overlapLen: Int): ListBuffer[(Int, FastaPartition)] = {
    val refContigInfo = refPartitionInfo.getRefContigInfo
    val contigIds = refContigInfo.getContigIds
    val fastaPartitions = ListBuffer[(Int, FastaPartition)]()
    for (contigId <- contigIds) {
      val contigContent: String = fastaContigMap.get(contigId).get
      val fastaPartitionsInContigTmp =
        refPartitionInfo.getPartitionRangesInContig(contigId)

      val fastaPartitionsInContig = fastaPartitionsInContigTmp
        .map(range => {
          val partitionId = range._1
          val startCoordinate = range._2
          val endCoordiname = range._3
          (partitionId, FastaPartition(
            partitionId, contigId, refContigInfo.getName(contigId),
            contigContent, startCoordinate, endCoordiname, overlapLen, SAFE_OVERLAP_LEN))
        })
      fastaPartitions ++= fastaPartitionsInContig
    }
    fastaPartitions
  }

  /**
    * 假设fasta中的contig不会超过int的表示范围，实际上int有21亿的表示范围，最长的contig只有2亿多
    *
    * @param lineIter
    * @param refContigInfo
    * @return
    */
  def parseFasta(lineIter: Iterator[String],
                 refContigInfo: RefContigInfo): List[(Int, String)] = {
    val contigs = ListBuffer[(Int, String)]()
    val buffer = new StringBuilder()
    var curContigId = SamRecordConst.FAKE_CONTIG_ID
    for (line <- lineIter) {
      if (line.startsWith(">")) {
        // description line
        if (curContigId != SamRecordConst.FAKE_CONTIG_ID) {
          contigs.append((curContigId, buffer.toString()))
        }
        val contigName = StringUtils.split(line, " ").next().substring(1)
        curContigId = refContigInfo.getId(contigName)
        if (curContigId == SamRecordConst.FAKE_CONTIG_ID) {
          throw new IllegalInputException("Error when parse fasta file, can't find contig name in dict")
        }
        buffer.clear()
      } else {
        buffer.append(line)
      }
    }
    if (curContigId != SamRecordConst.FAKE_CONTIG_ID) {
      contigs.append((curContigId, buffer.toString()))
    }
    contigs.toList
  }
}
