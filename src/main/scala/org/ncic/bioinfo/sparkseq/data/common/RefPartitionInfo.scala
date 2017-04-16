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

import org.ncic.bioinfo.sparkseq.const.{BinTools, PipelineConst, SamRecordConst}
import org.ncic.bioinfo.sparkseq.exceptions.{IllegalInputException, PipelineException}

import scala.collection.mutable
import scala.collection.mutable.ListBuffer

/**
  * Author: wbc
  */
object RefPartitionInfo {
  /**
    * 将reference按照scatterLen均切，如果一个contig的边缘或者整条contig都不足长度，则视为一条scatter
    *
    * @param refContigInfo reference上contig的id、name和长度信息
    * @param scatterLen    每一段scatter的长度
    * @return Reference切分后的Partition的信息
    */
  def apply(refContigInfo: RefContigInfo, scatterLen: Int): RefPartitionInfo = {

    val contigIds = refContigInfo.getContigIds
    var scatterIdIdx = 0
    val contigPartitionNoStartMap = mutable.HashMap[Int, Int]()
    val partitionRangeList = ListBuffer[(Int, Int, Int)]()
    contigIds.foreach(contigId => {
      val contigLen = refContigInfo.getLength(contigId)
      val scatterCountInContig = (contigLen + scatterLen - 1) / scatterLen
      val partitionRangesInContig = Range(0, scatterCountInContig)
        .map(idx => (idx * scatterLen + 1, Math.min(idx * scatterLen + scatterLen, contigLen), contigId))

      contigPartitionNoStartMap.put(contigId, scatterIdIdx)
      partitionRangeList.appendAll(partitionRangesInContig)
      scatterIdIdx += scatterCountInContig
    })
    new RefPartitionInfo(refContigInfo, scatterIdIdx, scatterLen,
      contigPartitionNoStartMap.toMap, partitionRangeList.toList, Set())
  }

  /**
    * 将reference按照默认的scatterLen均切
    *
    * @param refContigInfo reference上contig的id、name和长度信息
    * @return Reference切分后的Partition的信息
    */
  def apply(refContigInfo: RefContigInfo): RefPartitionInfo = {
    apply(refContigInfo, BinTools.DEFAULT_PARTITION_LENGTH)
  }

  def apply(old: RefPartitionInfo, skipPartitionIds: Set[Int]): RefPartitionInfo = {
    new RefPartitionInfo(old.refContigInfo, old.partitionCount, old.scatterLen,
      old.contigPartitionNoStartMap, old.partitionRangeList, skipPartitionIds)
  }
}

/**
  * 这个类和overlap没有任何关系，所以所有返回的值和判断逻辑都没有考虑过overlap
  */
class RefPartitionInfo(val refContigInfo: RefContigInfo,
                       val partitionCount: Int,
                       val scatterLen: Int,
                       val contigPartitionNoStartMap: Map[Int, Int],
                       val partitionRangeList: List[(Int, Int, Int)],
                       val skipPartitionIds: Set[Int]) extends Serializable {

  def getPartitionId(contigId: Int, position: Int): Int = {
    val scatterStartId = contigPartitionNoStartMap.getOrElse(contigId, -1)
    if (scatterStartId < 0) {
      return SamRecordConst.FAKE_PARTITION_ID
    }
    // 如果scatterLen是10000，那么表示从[1,10000]是一个partition
    val rawId = scatterStartId + (position - 1) / scatterLen
    if(skipPartitionIds.contains(rawId)) {
      return SamRecordConst.FAKE_PARTITION_ID
    } else {
      return rawId
    }
  }

  /**
    * 获取partition在contig中的起止position，以1为基
    * 当不存在对应partition时返回(-1,-1,FAKE_CONTIG_ID)
    *
    * @param partitionId
    * @return (start, end, contigId)
    */
  def getPartitionRange(partitionId: Int): (Int, Int, Int) = {
    if (partitionId < 0 || partitionId >= partitionCount) {
      (-1, -1, SamRecordConst.FAKE_CONTIG_ID)
    } else {
      partitionRangeList(partitionId)
    }
  }

  /**
    * 获得一个contig中的partition的range列表
    * 当不存在指定contigId时抛出异常
    *
    * @param contigId
    * @return (partitionId, start, end)的数组
    */
  def getPartitionRangesInContig(contigId: Int): List[(Int, Int, Int)] = {
    val startPartitionId = contigPartitionNoStartMap.getOrElse(contigId, -1)
    val len = refContigInfo.getLength(contigId)
    if (len == 0) {
      throw new IllegalInputException("Can't find contig with id:" + contigId)
    }
    val partitionCount = (len + scatterLen - 1) / scatterLen
    Range(0, partitionCount)
      .map(id => {
        val endTmp = id * scatterLen + scatterLen
        val end = if (endTmp > len) len else endTmp
        (startPartitionId + id, id * scatterLen + 1, end)
      }).toList
  }

  /**
    * 获得一个contig中的partition的range列表
    * 同时如果指定带有overlap，那么将range向两端扩展。
    * 如果overlap大于1/2 scatter的长度，那么实际的overlap是1/2 scatter的长度。
    * 当不存在指定contigId时抛出异常
    *
    * @param contigId
    * @return (partitionId, start, end)的数组
    */
  def getPartitionRangesInContig(contigId: Int, overlapLen: Int): List[(Int, Int, Int)] = {
    val finalOverlapLen = if (overlapLen > scatterLen / 2) scatterLen / 2 else overlapLen
    val len = refContigInfo.getLength(contigId)
    val originRanges = getPartitionRangesInContig(contigId)
    originRanges.map(range => {
      val partitionId = range._1
      val start = if (range._2 - finalOverlapLen < 1) 1 else range._2 - finalOverlapLen
      val end = if (range._3 + finalOverlapLen > len) len else range._3 + finalOverlapLen
      (partitionId, start, end)
    })
  }

  def getPartitionCount = partitionCount

  def getScatterLen = scatterLen

  def getRefContigInfo = refContigInfo
}
