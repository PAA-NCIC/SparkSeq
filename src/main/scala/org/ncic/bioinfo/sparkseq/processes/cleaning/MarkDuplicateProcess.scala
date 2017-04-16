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
package org.ncic.bioinfo.sparkseq.processes.cleaning

import org.apache.spark.rdd.RDD
import org.ncic.bioinfo.sparkseq.data.basic.BasicSamRecord
import org.ncic.bioinfo.sparkseq.data.bundle.SAMBundle
import org.ncic.bioinfo.sparkseq.engine.AbstractProcess
import org.ncic.bioinfo.sparkseq.exceptions.{ResourceNotSetException, ResourceSetException}
import org.ncic.bioinfo.sparkseq.resource.Resource
import org.ncic.bioinfo.sparkseq.utils.FlagUtils

import scala.collection.mutable.ListBuffer

/**
  * Author: wbc
  */
object MarkDuplicateProcess {
  def apply(name: String,
            inputSamBundle: SAMBundle,
            outputSamBundle: SAMBundle): MarkDuplicateProcess = {
    new MarkDuplicateProcess(name, inputSamBundle, outputSamBundle);
  }
}

class MarkDuplicateProcess(name: String,
                           inputSamBundle: SAMBundle,
                           outputSamBundle: SAMBundle) extends AbstractProcess(name) {

  override def getInputResourceList(): List[Resource] = List(inputSamBundle)

  override def getOutputResourceList(): List[Resource] = List(outputSamBundle)

  override def runProcess(): Unit = {
    if (!inputSamBundle.isSet) {
      throw new ResourceNotSetException(inputSamBundle.key)
    }
    if (outputSamBundle.isSet) {
      throw new ResourceSetException(outputSamBundle.key)
    }

    val inputSamRdd = inputSamBundle.samRecordRDD
    //inputSamRdd.cache()
    val dedupedSamRdd = picardMarkDuplicate(inputSamRdd)
    outputSamBundle.samRecordRDD = dedupedSamRdd
    outputSamBundle.setFlag = true
  }

  private def picardMarkDuplicate(recordRDD: RDD[BasicSamRecord]): RDD[BasicSamRecord] = {
    val supplementOrSecondaryFalg =
      FlagUtils.initialFlag
        .addFlag(FlagUtils.SECONDARY_ALIGNMENT)
        .addFlag(FlagUtils.SUPPLEMENTARY_ALIGNMENT)
        .intValue
    //val supplementReads = recordRDD.filter(f => (f.flag & supplementOrSecondaryFalg) != 0)

    val unmappedOrSecSupDupFlag =
      FlagUtils.initialFlag
        .addFlag(FlagUtils.SECONDARY_ALIGNMENT)
        .addFlag(FlagUtils.SUPPLEMENTARY_ALIGNMENT)
        .addFlag(FlagUtils.PCR_DUPLICATE)
        .addFlag(FlagUtils.SEGMENT_UNMAPPED)
        .intValue
    val dedupedRecords = recordRDD.groupBy(r => r.readName)
      .map(rg => {
        // get readEnds for a pair, may be empty
        val samRecordGroup = rg._2
        val result = ReadPairInfo()
        for (samStr <- samRecordGroup) {
          val readEnds = ReadPairInfo(samStr)
          if ((readEnds.getFlag1 & unmappedOrSecSupDupFlag) == 0) {
            result.join(readEnds)
          }
        }
        result
      })
      .groupBy(r => r.getReadPositionInfo)
      .flatMap(rg => {
        // vote the one with highest score, and make them all fragment
        var reads = ListBuffer[ReadPairInfo]()
        val codes = rg._1
        val readEndsGroup = rg._2
        if (codes._1 != 0) {
          var maxScore = -1
          var maxReadEnd: ReadPairInfo = null
          for (readEnd <- readEndsGroup) {
            if (readEnd.score > maxScore) {
              maxScore = readEnd.score
              maxReadEnd = readEnd
            }
          }
          if (maxReadEnd != null) {
            reads += ReadPairInfo(maxReadEnd.samRecord1)
            if (maxReadEnd.samRecord2 != null)
              reads += ReadPairInfo(maxReadEnd.samRecord2)
          }
        }
        reads.toList
      })
      .groupBy(r => r.getReadPositionInfo) // regroup them again to see if a fragment is with a pair
      .flatMap(rg => {
      var reads = ListBuffer[BasicSamRecord]()
      val readEndsGroup = rg._2
      val hasPair = readEndsGroup.exists(r => (r.getFlag1 & 0x08) != 0)
      for (readEnd <- readEndsGroup) {
        if (!hasPair || (readEnd.getFlag1 & 0x08) != 0)
          reads += readEnd.samRecord1
      }
      reads.toList
    })

    //dedupedRecords ++ supplementReads
    dedupedRecords
  }
}


object ReadPairInfo {

  def apply(): ReadPairInfo = {
    new ReadPairInfo
  }

  def apply(samRecord: BasicSamRecord): ReadPairInfo = {
    val readPairInfo = new ReadPairInfo()
    readPairInfo.samRecord1 = samRecord
    readPairInfo.score = samRecord.mapQ
    readPairInfo
  }
}

class ReadPairInfo() extends Serializable {
  var score = -1
  var samRecord1: BasicSamRecord = null
  var samRecord2: BasicSamRecord = null

  def getContigNo1: Int = {
    if (samRecord1 == null) 0 else samRecord1.contigId
  }

  def getContigNo2: Int = {
    if (samRecord2 == null) 0 else samRecord2.contigId
  }

  def getPosition1: Int = {
    if (samRecord1 == null) 0 else samRecord1.position
  }

  def getPosition2: Int = {
    if (samRecord2 == null) 0 else samRecord2.position
  }

  def getFlag1: Int = {
    if (samRecord1 == null) 0 else samRecord1.flag
  }

  def getFlag2: Int = {
    if (samRecord2 == null) 0 else samRecord2.flag
  }

  def join(anotherInfo: ReadPairInfo): Unit = {
    if (samRecord1 == null) {
      samRecord1 = anotherInfo.samRecord1
    }
    else {
      val contigNo = getContigNo1
      val position = getPosition1
      val newContigNo = anotherInfo.getContigNo1
      val newPosition = anotherInfo.getPosition1
      if (contigNo < newContigNo || (contigNo == newContigNo && position <= newPosition)) {
        samRecord2 = anotherInfo.samRecord1
      }
      else {
        samRecord2 = samRecord1
        samRecord1 = anotherInfo.samRecord1
      }
    }
    if (score < 0) // clean no read flag
      score = 0
    score += anotherInfo.score
  }

  def getReadPositionInfo: (Long, Long) = {
    var readCode1 = 0
    var readCode2 = 0
    if (samRecord1 != null && (getFlag1 & 3332) == 0) {
      // unmapped, duplicate, secondary, supplement
      readCode1 = getContigNo1 // contig number
      readCode1 <<= 32
      readCode1 += getPosition1 // position
      readCode1 <<= 8
      if ((getFlag1 & 0x10) > 0) // flags if negative, set 1
        readCode1 += 1
    }
    if (samRecord2 != null && (getFlag2 & 3332) == 0) {
      readCode2 = getContigNo2
      readCode2 <<= 32
      readCode2 += getPosition2
      readCode2 <<= 8
      if ((getFlag2 & 0x10) > 0)
        readCode2 += 1
    }
    (readCode1, readCode2)
  }
}