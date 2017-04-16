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

import org.apache.spark.SparkContext._
import org.apache.spark.rdd.RDD
import org.ncic.bioinfo.sparkseq.algorithms.adapters.IndelRealignAdapter
import org.ncic.bioinfo.sparkseq.const.SamRecordConst
import org.ncic.bioinfo.sparkseq.data.basic.BasicSamRecord
import org.ncic.bioinfo.sparkseq.data.bundle.{RefPartitionInfoBundle, SAMBundle}
import org.ncic.bioinfo.sparkseq.data.common.{RefContigInfo, RefPartitionInfo}
import org.ncic.bioinfo.sparkseq.data.partition.{BundlePartition, SamRecordPartition}
import org.ncic.bioinfo.sparkseq.debug.Dumper
import org.ncic.bioinfo.sparkseq.exceptions.PipelineException
import org.ncic.bioinfo.sparkseq.partitioner.SamRecordPartitioner
import org.ncic.bioinfo.sparkseq.utils.FlagUtils

import scala.collection.mutable.ListBuffer

/**
  * Author: wbc
  */
object PartitionMarkDuplicateProcess {
  def apply(name: String,
            referencePath: String,
            refPartitionInfoBundle: RefPartitionInfoBundle,
            inputSamBundleList: List[SAMBundle],
            outputSamBundleList: List[SAMBundle]): PartitionMarkDuplicateProcess = {
    if (inputSamBundleList.size != outputSamBundleList.size) {
      throw new PipelineException("Count of input and output samples must be equal")
    }
    val inputSamBundleMap = inputSamBundleList.map(samBundle => (samBundle.key, samBundle)).toMap
    val outputSamBundleMap = outputSamBundleList.map(samBundle => (samBundle.key, samBundle)).toMap

    val inOutSampleMap = inputSamBundleList.zip(outputSamBundleList)
      .map(bundlePair => (bundlePair._1.key, bundlePair._2.key)).toMap

    new PartitionMarkDuplicateProcess(name, referencePath, refPartitionInfoBundle,
      inputSamBundleMap, outputSamBundleMap, inOutSampleMap)
  }
}

class PartitionMarkDuplicateProcess(name: String,
                                    referencePath: String,
                                    refPartitionInfoBundle: RefPartitionInfoBundle,
                                    inputSamBundleMap: Map[String, SAMBundle],
                                    outputSamBundleMap: Map[String, SAMBundle],
                                    inOutSampleMap: Map[String, String])
  extends DataCleanProcess(name, referencePath, Map[String, String](), refPartitionInfoBundle, inputSamBundleMap, outputSamBundleMap) {

  override def getCleanedBundleRDD(bundleRDD: RDD[BundlePartition]): RDD[BundlePartition] = {
    val refPartitionInfo = refPartitionInfoBundle.refPartitionInfo
    val inOutSampleMapBD = sc.broadcast(inOutSampleMap).value
    val refPartitionInfoBD = sc.broadcast(refPartitionInfo).value
    val refContigInfoBD = sc.broadcast(refPartitionInfo.getRefContigInfo).value
    bundleRDD.map(bundle => {
      try {
        val partitionId = bundle.partitionId
        var resultSAMPartitionMap = Map[String, SamRecordPartition]()

        bundle.samRecordPartitionMap.foreach(samRecordPartitionWithKey => {
          val key = samRecordPartitionWithKey._1
          val samRecordPartition = samRecordPartitionWithKey._2
          val records = samRecordPartition.records

          // 带有这种flag的最后需要保留
          val supplementOrSecondaryFalg =
            FlagUtils.initialFlag
              .addFlag(FlagUtils.SECONDARY_ALIGNMENT)
              .addFlag(FlagUtils.SUPPLEMENTARY_ALIGNMENT)
              .intValue

          // 不参与dedup的read带有的flag
          val unmappedOrSecSupDupFlag =
            FlagUtils.initialFlag
              .addFlag(FlagUtils.SECONDARY_ALIGNMENT)
              .addFlag(FlagUtils.SUPPLEMENTARY_ALIGNMENT)
              .addFlag(FlagUtils.PCR_DUPLICATE)
              .addFlag(FlagUtils.SEGMENT_UNMAPPED)
              .intValue

          val supplementReads = records.filter(f => (f.flag & supplementOrSecondaryFalg) != 0)
          val dedupedRecords = records.groupBy(r => r.readName)
            .map(rg => {
              // get readEnds for a pair, may be empty
              val samRecordGroup = rg._2
              val result = ReadPairInfo2()
              for (samStr <- samRecordGroup) {
                if ((samStr.flag & unmappedOrSecSupDupFlag) == 0) {
                  result.join(ReadPairInfo2(samStr))
                }
              }
              result
            })
            .groupBy(r => r.getPairReadPositionInfo(refContigInfoBD))
            .flatMap(rg => {
              // vote the one with highest score, and make them all fragment
              var reads = ListBuffer[ReadPairInfo2]()
              val codes = rg._1
              val readEndsGroup = rg._2
              if (codes._1 != 0) {
                var maxScore = -1
                var maxReadEnd: ReadPairInfo2 = null
                for (readEnd <- readEndsGroup) {
                  if (readEnd.score > maxScore) {
                    maxScore = readEnd.score
                    maxReadEnd = readEnd
                  }
                }
                if (maxReadEnd != null) {
                  reads += ReadPairInfo2(maxReadEnd.samRecord1)
                  if (maxReadEnd.samRecord2 != null)
                    reads += ReadPairInfo2(maxReadEnd.samRecord2)
                }
              }
              reads.toList
            })
            .groupBy(r => r.getSingleReadPositionInfo) // regroup them again to see if a fragment is with a pair
            .flatMap(rg => {
            var reads = ListBuffer[BasicSamRecord]()
            val readEndsGroup = rg._2
            val hasPair = readEndsGroup.exists(r => (r.getFlag1 & 0x08) != 0)
            for (readEnd <- readEndsGroup) {
              if (!hasPair || (readEnd.getFlag1 & 0x08) != 0)
                reads += readEnd.samRecord1
            }
            reads
          })

          val finalPartitionReads = (dedupedRecords ++ supplementReads)
            .filter(record => refPartitionInfoBD.getPartitionId(record.contigId, record.position) == partitionId)
            .toList

          val resultSAMPartiton = new SamRecordPartition(partitionId, samRecordPartition.contigId,
            finalPartitionReads, samRecordPartition.samHeaderInfo)
          resultSAMPartitionMap += (inOutSampleMapBD(key) -> resultSAMPartiton)
        })

        new BundlePartition(partitionId, bundle.refContigInfo, bundle.fastaPartition, resultSAMPartitionMap, bundle.rodPartitionMap)
      } catch {
        case e: Exception => {
          Dumper.dumpBundle(bundle, Dumper.defaultString)
          throw e
        }
      }
    })
  }
}

object ReadPairInfo2 {

  def apply(): ReadPairInfo2 = {
    new ReadPairInfo2
  }

  def apply(samRecord: BasicSamRecord): ReadPairInfo2 = {
    val readPairInfo = new ReadPairInfo2()
    readPairInfo.samRecord1 = samRecord
    readPairInfo.score = samRecord.mapQ
    readPairInfo
  }
}

class ReadPairInfo2() extends Serializable {
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

  def join(anotherInfo: ReadPairInfo2): Unit = {
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

  def getPairReadPositionInfo(refContigInfo: RefContigInfo): (Long, Long) = {
    var readCode1 = 0
    var readCode2 = 0
    // 如果samRecord1有效，且samRecord2 map上了却没有出现在readPairInfo中，
    // 则认为是被划在了其他的partition中。
    // 此时生成的code按照samRecord1中对自己和mate的信息，质量只计算samRecord1中的质量
    if(samRecord1 != null && (getFlag1 & 3332) == 0
      && samRecord2 == null && (getFlag1 & 0x8) == 0) {
      // read本身
      var code1 = getContigNo1 // contig number
      code1 <<= 32
      code1 += getPosition1 // position
      code1 <<= 8
      // flags if negative, set 1
      if ((getFlag1 & 0x10) > 0) {
        code1 += 1
      }
      // mate read
      val mateContigId = samRecord1.mateContigId
      val matePosition = samRecord1.matePosition
      var code2 = mateContigId
      code2 <<= 32
      code2 += matePosition
      code2 <<= 8
      if ((getFlag1 & 0x20) > 0){
        code2 += 1
      }
      // 实际上顺序不重要
      if(getContigNo1 < mateContigId
        || (getContigNo1 == mateContigId && getPosition1 <= matePosition)) {
        readCode1 = code1
        readCode2 = code2
      } else {
        readCode1 = code2
        readCode2 = code1
      }
    } else {
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
    }

    (readCode1, readCode2)
  }

  def getSingleReadPositionInfo: Long = {
    var readCode1 = 0
    if (samRecord1 != null && (getFlag1 & 3332) == 0) {
      readCode1 = getContigNo1 // contig number
      readCode1 <<= 32
      readCode1 += getPosition1 // position
      readCode1 <<= 8
      if ((getFlag1 & 0x10) > 0) // flags if negative, set 1
        readCode1 += 1
    }
    readCode1
  }
}
