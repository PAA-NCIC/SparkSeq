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
import org.ncic.bioinfo.sparkseq.data.bundle.{RefPartitionInfoBundle, SAMBundle}
import org.ncic.bioinfo.sparkseq.data.common.{RefPartitionInfo, SamHeaderInfo}
import org.ncic.bioinfo.sparkseq.data.partition.{BundlePartition, SamRecordPartition}
import org.ncic.bioinfo.sparkseq.engine.{PartitionConsumer, PartitionGenerator}
import org.ncic.bioinfo.sparkseq.processes.PartitionProcess
import org.ncic.bioinfo.sparkseq.resource.Resource

import scala.collection.mutable.ListBuffer

/**
  * Author: wbc
  */
abstract class DataCleanProcess(name: String,
                                referencePath: String,
                                rodMap: Map[String, String],
                                refPartitionInfoBundle: RefPartitionInfoBundle,
                                inputSamBundleMap: Map[String, SAMBundle],
                                outputSamBundleMap: Map[String, SAMBundle])
  extends PartitionProcess(name, referencePath, rodMap, refPartitionInfoBundle, inputSamBundleMap)
    with PartitionGenerator {

  override def getInputResourceList(): List[Resource] = {
    val inputList = ListBuffer[Resource]()
    inputList ++= inputSamBundleMap.values
    inputList.append(refPartitionInfoBundle)
    inputList.toList
  }

  override def getOutputResourceList(): List[Resource] = outputSamBundleMap.values.toList

  override def generateBundlePartition(bundleRDD: RDD[BundlePartition]): RDD[BundlePartition] = {
    getCleanedBundleRDD(bundleRDD)
  }

  override def consumeBundlePartition(bundleRDD: RDD[BundlePartition]): Unit = {
    val resultBundleRDD = getCleanedBundleRDD(bundleRDD)
    val refContigInfoValue = sc.broadcast(refPartitionInfoBundle.refPartitionInfo.getRefContigInfo).value

    // 写入结果
    val samMapRDD = resultBundleRDD.map(bundle => bundle.samRecordPartitionMap)
    outputSamBundleMap.foreach(samBundleWithKey => {
      val key = samBundleWithKey._1
      val samBundle = samBundleWithKey._2
      val keyBD = sc.broadcast(key).value
      val samRDD = resultBundleRDD.flatMap(bundle => {
        val refPartition = bundle.fastaPartition
        val samPartition = bundle.samRecordPartitionMap.get(keyBD)
        samPartition.getOrElse(SamRecordPartition.empty(refPartition.partitionId,  refPartition.contigId, SamHeaderInfo.sortedHeader(refContigInfoValue, null)))
          .records.filter(record =>
          (record.position >= refPartition.originStart && record.position <= refPartition.originEnd))
      })
      samBundle.setFlag = true
      samBundle.samRecordRDD = samRDD
    })
  }

  def getCleanedBundleRDD(bundleRDD: RDD[BundlePartition]): RDD[BundlePartition]
}
