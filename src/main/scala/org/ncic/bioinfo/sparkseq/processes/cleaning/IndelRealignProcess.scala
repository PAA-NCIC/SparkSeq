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
import org.ncic.bioinfo.sparkseq.data.common.RefPartitionInfo
import org.ncic.bioinfo.sparkseq.data.partition.{BundlePartition, FastaPartition, SamRecordPartition, VcfRecordPartition}
import org.ncic.bioinfo.sparkseq.algorithms.adapters.IndelRealignAdapter
import org.apache.spark.SparkContext._
import org.ncic.bioinfo.sparkseq.debug.Dumper
import org.ncic.bioinfo.sparkseq.exceptions.PipelineException

import collection.JavaConversions._
import scala.collection.mutable.ListBuffer

/**
  * Author: wbc
  */
object IndelRealignProcess {
  def apply(name: String,
            referencePath: String,
            rodMap: Map[String, String],
            refPartitionInfoBundle: RefPartitionInfoBundle,
            inputSamBundleList: List[SAMBundle],
            outputSamBundleList: List[SAMBundle]): IndelRealignProcess = {
    if (inputSamBundleList.size != outputSamBundleList.size) {
      throw new PipelineException("Count of input and output samples must be equal")
    }
    val inputSamBundleMap = inputSamBundleList.map(samBundle => (samBundle.key, samBundle)).toMap
    val outputSamBundleMap = outputSamBundleList.map(samBundle => (samBundle.key, samBundle)).toMap

    val inOutSampleMap = inputSamBundleList.zip(outputSamBundleList)
      .map(bundlePair => (bundlePair._1.key, bundlePair._2.key)).toMap

    new IndelRealignProcess(name, referencePath, rodMap, refPartitionInfoBundle,
      inputSamBundleMap, outputSamBundleMap, inOutSampleMap)
  }
}

class IndelRealignProcess(name: String,
                          referencePath: String,
                          rodMap: Map[String, String],
                          refPartitionInfoBundle: RefPartitionInfoBundle,
                          inputSamBundleMap: Map[String, SAMBundle],
                          outputSamBundleMap: Map[String, SAMBundle],
                          inOutSampleMap: Map[String, String])
  extends DataCleanProcess(name, referencePath, rodMap, refPartitionInfoBundle, inputSamBundleMap, outputSamBundleMap) {

  override def getCleanedBundleRDD(bundleRDD: RDD[BundlePartition]): RDD[BundlePartition] = {
    val refPartitionInfo = refPartitionInfoBundle.refPartitionInfo
    val rodKeysBD = sc.broadcast(rodMap.keys).value
    val inOutSampleMapBD = sc.broadcast(inOutSampleMap).value

    bundleRDD.map(bundle => {
      try {
        val partitionId = bundle.partitionId
        var resultSAMPartitionMap = Map[String, SamRecordPartition]()
        val rodList = rodKeysBD.map(key => bundle.rodPartitionMap.get(key).get).toList

        bundle.samRecordPartitionMap.foreach(samRecordPartitionWithKey => {
          val key = samRecordPartitionWithKey._1
          val samRecordPartition = samRecordPartitionWithKey._2
          val resultRecords = IndelRealignAdapter.realign(
            bundle.refContigInfo, samRecordPartition, bundle.fastaPartition, rodList)

          val resultSAMPartiton = new SamRecordPartition(partitionId, samRecordPartition.contigId,
            resultRecords, samRecordPartition.samHeaderInfo)
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
