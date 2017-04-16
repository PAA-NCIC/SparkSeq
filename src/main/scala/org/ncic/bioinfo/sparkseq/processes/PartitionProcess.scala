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
package org.ncic.bioinfo.sparkseq.processes

import org.apache.spark.SparkContext._
import org.apache.spark.rdd.RDD
import org.ncic.bioinfo.sparkseq.const.{PipelineConst, ResourceKeys}
import org.ncic.bioinfo.sparkseq.data.bundle.{FASTAPartitionBundle, RefPartitionInfoBundle, SAMBundle, VCFPartitionBundle}
import org.ncic.bioinfo.sparkseq.data.common.RefPartitionInfo
import org.ncic.bioinfo.sparkseq.data.partition.{BundlePartition, FastaPartition, SamRecordPartition, VcfRecordPartition}
import org.ncic.bioinfo.sparkseq.engine.{PartitionConsumer, PartitionGenerator, PartitionOptimizedProcess}
import org.ncic.bioinfo.sparkseq.exceptions.ResourceNotSetException
import org.ncic.bioinfo.sparkseq.fileio.NormalFileLoader
import org.ncic.bioinfo.sparkseq.partitioner.{FastaPartitioner, SamRecordPartitioner, VcfPartitioner}

/**
  * Author: wbc
  */
abstract class PartitionProcess(name: String,
                                referencePath: String,
                                val rodMap: Map[String, String],
                                refPartitionInfoBundle: RefPartitionInfoBundle,
                                inputSamBundleMap: Map[String, SAMBundle]) extends PartitionOptimizedProcess(name) {

  protected def getReferencePartition(overlapLength: Int): RDD[(Int, FastaPartition)] = {
    val refKey = ResourceKeys.REFERENCE_BUNDLE + referencePath
    val refPartitionInfo = refPartitionInfoBundle.refPartitionInfo
    if (resourcePool.containsResourceKey(refKey)) {
      val fastaBundle = resourcePool.getResourceByKey(refKey)
        .get.asInstanceOf[FASTAPartitionBundle]
      fastaBundle.fastaPartitionRDD
    } else {
      val fastaPartition = FastaPartitioner
        .partition(sc, refPartitionInfo, referencePath, overlapLength)
      val refContigInfo = refPartitionInfo.getRefContigInfo
      val fastaBundle = new FASTAPartitionBundle(refKey, refContigInfo, fastaPartition)
      resourcePool.addResource(fastaBundle)
      fastaPartition
    }
  }

  protected def getRODMap(): Map[String, String] = {
    // 如果是chain的开头，则需要向后收集所有的ROD
    if (this.isInstanceOf[PartitionGenerator] && this.asInstanceOf[PartitionGenerator].setAsChainGenerator) {
      var allMap = Map[String, String]()
      var processPtr = this.asInstanceOf[PartitionConsumer]
      while (processPtr != null && processPtr.isInstanceOf[PartitionGenerator]) {
        allMap ++= processPtr.asInstanceOf[PartitionProcess].rodMap
        processPtr = processPtr.asInstanceOf[PartitionGenerator].childProcess
      }
      if(processPtr != null) {
        allMap ++= processPtr.asInstanceOf[PartitionProcess].rodMap
      }

      allMap
    } else {
      rodMap
    }
  }

  protected def getRODPartitions(overlapLength: Int): RDD[(Int, Map[String, VcfRecordPartition])] = {
    val refPartitionInfo = refPartitionInfoBundle.refPartitionInfo
    val refContigInfo = refPartitionInfo.getRefContigInfo
    var rodPartitionListRDD = sc.makeRDD(Range(0, refPartitionInfo.getPartitionCount)
      .map(partitionId => (partitionId, Map[String, VcfRecordPartition]())))
    getRODMap().foreach(bundle => {
      val key = bundle._1
      val path = bundle._2

      // 优先从resource pool中获取rdd
      val vcfHeader = NormalFileLoader.loadVcfHeader(path)
      val resourceKey = path
      val vcfRecordPartitionPairRDD = if (resourcePool.containsResourceKey(resourceKey)) {
        val vcfPartitionBundle = resourcePool.getResourceByKey(resourceKey).asInstanceOf[VCFPartitionBundle]
        vcfPartitionBundle.vcfPartitionPairRDD
      } else {
        val vcfRDD = NormalFileLoader.loadVcfToRdd(sc, path, refContigInfo)
        VcfPartitioner.partition(sc, refPartitionInfo, key, vcfHeader, vcfRDD, overlapLength)
      }

      // 将rod join成一个RDD
      val rodKeyBD = sc.broadcast(key).value
      val vcfHeaderBD = sc.broadcast(vcfHeader).value
      val refPartitionInfoBD = sc.broadcast(refPartitionInfo).value
      rodPartitionListRDD = rodPartitionListRDD.leftOuterJoin(vcfRecordPartitionPairRDD)
        .map(bundle => {
          val partitionId = bundle._1
          var rodMap = bundle._2._1
          val partitionOption = bundle._2._2
          rodMap += (rodKeyBD -> partitionOption.getOrElse(
            VcfRecordPartition.empty(partitionId, rodKeyBD, refPartitionInfoBD.getPartitionRange(partitionId)._3, vcfHeaderBD)))
          (partitionId, rodMap)
        })
    })
    rodPartitionListRDD
  }

  protected def getSamPartitions(overlapLength: Int): RDD[(Int, Map[String, SamRecordPartition])] = {
    //var samRecordPartitionListRDD = sc.makeRDD(Range(0, refPartitionInfo.getPartitionCount)
    //  .map(partitionId => (partitionId, Map[String, SamRecordPartition]())))

    val refPartitionInfo = refPartitionInfoBundle.refPartitionInfo
    var samRecordPartitionListRDD:RDD[(Int, Map[String, SamRecordPartition])] = null

    inputSamBundleMap.foreach(samBundleWithSampleName => {
      val sampleName = samBundleWithSampleName._1
      val samBundle = samBundleWithSampleName._2
      val inputSamRdd = samBundle.samRecordRDD
      val headerInfo = samBundle.samHeaderInfo
      headerInfo.sorted = true
      val samRecordPartitionRDD = if (overlapLength != 0) {
        SamRecordPartitioner.partition(sc, refPartitionInfo, headerInfo, inputSamRdd, overlapLength)
      } else {
        SamRecordPartitioner.partition(sc, refPartitionInfo, headerInfo, inputSamRdd)
      }

      val sampleNameBD = sc.broadcast(sampleName).value
      val headerInfoBD = sc.broadcast(headerInfo).value
      val refPartitionInfoBD = sc.broadcast(refPartitionInfo).value
      if(samRecordPartitionListRDD == null) {
        samRecordPartitionListRDD = samRecordPartitionRDD
          .map(bundle => {
            val partitionId = bundle._1
            val samRecordPartition = bundle._2
            val samMap = Map(sampleNameBD -> samRecordPartition)
            (partitionId, samMap)
          })
      } else {
        samRecordPartitionListRDD = samRecordPartitionListRDD.fullOuterJoin(samRecordPartitionRDD)
          .map(bundle => {
            val partitionId = bundle._1
            val samMapOption = bundle._2._1
            val partitionOption = bundle._2._2
            var samMap = samMapOption.getOrElse(Map[String, SamRecordPartition]())
            samMap += (sampleNameBD -> partitionOption.getOrElse(
              SamRecordPartition.empty(partitionId, refPartitionInfoBD.getPartitionRange(partitionId)._3, headerInfoBD)))
            (partitionId, samMap)
          })
      }
    })

    samRecordPartitionListRDD
  }

  protected def getBundlePartition(): RDD[BundlePartition] = {
    inputSamBundleMap.values.foreach(inputSamBundle =>
      if (!inputSamBundle.isSet) {
        throw new ResourceNotSetException(inputSamBundle.key)
      })

    // overlap length
    val overlapLength = PipelineConst.OVERLAP_LENGTH

    // reference
    val referencePartitionRDD = getReferencePartition(overlapLength)

    // rods
    val rodPartitionRDD = getRODPartitions(overlapLength)

    // input sam
    val samPartitionRDD = getSamPartitions(overlapLength)

    val refPartitionInfo = refPartitionInfoBundle.refPartitionInfo
    val refContigInfoBD = sc.broadcast(refPartitionInfo.getRefContigInfo).value

    referencePartitionRDD.leftOuterJoin(samPartitionRDD).join(rodPartitionRDD)
      .map(bundle => {
        val partitionId = bundle._1
        val referencePartition: FastaPartition = bundle._2._1._1
        val samPartitionMap: Map[String, SamRecordPartition] = bundle._2._1._2.getOrElse(Map[String, SamRecordPartition]())
        val rodPartitionMap: Map[String, VcfRecordPartition] = bundle._2._2
        new BundlePartition(partitionId, refContigInfoBD, referencePartition, samPartitionMap, rodPartitionMap)
      })
  }
}
