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
package org.ncic.bioinfo.sparkseq.processes.variantcalling

import org.apache.spark.rdd.RDD
import org.ncic.bioinfo.sparkseq.data.basic.VcfRecord
import org.ncic.bioinfo.sparkseq.data.bundle.{RefPartitionInfoBundle, SAMBundle, VCFBundle}
import org.ncic.bioinfo.sparkseq.data.common.RefPartitionInfo
import org.ncic.bioinfo.sparkseq.data.partition.BundlePartition
import org.ncic.bioinfo.sparkseq.engine.PartitionConsumer
import org.ncic.bioinfo.sparkseq.exceptions.ResourceNotSetException
import org.ncic.bioinfo.sparkseq.processes.PartitionProcess
import org.ncic.bioinfo.sparkseq.resource.Resource

import scala.collection.mutable.ListBuffer

/**
  * Author: wbc
  */
abstract class VariantCallingProcess(name: String,
                                     referencePath: String,
                                     rodMap: Map[String, String],
                                     refPartitionInfoBundle: RefPartitionInfoBundle,
                                     inputSamBundleMap: Map[String, SAMBundle],
                                     outputVCFBundle: VCFBundle)
  extends PartitionProcess(name, referencePath, rodMap, refPartitionInfoBundle, inputSamBundleMap)
    with PartitionConsumer {

  override def getInputResourceList(): List[Resource] = {
    val inputList = ListBuffer[Resource]()
    inputList ++= inputSamBundleMap.values
    inputList.append(refPartitionInfoBundle)
    inputList.toList
  }

  override def getOutputResourceList(): List[Resource] = List(outputVCFBundle)

  override def consumeBundlePartition(bundleRDD: RDD[BundlePartition]): Unit = {
    val vcfRDDs = getVCFRecords(bundleRDD)

    // 写入结果
    outputVCFBundle.vcfRecordRDD = vcfRDDs
    outputVCFBundle.setFlag = true
  }

  def getVCFRecords(bundleRDD: RDD[BundlePartition]): RDD[VcfRecord]
}
