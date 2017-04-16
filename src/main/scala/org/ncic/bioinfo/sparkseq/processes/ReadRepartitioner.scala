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

import org.ncic.bioinfo.sparkseq.data.bundle.{RefPartitionInfoBundle, SAMBundle}
import org.ncic.bioinfo.sparkseq.data.common.RefPartitionInfo
import org.ncic.bioinfo.sparkseq.engine.AbstractProcess
import org.ncic.bioinfo.sparkseq.exceptions.{ResourceNotSetException, ResourceSetException}
import org.ncic.bioinfo.sparkseq.resource.Resource
import org.apache.spark.SparkContext._
import org.apache.spark.storage.StorageLevel
import org.ncic.bioinfo.sparkseq.const.BinTools

/**
  * @author wbc
  */
object ReadRepartitioner {
  def apply(name: String,
            inputSamBundle: SAMBundle,
            outputRefPartitionInfo: RefPartitionInfoBundle,
            refPartitionInfo: RefPartitionInfo): ReadRepartitioner =
    new ReadRepartitioner(name, inputSamBundle, outputRefPartitionInfo, refPartitionInfo)
}

class ReadRepartitioner(name: String,
                        inputSamBundle: SAMBundle,
                        outputRefPartitionInfo: RefPartitionInfoBundle,
                        refPartitionInfo: RefPartitionInfo) extends AbstractProcess(name) {

  override def getInputResourceList(): List[Resource] = List(inputSamBundle)

  override def getOutputResourceList(): List[Resource] = List(outputRefPartitionInfo)

  override def runProcess(): Unit = {
    if (inputSamBundle != null && !inputSamBundle.isSet) {
      throw new ResourceNotSetException(inputSamBundle.key)
    }
    if (outputRefPartitionInfo.isSet) {
      throw new ResourceSetException(outputRefPartitionInfo.key)
    }

    val samRDD = inputSamBundle.samRecordRDD
    samRDD.persist(StorageLevel.DISK_ONLY)
    val refPartitionInfoValue = sc.broadcast(refPartitionInfo).value
    val readStatistic = samRDD.map(record =>
      (refPartitionInfoValue.getPartitionId(record.contigId, record.position), 1))
      .reduceByKey(_ + _)
      .collect()

    val averageLength = readStatistic.map(_._2).sum / readStatistic.length
    val skipPartitions = readStatistic.filter(pair => pair._2 > BinTools.splitPartitionThres * averageLength)
      .map(_._1).toSet

    println("*****Skip partition count:" + skipPartitions.size)
    outputRefPartitionInfo.refPartitionInfo = RefPartitionInfo(refPartitionInfo, skipPartitions)
    outputRefPartitionInfo.setFlag = true
  }

}
