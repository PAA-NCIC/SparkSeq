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
package org.ncic.bioinfo.sparkseq.processes.mapping

import org.apache.spark.rdd.RDD
import org.ncic.bioinfo.sparkseq.const.BinTools
import org.ncic.bioinfo.sparkseq.data.bundle.{FASTQBundle, FASTQPairBundle, SAMBundle}
import org.ncic.bioinfo.sparkseq.data.common.{ReadGroupInfo, RefContigInfo}
import org.ncic.bioinfo.sparkseq.engine.AbstractProcess
import org.ncic.bioinfo.sparkseq.exceptions.{ResourceNotSetException, ResourceSetException}
import org.ncic.bioinfo.sparkseq.fileio.{NormalFileLoader, NormalFileWriter}
import org.ncic.bioinfo.sparkseq.resource.Resource
import org.ncic.bioinfo.sparkseq.utils.FileUtils

import scala.collection.mutable
import scala.collection.mutable.ListBuffer
import scala.sys.process._

/**
  * Author: wbc
  */
object BwaMappingProcess {
  def singleEnd(name: String,
                partitionCount: Int,
                refContigInfo: RefContigInfo,
                referencePath: String,
                fastqBundle: FASTQBundle,
                outputSamBundle: SAMBundle): BwaMappingProcess = {
    val process = new BwaMappingProcess(name,
      partitionCount, refContigInfo, referencePath, outputSamBundle,
      fastqBundle, null)
    process
  }

  def pairEnd(name: String,
              partitionCount: Int,
              refContigInfo: RefContigInfo,
              referencePath: String,
              fastqPairBundle: FASTQPairBundle,
              outputSamBundle: SAMBundle): BwaMappingProcess = {
    val process = new BwaMappingProcess(
      name, partitionCount, refContigInfo, referencePath, outputSamBundle,
      null, fastqPairBundle)
    process
  }
}

/**
  * 进行bwa mapping操作，输入fastq或者fastq pair，返回samRecord的RDD
  * 在mapping过程中会占用tmp目录中的空间作为中间文件，当前版本中需要手动清除
  * 在mapping最后会去除掉unmap的reads
  *
  * @param name
  * @param outputSamBundle
  * @param referencePath
  */
class BwaMappingProcess(name: String,
                        partitionCount: Int,
                        refContigInfo: RefContigInfo,
                        referencePath: String,
                        outputSamBundle: SAMBundle,
                        fastqBundle: FASTQBundle,
                        fastqPairBundle: FASTQPairBundle) extends AbstractProcess(name) {

  final val THREAD_COUNT = 12

  val arguments = mutable.HashMap[String, String]()

  def addArguments(key: String, value: String) = {
    arguments.put(key, value)
  }

  def addArguments(key: String) = {
    arguments.put(key, "")
  }

  override def getInputResourceList(): List[Resource] = {
    val resourceList = ListBuffer[Resource]()
    if (fastqBundle != null) {
      resourceList += fastqBundle
    }
    if (fastqPairBundle != null) {
      resourceList += fastqPairBundle
    }
    resourceList.toList
  }

  override def getOutputResourceList(): List[Resource] = List(outputSamBundle)

  private def buildArgString(arguments: mutable.HashMap[String, String]): String = {
    arguments.put("-t", THREAD_COUNT.toString)
    val sBuilder = new StringBuilder
    arguments.foreach(entry => {
      sBuilder.append(" %s %s".format(entry._1, entry._2))
    })
    sBuilder.toString()
  }

  private def transeReadGroupInfoIntoArg(readGroupInfo: ReadGroupInfo): String = {
    "@RG\\tID:%s\\tSM:%s\\tPL:%s\\tLB:%s\\tPU:%s".format(
      readGroupInfo.id, readGroupInfo.sample, readGroupInfo.platform,
      readGroupInfo.lib, readGroupInfo.platformUnit
    )
  }

  private def runSingle(): Unit = {
    // 输入
    val fastqRecordRdd = fastqBundle.fastqRecordRDD
    // 将readGroupInfo的信息加入arguments
    // 如果之前的sam header中不包含read group，则创建一个
    val readGroupInfo = {
      if (outputSamBundle.samHeaderInfo.getReadGroupInfos().isEmpty) {
        val readGroup = ReadGroupInfo.apply("rg1", "sample1")
        outputSamBundle.samHeaderInfo.addReadGroupInfo(readGroup)
        readGroup
      } else {
        outputSamBundle.samHeaderInfo.getReadGroupInfos().head
      }
    }
    addArguments("-R", transeReadGroupInfoIntoArg(readGroupInfo))

    // 将fatsq record划分成多个partition
    val fastqRepartitionRdd = {
      if (partitionCount > 0) {
        // bwa进程数 = partition数 / bwa线程数的向上取整
        val bwaPartitionCount = (partitionCount + THREAD_COUNT - 1) / THREAD_COUNT
        fastqRecordRdd.repartition(bwaPartitionCount)
      } else {
        fastqRecordRdd
      }
    }

    val bwaPath = BinTools.bwaPath
    val tmpSamDirPath = FileUtils.join(BinTools.publicTmpFilePath, "tmpSam")
    // 广播变量
    val bwaPathBD = sc.broadcast(bwaPath).value
    val tmpSamDirPathBD = sc.broadcast(tmpSamDirPath).value
    val localTmpPathBD = sc.broadcast(BinTools.localTmpFilePath).value
    val argStrBD = sc.broadcast(buildArgString(arguments)).value
    val refContigInfoBD = sc.broadcast(refContigInfo).value
    val referencePathBD = sc.broadcast(referencePath).value
    val sampleNameValue = sc.broadcast(outputSamBundle.key).value

    //运行
    if (!FileUtils.isExists(tmpSamDirPath)) {
      "mkdir -p %s".format(tmpSamDirPath) !;
    }

    val samRddLinesRDD: RDD[String] = fastqRepartitionRdd.mapPartitionsWithIndex((partitionId, recordIterator) => {
      val tmpFastqFilePath = FileUtils.join(
        localTmpPathBD, "%s_%d.fastq".format(sampleNameValue, partitionId))
      NormalFileWriter.writeFastq(recordIterator.toIterable, tmpFastqFilePath)
      // tmpSam存放在公共目录下
      val tmpSamFilePath = FileUtils.join(
        tmpSamDirPathBD, "%s_%d.sam".format(sampleNameValue, partitionId))
      // trans to commands
      val command = "%s mem %s %s %s"
        .format(bwaPathBD, argStrBD, referencePathBD, tmpFastqFilePath)
      command #> new java.io.File(tmpSamFilePath) !;
      "rm %s".format(tmpFastqFilePath) !;
      Seq("Empty").iterator
    })

    // 提前collect生成一个action，虽然没有意义，但是可以让其他的job不会先于bwa执行
    // 这种方法使得bwa的进程可以均匀地划分到各个节点上。
    samRddLinesRDD.collect()

    val samRdd = NormalFileLoader.loadSamToRdd(sc, FileUtils.join(tmpSamDirPath, outputSamBundle.key + "*"), refContigInfo)

    // set result
    outputSamBundle.samRecordRDD = samRdd
    outputSamBundle.setFlag = true
  }

  private def runPair(): Unit = {
    // 输入
    val fastqPairRecordRdd = fastqPairBundle.fastqPairRecordRDD

    // 将readGroupInfo的信息加入arguments
    // 如果之前的sam header中不包含read group，则创建一个
    val readGroupInfo = {
      if (outputSamBundle.samHeaderInfo.getReadGroupInfos().isEmpty) {
        val readGroup = ReadGroupInfo.apply("rg1", "sample1")
        outputSamBundle.samHeaderInfo.addReadGroupInfo(readGroup)
        readGroup
      } else {
        outputSamBundle.samHeaderInfo.getReadGroupInfos().head
      }
    }
    addArguments("-R", transeReadGroupInfoIntoArg(readGroupInfo))

    // 将fatsq record划分成多个partition
    val fastqPairRepartitionRdd = {
      if (partitionCount > 0) {
        // bwa进程数 = partition数 / bwa线程数的向上取整
        val bwaPartitionCount = (partitionCount + THREAD_COUNT - 1) / THREAD_COUNT
        fastqPairRecordRdd.repartition(bwaPartitionCount)
      } else {
        fastqPairRecordRdd
      }
    }

    val bwaPath = BinTools.bwaPath
    val tmpSamDirPath = FileUtils.join(BinTools.publicTmpFilePath, "tmpSam")
    // 广播变量
    val bwaPathBD = sc.broadcast(bwaPath).value
    val tmpSamDirPathBD = sc.broadcast(tmpSamDirPath).value
    val localTmpPathBD = sc.broadcast(BinTools.localTmpFilePath).value
    val argStrBD = sc.broadcast(buildArgString(arguments)).value
    val refContigInfoBD = sc.broadcast(refContigInfo).value
    val referencePathBD = sc.broadcast(referencePath).value
    val sampleNameValue = sc.broadcast(outputSamBundle.key).value

    //运行
    if (!FileUtils.isExists(tmpSamDirPath)) {
      "mkdir -p %s".format(tmpSamDirPath) !;
    }

    val samRddLinesRDD: RDD[String] = fastqPairRepartitionRdd.mapPartitionsWithIndex((partitionId, recordIterator) => {
      val tmpFastqFilePath1 = FileUtils.join(
        localTmpPathBD, "%s_%d_1.fastq".format(sampleNameValue, partitionId))
      val tmpFastqFilePath2 = FileUtils.join(
        localTmpPathBD, "%s_%d_2.fastq".format(sampleNameValue, partitionId))
      NormalFileWriter.writeFastqPair(recordIterator.toIterable, tmpFastqFilePath1, tmpFastqFilePath2)
      // tmpSam存放在公共目录下
      val tmpSamFilePath = FileUtils.join(
        tmpSamDirPathBD, "%s_%d.sam".format(sampleNameValue, partitionId))
      // trans to commands
      val command = "%s mem %s %s %s %s"
        .format(bwaPathBD, argStrBD, referencePathBD,
          tmpFastqFilePath1, tmpFastqFilePath2)
      command #> new java.io.File(tmpSamFilePath) !;
      "rm %s %s".format(tmpFastqFilePath1, tmpFastqFilePath2) !;
      Seq("Empty").iterator
    })

    // 提前collect生成一个action，虽然没有意义，但是可以让其他的job不会先于bwa执行
    // 这种方法使得bwa的进程可以均匀地划分到各个节点上。
    samRddLinesRDD.collect()

    val samRdd = NormalFileLoader.loadSamToRdd(sc, FileUtils.join(tmpSamDirPath, outputSamBundle.key + "*"), refContigInfo)

    // set result
    outputSamBundle.samRecordRDD = samRdd
    outputSamBundle.setFlag = true
  }

  override def runProcess(): Unit = {

    if (fastqBundle != null && !fastqBundle.isSet) {
      throw new ResourceNotSetException(fastqBundle.key)
    }
    if (fastqPairBundle != null && !fastqPairBundle.isSet) {
      throw new ResourceNotSetException(fastqPairBundle.key)
    }
    if (outputSamBundle.isSet) {
      throw new ResourceSetException(outputSamBundle.key)
    }

    if (fastqBundle != null) {
      runSingle()
    }
    if (fastqPairBundle != null) {
      runPair()
    }
  }
}
