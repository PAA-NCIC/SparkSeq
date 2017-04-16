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
package org.ncic.bioinfo.sparkseq.fileio

import java.io.{BufferedReader, File}

import org.apache.hadoop.io.Text
import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD
import org.ncic.bioinfo.sparkseq.compress.{BaseCompressTools, QualityCompressTools}
import org.ncic.bioinfo.sparkseq.const.{BinTools, SamRecordConst}
import org.ncic.bioinfo.sparkseq.data.basic.{BasicSamRecord, FastqPairRecord, FastqRecord, VcfRecord}
import org.ncic.bioinfo.sparkseq.data.common.{Locus, RefContigInfo, VcfHeaderInfo}
import org.ncic.bioinfo.sparkseq.fileio.format.SingleFastqInputFormat
import org.ncic.bioinfo.sparkseq.utils.StringUtils

import scala.collection.mutable.ListBuffer
import scala.io.Source

/**
  * Author: wbc
  */
object NormalFileLoader extends FileLoader {

  def transFilePath(filePath: String): String = {
    if (filePath.startsWith("hdfs://")) {
      filePath
    } else if (filePath.startsWith("file://"))
      filePath
    else
      "file://" + filePath
  }

  private def getReadNameFromDescptionLine(descriptionLine: String): String = {
    val nameEndIdx = descriptionLine.indexOf(' ')
    if (nameEndIdx < 0) {
      descriptionLine
    } else {
      val baseName = descriptionLine.substring(0, nameEndIdx)
      //兼容老版本，后面可能带有/1或/2
      if (baseName.endsWith("/1") || baseName.endsWith("/2")) {
        baseName.substring(0, baseName.length - 2)
      } else {
        baseName
      }
    }
  }

  override def loadFastq(filePath: String): List[FastqRecord] = {
    val lines = Source.fromFile(new File(filePath)).getLines()
    val fastqRecords = ListBuffer[FastqRecord]()
    var idx = 0
    var descriptionLine: String = null
    var readName: String = null
    var sequenceLine: String = null
    var qualityLine: String = null
    while (lines.hasNext) {
      if (idx % 4 == 0) {
        descriptionLine = lines.next()
        readName = getReadNameFromDescptionLine(descriptionLine)
      } else if (idx % 4 == 1) {
        sequenceLine = lines.next()
      } else if (idx % 4 == 2) {
        lines.next()
      } else {
        qualityLine = lines.next()
        fastqRecords += FastqRecord(readName, sequenceLine.getBytes, qualityLine.getBytes, false)
      }
      idx += 1
    }
    fastqRecords.toList
  }

  override def loadFastqToRdd(sc: SparkContext, filePath: String): RDD[FastqRecord] = {
    val records = sc.newAPIHadoopFile(
      transFilePath(filePath),
      classOf[SingleFastqInputFormat],
      classOf[Void],
      classOf[Text]
    )

    // 根据压缩设置选择是否在读入之后压缩
    val compressFlagValue = sc.broadcast(BinTools.shuffleCompress).value
    // convert records
    records.map(record => {
      val strIter = StringUtils.split(record._2.toString, "\n")
      val descriptionLine = strIter.next()
      val readName = getReadNameFromDescptionLine(descriptionLine)
      var sequenceBytes = strIter.next().getBytes()
      strIter.next() // abandon second description line
      var qualityBytes = strIter.next().getBytes
      if (compressFlagValue) {
        sequenceBytes = BaseCompressTools.compressBase(sequenceBytes, qualityBytes)
        qualityBytes = QualityCompressTools.compressQual(qualityBytes)
      }
      FastqRecord(readName, sequenceBytes, qualityBytes, compressFlagValue)
    })
  }

  override def loadFastqPair(filePath1: String, filePath2: String): List[FastqPairRecord] = {
    val records1 = loadFastq(filePath1)
    val records2 = loadFastq(filePath2)
    records1.zip(records2)
      .map(pair => FastqPairRecord(pair._1, pair._2))
      .toList
  }

  override def loadFastqPairToRdd(sc: SparkContext,
                                  filePath1: String, filePath2: String): RDD[FastqPairRecord] = {
    val rdd1 = loadFastqToRdd(sc, filePath1)
    val rdd2 = loadFastqToRdd(sc, filePath2)

    val allRdd = rdd1 ++ rdd2
    allRdd
      .groupBy(record => record.descriptionLine) // description line中只保留了read name，因此是按照read name group
      .map(bundle => {
      val iter = bundle._2.iterator
      val record1 = iter.next()
      val record2 = iter.next()
      FastqPairRecord(record1, record2)
    })
  }

  override def loadSam(filePath: String, refContigInfo: RefContigInfo): List[BasicSamRecord] = {
    Source.fromFile(new File(filePath)).getLines()
      .filter(_.charAt(0) != '@')
      .map(line => BasicSamRecord(line, refContigInfo, false))
      .filter(record => record.contigId != SamRecordConst.FAKE_CONTIG_ID)
      .toList
  }

  override def loadSamToRdd(sc: SparkContext, filePath: String,
                            refContigInfo: RefContigInfo): RDD[BasicSamRecord] = {
    val refContigInfoValue = sc.broadcast(refContigInfo).value

    // 根据压缩设置选择是否在读入之后压缩
    val compressFlagValue = sc.broadcast(BinTools.shuffleCompress).value
    sc.textFile(transFilePath(filePath))
      .filter(_.charAt(0) != '@')
      .map(line => BasicSamRecord(line, refContigInfoValue, compressFlagValue))
      .filter(record => record.contigId != SamRecordConst.FAKE_CONTIG_ID)
  }

  override def loadVcf(filePath: String, refContigInfo: RefContigInfo): List[VcfRecord] = {
    Source.fromFile(new File(filePath)).getLines()
      .filter(_.charAt(0) != '#')
      .map(line => VcfRecord(line, refContigInfo))
      .toList
  }

  override def loadVcfToRdd(sc: SparkContext, filePath: String,
                            refContigInfo: RefContigInfo): RDD[VcfRecord] = {
    sc.textFile(transFilePath(filePath))
      .filter(_.charAt(0) != '#')
      .map(line => VcfRecord(line, refContigInfo))
  }

  override def loadVcfHeader(filePath: String): VcfHeaderInfo = {
    val headerLines = ListBuffer[String]()

    if (filePath.startsWith("hdfs://")) {
      val reader = HDFSReader(filePath)
      var line = reader.readLine()
      while (line != null && line.startsWith("##")) {
        headerLines += line
        line = reader.readLine()
      }
    } else {
      val lines = Source.fromFile(filePath).getLines()
      while (lines.hasNext) {
        val line = lines.next()
        if (line.startsWith("##")) {
          headerLines += line
        } else {
          return VcfHeaderInfo.newHeader(headerLines.toList)
        }
      }
    }
    VcfHeaderInfo.newHeader(headerLines.toList)
  }

  def loadIntervals(filePath: String, refContigInfo: RefContigInfo): List[Locus] = {
    val lines = Source.fromFile(filePath).getLines()
    val locus = ListBuffer[Locus]()
    while (lines.hasNext) {
      val line = lines.next()
      if (!line.startsWith("@")) {
        locus += Locus(line, refContigInfo)
      }
    }
    locus.toList
  }
}
