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

import java.io._

import net.java.truecommons.io.Loan._
import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD
import org.ncic.bioinfo.sparkseq.data.basic.{BasicSamRecord, FastqPairRecord, FastqRecord, VcfRecord}
import org.ncic.bioinfo.sparkseq.data.common.{SamHeaderInfo, VcfHeaderInfo}

/**
  * Author: wbc
  */
object NormalFileWriter extends FileWriter {

  private def getBufferedReader(filePath: String): BufferedWriter = {
    new BufferedWriter(new java.io.FileWriter(new File(filePath)))
  }

  override def writeFastq(records: Iterable[FastqRecord], filePath: String) = {
    loan(getBufferedReader(filePath)) to (writer => {
      records.foreach(record => {
        writer.write(record.toString())
        writer.newLine()
      })
    })
  }

  override def writeFastqRdd(sc: SparkContext, rdd: RDD[FastqRecord], filePath: String) = {
    //TODO 会导致master上的内存需求特别大。目前该方法只能用于调试。
    rdd.repartition(1).saveAsTextFile(filePath)
  }

  def writeFastqPair(records: Iterable[FastqPairRecord], filePath1: String, filePath2: String) = {
    // write file 1
    loan(getBufferedReader(filePath1)) to (writer => {
      records.foreach(record => {
        writer.write(record.descriptionLine)
        writer.newLine()
        writer.write(new String(record.sequence1))
        writer.newLine()
        writer.write(new String(record.quality1))
        writer.newLine()
      })
    })
    // write file 2
    loan(getBufferedReader(filePath2)) to (writer => {
      records.foreach(record => {
        writer.write(record.descriptionLine)
        writer.newLine()
        writer.write(new String(record.sequence2))
        writer.newLine()
        writer.write(new String(record.quality2))
        writer.newLine()
      })
    })
  }

  def writeFastqPairRdd(sc: SparkContext, rdd: RDD[FastqPairRecord],
                        filePath1: String, filePath2: String) = {
    //TODO 会导致master上的内存需求特别大。目前该方法只能用于调试。
    val records = rdd.collect().toList
    writeFastqPair(records, filePath1, filePath2)
  }

  override def writeSam(headerInfo: SamHeaderInfo, records: Iterable[BasicSamRecord], filePath: String) = {
    val sortedSAMRecords = records.toList.sortWith((record1, record2) => {
      if (record1.contigId != record2.contigId){
        record1.contigId < record2.contigId
      } else {
        record1.position < record2.position
      }
    })

    loan(getBufferedReader(filePath)) to (writer => {
      headerInfo.getHeaderLines.foreach(line => {
        writer.write(line)
        writer.newLine()
      })
      sortedSAMRecords.foreach(record => {
        writer.write(record.toString())
        writer.newLine()
      })
    })
  }

  override def writeSamRdd(sc: SparkContext, headerInfo: SamHeaderInfo,
                           rdd: RDD[BasicSamRecord], filePath: String) = {
    //TODO collect会导致master上的内存需求特别大。目前该方法只能用于调试。
    val records = rdd.collect().toList
    writeSam(headerInfo, records, filePath)
  }

  override def writeVcf(headerInfo: VcfHeaderInfo, records: Iterable[VcfRecord], filePath: String) = {
    val sortedVcfRecords = records.toList.sortWith((record1, record2) => {
      if (record1.contigId != record2.contigId){
        record1.contigId < record2.contigId
      } else {
        record1.position < record2.position
      }
    })

    loan(getBufferedReader(filePath)) to (writer => {
      headerInfo.getHeaderLines.foreach(line => {
        writer.write(line)
        writer.newLine()
      })
      sortedVcfRecords.foreach(record => {
        writer.write(record.toString())
        writer.newLine()
      })
    })
  }

  override def writeVcfRdd(sc: SparkContext, headerInfo: VcfHeaderInfo,
                           rdd: RDD[VcfRecord], filePath: String) = {
    val records = rdd.collect().toList
    writeVcf(headerInfo, records, filePath)
  }
}
