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
package org.ncic.bioinfo.sparkseq.debug

import java.io.{BufferedWriter, File, FileWriter}

import net.java.truecommons.io.Loan._
import org.ncic.bioinfo.sparkseq.const.BinTools
import org.ncic.bioinfo.sparkseq.data.basic.{BasicSamRecord, VcfRecord}
import org.ncic.bioinfo.sparkseq.data.common.{ReadGroupInfo, RefContigInfo, SamHeaderInfo, VcfHeaderInfo}
import org.ncic.bioinfo.sparkseq.data.partition.{BundlePartition, FastaPartition, SamRecordPartition, VcfRecordPartition}
import org.ncic.bioinfo.sparkseq.utils.FileUtils

import scala.collection.mutable.ListBuffer
import scala.io.Source
import scala.sys.process._
import collection.JavaConversions._

/**
  * Author: wbc
  */
object Dumper {
  val defaultString = "/tmp"

  def dumpSamRecordPartition(samRecordPartition: SamRecordPartition, filePath: String): Unit = {
    loan(new BufferedWriter(new FileWriter(new File(filePath)))) to (writer => {
      writer.write(samRecordPartition.partitionId.toString)
      writer.newLine()
      writer.write(samRecordPartition.contigId.toString)
      writer.newLine()
      writer.write(samRecordPartition.samHeaderInfo.readGroupInfos.head.id)
      writer.newLine()
      writer.write(samRecordPartition.samHeaderInfo.readGroupInfos.head.sample)
      writer.newLine()
      samRecordPartition.records.foreach(record => {
        writer.write(record.toString())
        writer.newLine()
      })
    })
  }

  def dedumpSamRecordPartition(filePath: String, refContigInfo: RefContigInfo): SamRecordPartition = {

    val lines = Source.fromFile(filePath).getLines()
    val partitionId = lines.next().toInt
    val contigId = lines.next().toInt
    val rgId = lines.next()
    val sample = lines.next()

    val samHeaderInfo = SamHeaderInfo.sortedHeader(refContigInfo, List())
    samHeaderInfo.addReadGroupInfo(ReadGroupInfo(rgId, sample))

    val records = ListBuffer[BasicSamRecord]()
    while (lines.hasNext) {
      val line = lines.next
      if (line.length > 0) {
        records += BasicSamRecord(line, refContigInfo, false)
      }
    }

    new SamRecordPartition(partitionId, contigId, records, samHeaderInfo)
  }

  def dumpRefPartition(fastaPartition: FastaPartition, filePath: String): Unit = {
    loan(new BufferedWriter(new FileWriter(new File(filePath)))) to (writer => {
      writer.write(fastaPartition.partitionId.toString)
      writer.newLine()
      writer.write(fastaPartition.contigId.toString)
      writer.newLine()
      writer.write(fastaPartition.contigName)
      writer.newLine()
      writer.write(fastaPartition.content)
      writer.newLine()
      writer.write(fastaPartition.originStart.toString)
      writer.newLine()
      writer.write(fastaPartition.originEnd.toString)
      writer.newLine()
      writer.write(fastaPartition.overlappedStart.toString)
      writer.newLine()
      writer.write(fastaPartition.overlappedEnd.toString)
      writer.newLine()
      writer.write(fastaPartition.safeOverlappedStart.toString)
      writer.newLine()
      writer.write(fastaPartition.safeOverlappedEnd.toString)
      writer.newLine()
    })
  }

  def deDumpRefPartition(filePath: String): FastaPartition = {
    val lines = Source.fromFile(filePath).getLines()
    val partitionId = lines.next.toInt
    val contigId = lines.next.toInt
    val contigName = lines.next
    val content = lines.next
    val originStart = lines.next.toInt
    val originEnd = lines.next.toInt
    val overlappedStart = lines.next.toInt
    val overlappedEnd = lines.next.toInt
    val safeOverlapppedStart = lines.next.toInt
    val safeOverlappedEnd = lines.next.toInt

    new FastaPartition(partitionId, contigId, contigName, content,
      safeOverlapppedStart, safeOverlappedEnd, overlappedStart, overlappedEnd, originStart, originEnd)
  }

  def dumpRODPartitions(vcfRecordPartition: VcfRecordPartition, filePath: String): Unit = {
    loan(new BufferedWriter(new FileWriter(new File(filePath)))) to (writer => {
      writer.write(vcfRecordPartition.partitionId.toString)
      writer.newLine()
      writer.write(vcfRecordPartition.key)
      writer.newLine()
      writer.write(vcfRecordPartition.contigId.toString)
      writer.newLine()
      vcfRecordPartition.vcfHeader.getHeaderLines.foreach(line => {
        writer.write(line)
        writer.newLine()
      })
      vcfRecordPartition.records.foreach(record => {
        writer.write(record.toString())
        writer.newLine()
      })
    })
  }

  def deDumpRODPartitions(filePath: String, refContigInfo: RefContigInfo): VcfRecordPartition = {
    val lines = Source.fromFile(filePath).getLines()
    val partitionId = lines.next.toInt
    val key = lines.next
    val contigId = lines.next.toInt
    val headerLines = ListBuffer[String]()
    val records = ListBuffer[VcfRecord]()

    while (lines.hasNext) {
      val line = lines.next()
      if (line.length > 0) {
        if (line.startsWith("#")) {
          headerLines += line
        } else {
          records += VcfRecord(line, refContigInfo)
        }
      }
    }

    new VcfRecordPartition(partitionId, key, contigId, VcfHeaderInfo.newHeader(headerLines.toList), records)
  }

  def dumpBundle(bundlePartition: BundlePartition, baseDir: String): Unit = {
    //"mkdir -p %s".format(baseDir) !;
    val refDumpPath = FileUtils.join(baseDir, bundlePartition.partitionId + "_ref")
    dumpRefPartition(bundlePartition.fastaPartition, refDumpPath)
    bundlePartition.samRecordPartitionMap.foreach(pair => {
      val samDumpPath = FileUtils.join(baseDir, bundlePartition.partitionId + "_sam_" + pair._1)
      dumpSamRecordPartition(pair._2, samDumpPath)
    })
    bundlePartition.rodPartitionMap.foreach(pair => {
      val rodDumpPath = FileUtils.join(baseDir, bundlePartition.partitionId + "_rod_" + pair._1)
      dumpRODPartitions(pair._2, rodDumpPath)
    })
  }

  def dumpBQSRTableMap(tableMap: Map[String, java.util.List[String]], baseDir: String = defaultString): Unit = {
    tableMap.foreach(pair => {
      val key = pair._1
      val table = pair._2
      val filePath = FileUtils.join(baseDir, "bqsrTable_" + key)
      loan(new BufferedWriter(new FileWriter(new File(filePath)))) to (writer => {
        table.foreach(line => {
          writer.write(line)
          writer.newLine()
        })
      })
    })
  }
}
