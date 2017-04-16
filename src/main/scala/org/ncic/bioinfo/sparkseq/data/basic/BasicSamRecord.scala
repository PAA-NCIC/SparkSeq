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
package org.ncic.bioinfo.sparkseq.data.basic

import java.io.Serializable

import org.ncic.bioinfo.sparkseq.compress.{BaseCompressTools, QualityCompressTools}
import org.ncic.bioinfo.sparkseq.const.SamRecordConst
import org.ncic.bioinfo.sparkseq.data.common.RefContigInfo
import org.ncic.bioinfo.sparkseq.utils.StringUtils

import collection.JavaConversions._

/**
  * Author: wbc
  */
object BasicSamRecord extends Serializable {
  def apply(samLine: String, refContigInfo: RefContigInfo):BasicSamRecord = {
    apply(samLine, refContigInfo, false)
  }

  def apply(samLine: String, refContigInfo: RefContigInfo, compressFlag: Boolean): BasicSamRecord = {
    val splits = org.apache.commons.lang3.StringUtils.split(samLine, '\t')
    val readName = splits(0)
    val flag = splits(1).toInt
    val contigName = splits(2)
    val contigId = refContigInfo.getId(contigName)
    val position = splits(3).toInt
    val mapQ = splits(4).toInt
    val cigar = splits(5)
    val rawMateContigName = splits(6)
    var mateContigName = ""
    var mateContigId = SamRecordConst.FAKE_CONTIG_ID
    var matePosition = 0
    var infferdSize = 0
    if (!rawMateContigName.equals("*")) {
      if (rawMateContigName.equals("=")) {
        mateContigName = contigName
        mateContigId = contigId
      }
      else {
        mateContigName = rawMateContigName
        mateContigId = refContigInfo.getId(mateContigName)
      }
      matePosition = splits(7).toInt
      infferdSize = splits(8).toInt
    }
    var sequence = splits(9).getBytes
    var quality = splits(10).getBytes
    if (compressFlag) {
      sequence = BaseCompressTools.compressBase(sequence, quality)
      quality = QualityCompressTools.compressQual(quality)
    }
    var attributeList = splits.takeRight(splits.length - 11).toList

    new BasicSamRecord(compressFlag, readName, flag, contigId, contigName, position, mapQ,
      cigar, mateContigId, mateContigName, matePosition, infferdSize, sequence, quality, attributeList)
  }

  def apply(compressFlag: Boolean, readName: String, flag: Int, contigId: Int, contigName: String, position: Int, mapQ: Int,
            cigar: String, mateContigId: Int, mateContigName: String, matePosition: Int, infferdSize: Int,
            sequence: Array[Byte], quality: Array[Byte], attributeList: java.util.List[String]): BasicSamRecord = {

    new BasicSamRecord(compressFlag, readName, flag, contigId, contigName, position, mapQ,
      cigar, mateContigId, mateContigName, matePosition, infferdSize,
      sequence, quality, attributeList.toList)
  }

  def apply(compressFlag: Boolean, readName: String, flag: Int, contigId: Int, contigName: String, position: Int, mapQ: Int,
            cigar: String, mateContigId: Int, mateContigName: String, matePosition: Int, infferdSize: Int,
            sequence: Array[Byte], quality: Array[Byte], attributeList: List[String]): BasicSamRecord = {

    new BasicSamRecord(compressFlag, readName, flag, contigId, contigName, position, mapQ,
      cigar, mateContigId, mateContigName, matePosition, infferdSize,
      sequence, quality, attributeList)
  }
}

class BasicSamRecord(val compressFlag: Boolean,
                     val readName: String,
                     val flag: Int,
                     val contigId: Int,
                     val contigName: String,
                     val position: Int, // 以1为base的偏移
                     val mapQ: Int,
                     val cigar: String,
                     val mateContigId: Int,
                     val mateContigName: String,
                     val matePosition: Int,
                     val infferdSize: Int,
                     val sequence: Array[Byte],
                     val quality: Array[Byte],
                     val attributeList: List[String]) extends Serializable {

  override def toString(): String = {
    val displayMateName = if (mateContigId == contigId) "=" else mateContigName
    val qualityBytes = if(compressFlag) QualityCompressTools.deCompressQual(quality) else quality
    val sequenceBytes = if(compressFlag) BaseCompressTools.decompressBase(sequence, qualityBytes) else sequence
    "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s".format(
      readName, flag, contigName, position, mapQ, cigar, displayMateName, matePosition, infferdSize,
      new String(sequenceBytes), new String(qualityBytes), StringUtils.join(attributeList, '\t'))
  }

  def compress: BasicSamRecord = {
    if(compressFlag) {
      this
    } else {
      val compressSequence = BaseCompressTools.compressBase(sequence, quality)
      val compressQuality = QualityCompressTools.compressQual(quality)
      BasicSamRecord(true, readName, flag, contigId, contigName, position, mapQ,
        cigar, mateContigId, mateContigName, matePosition, infferdSize,
        compressSequence, compressQuality, attributeList)
    }
  }

}
