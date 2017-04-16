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
package org.ncic.bioinfo.sparkseq.data.common

import java.io.File

import org.ncic.bioinfo.sparkseq.const.SamRecordConst
import org.ncic.bioinfo.sparkseq.exceptions.IllegalInputException
import org.ncic.bioinfo.sparkseq.utils.StringUtils

import scala.collection.mutable
import scala.io.Source
import collection.JavaConversions._

/**
  * Author: wbc
  */
object RefContigInfo {
  def apply(dictFilePath: String): RefContigInfo = {
    val info = new RefContigInfo()

    // unmap的contig位置是*
    info.name2IdMap += (SamRecordConst.FAKE_CONTIG_STR -> SamRecordConst.FAKE_CONTIG_ID)
    info.id2NameMap += (SamRecordConst.FAKE_CONTIG_ID -> SamRecordConst.FAKE_CONTIG_STR)
    info.id2LengthMap += (SamRecordConst.FAKE_CONTIG_ID -> 0)

    // 获取@SQ开头的line，描述了所有的contig
    val lines = Source.fromFile(new File(dictFilePath))
      .getLines()
      .filter(line => line.startsWith("@SQ")).toList

    if (lines.size == 0) {
      throw new IllegalInputException("非法的Dict文件：@SQ行数为0")
    }

    // 遍历lines，
    var idx = 0
    for (line <- lines) {
      val splitIter = StringUtils.split(line, "\t")
      val sqMark = splitIter.next()
      val contigName = splitIter.next().substring(3)
      val contigLength = splitIter.next().substring(3).toInt

      info.name2IdMap += (contigName -> idx)
      info.id2NameMap += (idx -> contigName)
      info.id2LengthMap += (idx -> contigLength)

      idx += 1
    }

    info.refType = getReferenceVersion(info.name2IdMap.keys.head)

    info
  }

  private def getReferenceVersion(contigName: String): String = {
    if (contigName.toLowerCase.startsWith("chr")) {
      "hg19"
    } else {
      "b37"
    }
  }
}

/**
  * id是从0开始的
  */
class RefContigInfo extends Serializable {

  private var refType: String = null

  private val id2NameMap = mutable.HashMap[Int, String]()

  private val name2IdMap = mutable.HashMap[String, Int]()

  private val id2LengthMap = mutable.HashMap[Int, Int]()

  def getName(id: Int): String = id2NameMap.getOrElse(id, SamRecordConst.FAKE_CONTIG_STR)

  def getId(name: String): Int = name2IdMap.getOrElse(name, SamRecordConst.FAKE_CONTIG_ID)

  def getLength(id: Int): Int = id2LengthMap.getOrElse(id, 0)

  def getLength(name: String): Int = getLength(getId(name))

  def getContigIds: List[Int] = Range(0, id2NameMap.size - 1).toList //减一是为了去掉fake_id

  def getContigIdsInteger: List[java.lang.Integer] = getContigIds.asInstanceOf[List[java.lang.Integer]]

  def getRefType: String = refType

}
