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

import scala.collection.mutable.ListBuffer

/**
  * Author: wbc
  */

object SamHeaderInfo extends Serializable {

  def sortedHeader(refContigInfo: RefContigInfo, initExtLines: List[String]): SamHeaderInfo = {
    val attributeLines = List[String](getHDHeaderline("1.4"))
    var extLines = List[String]()
    if (initExtLines != null) {
      extLines = initExtLines
    }
    new SamHeaderInfo(true, attributeLines, refContigInfo, extLines)
  }

  def unsortedHeader(refContigInfo: RefContigInfo, initExtLines: List[String]): SamHeaderInfo = {
    val attributeLines = List[String]()
    var extLines = List[String]()
    if (initExtLines != null) {
      extLines = initExtLines
    }
    new SamHeaderInfo(false, attributeLines, refContigInfo, extLines)
  }

  private def getHDHeaderline(vn: String): String = "@HD\tVN:%s\tSO:coordinate".format(vn)

}

class SamHeaderInfo(var sorted: Boolean,
                    attributeLines: List[String],
                    val refContigInfo: RefContigInfo,
                    var extLines: List[String]) extends Serializable {
  var readGroupInfos = List[ReadGroupInfo]()

  private def getSQHeaderline(contigName: String, contigLength: Int): String =
    "@SQ\tSN:%s\tLN:%d".format(contigName, contigLength)

  private def getSqLines(refContigInfo: RefContigInfo): List[String] = {
    refContigInfo.getContigIds.map(contigId => {
      val name = refContigInfo.getName(contigId)
      val length = refContigInfo.getLength(contigId)
      getSQHeaderline(name, length)
    })
  }

  private def getReadGroupLine(readGroupInfo: ReadGroupInfo): String = {
    "@RG\tID:%s\tSM:%s\tPL:%s\tLB:%s\tPU:%s".format(
      readGroupInfo.id,
      readGroupInfo.sample,
      readGroupInfo.platform,
      readGroupInfo.lib,
      readGroupInfo.platformUnit
    )
  }

  def addReadGroupInfo(readGroupInfo: ReadGroupInfo): Unit = {
    readGroupInfos ::= readGroupInfo
  }

  def getHeaderLines: List[String] = {
    val buffer = ListBuffer[String]()
    buffer ++= attributeLines
    buffer ++= getSqLines(refContigInfo)
    // read group行
    readGroupInfos.foreach(info => buffer += getReadGroupLine(info))
    // 检查extern line是否都以'@'开头，如果不是就加上
    extLines.foreach(line => {
      if (line.charAt(0) != '@') {
        buffer += "@" + line
      } else {
        buffer += line
      }
    })

    buffer.toList
  }

  // 已经是public了仍然提供get方法，是为了在java中调用
  def getRefContigInfo() = refContigInfo

  def getReadGroupInfos() = readGroupInfos
}
