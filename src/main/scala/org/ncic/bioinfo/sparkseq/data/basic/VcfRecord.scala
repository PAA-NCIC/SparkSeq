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

import java.util.NoSuchElementException

import org.ncic.bioinfo.sparkseq.data.common.RefContigInfo
import org.ncic.bioinfo.sparkseq.utils.StringUtils

import scala.collection.mutable.ListBuffer

/**
  * Author: wbc
  */
object VcfRecord extends Serializable {

  def apply(vcfLine: String, refContigInfo: RefContigInfo): VcfRecord = {
    try {
      val strIter = StringUtils.split(vcfLine, "\t")
      val contigName = strIter.next()
      val contigId = refContigInfo.getId(contigName)
      val position = strIter.next().toInt
      val id = strIter.next()
      val ref = strIter.next()
      val alt = strIter.next()
      val qualToken: String = strIter.next()
      val qual = if (qualToken.equals(".")) 0 else qualToken.toDouble
      val filter = strIter.next()

      val rest = new StringBuilder()
      while (strIter.hasNext) {
        rest.append(strIter.next);
        if (strIter.hasNext) {
          rest.append('\t')
        }
      }

      /*// attrs
      val attributs = ListBuffer[(String, String)]()
      if (strIter.hasNext) {
        val attributeLine = strIter.next()
        val attrIter = StringUtils.split(attributeLine, ";")
        while (attrIter.hasNext) {
          val keyValueStr = attrIter.next()
          val iter = StringUtils.split(keyValueStr, "=")
          attributs += ((iter.next(), iter.next()))
        }
      }

      // extends attrs
      val extendsAttrs = ListBuffer[(String, String)]()
      if (strIter.hasNext) {
        val extendsAttrKeys = strIter.next()
        val extendsAttrVals = strIter.next()
        val eKeyIter = StringUtils.split(extendsAttrKeys, ":")
        val eValIter = StringUtils.split(extendsAttrVals, ":")
        while (eKeyIter.hasNext) {
          val key = eKeyIter.next()
          val value = eValIter.next()
          extendsAttrs += ((key, value))
        }
      }*/

      new VcfRecord(contigName, contigId, position, id, ref, alt, qual, filter, rest.toString())
    } catch {
      case ex: NoSuchElementException => throw new RuntimeException("Illegal vcf format")
    }
  }
}

class VcfRecord(val contigName: String,
                val contigId: Int,
                val position: Int, // 以1为base的contig内偏移
                val id: String,
                val ref: String,
                val alt: String,
                val qual: Double,
                val filter: String,
                val rest: String)
  extends Serializable {

  def getQualString():String = {
    if ((qual == 0)) return "."
    val tmp = "%.2f".format(qual)
    if (tmp.endsWith(".00")) tmp.substring(0, tmp.length - 3)
    //else if(tmp.endsWith("0")) tmp.substring(0, tmp.length - 1)
    else tmp
  }

  override def toString(): String = {

    var resTmp = "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s".format(
      contigName, position, id, ref, alt, getQualString(),
      filter, rest)

    /*if (attributes != null && attributes.size > 0) {
      val attriSb = new StringBuilder
      attributes.sortBy(_._1)
      val attrIter = attributes.iterator
      while (attrIter.hasNext) {
        val keyVal = attrIter.next()
        attriSb.append(keyVal._1)
        attriSb.append("=")
        attriSb.append(keyVal._2)
        if (attrIter.hasNext) {
          attriSb.append(";")
        }
      }
      resTmp += "\t" + attriSb.toString()
    }

    if (extendAttrs != null && extendAttrs.size > 0) {
      val sbKeySb = new StringBuilder
      val sbValSb = new StringBuilder
      val extIter = extendAttrs.iterator
      while (extIter.hasNext) {
        val keyVal = extIter.next()
        sbKeySb.append(keyVal._1)
        sbValSb.append(keyVal._2)
        if (extIter.hasNext) {
          sbKeySb.append(":")
          sbValSb.append(":")
        }
      }
      resTmp += "\t" + sbKeySb.toString() + "\t" + sbValSb.toString()
    }*/

    resTmp
  }
}
