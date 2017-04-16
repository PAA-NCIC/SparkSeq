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

import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD
import org.ncic.bioinfo.sparkseq.data.basic.{BasicSamRecord, FastqPairRecord, FastqRecord, VcfRecord}
import org.ncic.bioinfo.sparkseq.data.common.{SamHeaderInfo, VcfHeaderInfo}

/**
  * Author: wbc
  */
trait FileWriter {

  def writeFastq(records: Iterable[FastqRecord], filePath: String)

  def writeFastqRdd(sc: SparkContext, rdd: RDD[FastqRecord], filePath: String)

  def writeFastqPair(records: Iterable[FastqPairRecord], filePath1: String, filePath2: String)

  def writeFastqPairRdd(sc: SparkContext, rdd: RDD[FastqPairRecord], filePath1: String, filePath2: String)

  def writeSam(headerInfo: SamHeaderInfo, records: Iterable[BasicSamRecord], filePath: String)

  def writeSamRdd(sc: SparkContext, headerInfo: SamHeaderInfo, rdd: RDD[BasicSamRecord], filePath: String)

  def writeVcf(headerInfo: VcfHeaderInfo, records: Iterable[VcfRecord], filePath: String)

  def writeVcfRdd(sc: SparkContext, headerInfo: VcfHeaderInfo, rdd: RDD[VcfRecord], filePath: String)
}
