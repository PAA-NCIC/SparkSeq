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
package org.ncic.bioinfo.sparkseq.data.bundle

import org.apache.spark.rdd.RDD
import org.ncic.bioinfo.sparkseq.data.basic.{FastqPairRecord, FastqRecord}
import org.ncic.bioinfo.sparkseq.data.common.RefContigInfo
import org.ncic.bioinfo.sparkseq.data.partition.FastaPartition
import org.ncic.bioinfo.sparkseq.resource.AbstractResource

/**
  * Author: wbc
  */
object FASTQBundle {

  /**
    * 未赋值的bundle
    *
    * @param key
    * @return
    */
  def undefined(key: String): FASTQBundle = {
    val bundle = new FASTQBundle(key, null)
    bundle.setFlag = false
    bundle
  }

  /**
    * 已赋值的bundle
    *
    * @param key
    * @param fastqRecordRDD
    * @return
    */
  def defined(key: String,
              fastqRecordRDD: RDD[FastqRecord]): FASTQBundle = {
    val bundle = new FASTQBundle(key, fastqRecordRDD)
    bundle.setFlag = true
    bundle
  }
}

class FASTQBundle(key: String, var fastqRecordRDD: RDD[FastqRecord])
  extends AbstractResource(key) {

}
