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
import org.ncic.bioinfo.sparkseq.data.basic.VcfRecord
import org.ncic.bioinfo.sparkseq.resource.AbstractResource

/**
  * Author: wbc
  */
object RODBundle {

  /**
    * 未赋值的bundle
    *
    * @param key
    * @param rodName
    * @return
    */
  def undefined(key: String, rodName: String): RODBundle = {
    val bundle = new RODBundle(key, rodName, null);
    bundle.setFlag = false
    bundle
  }

  /**
    * 赋值的bundle
    *
    * @param key
    * @param rodName
    * @param vcfRcordRDD
    * @return
    */
  def defined(key: String, rodName: String, vcfRcordRDD: RDD[VcfRecord]): RODBundle = {
    val bundle = new RODBundle(key, rodName, vcfRcordRDD);
    bundle.setFlag = true
    bundle
  }

}

class RODBundle(key: String,
                val rodName: String,
                var vcfRecordRDD: RDD[VcfRecord]) extends AbstractResource(key) {

}
