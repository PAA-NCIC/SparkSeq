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
import org.ncic.bioinfo.sparkseq.data.basic.BasicSamRecord
import org.ncic.bioinfo.sparkseq.data.common.SamHeaderInfo
import org.ncic.bioinfo.sparkseq.resource.AbstractResource

/**
  * Author: wbc
  */
object SAMBundle {

  /**
    * 未赋值的bundle
    *
    * @param key
    * @param samHeaderInfo
    * @return
    */
  def undefined(key: String,
                samHeaderInfo: SamHeaderInfo): SAMBundle = {
    val bundle = new SAMBundle(key, samHeaderInfo, null)
    bundle.setFlag = false
    bundle
  }

  /**
    * 从父亲sam bundle那里继承header信息，数据部分为空，是undefined的
    *
    * @param key resource key
    * @param parentSamBundle 父亲bundle
    * @return
    */
  def undefinedFromParent(key: String,
                parentSamBundle: SAMBundle): SAMBundle = {
    val bundle = new SAMBundle(key, parentSamBundle.samHeaderInfo, null)
    bundle.setFlag = false
    bundle
  }

  /**
    * 赋值的bundle
    *
    * @param key
    * @param samHeaderInfo
    * @param samRecordRDD
    * @return
    */
  def defined(key: String,
              samHeaderInfo: SamHeaderInfo,
              samRecordRDD: RDD[BasicSamRecord]): SAMBundle = {
    val bundle = new SAMBundle(key, samHeaderInfo, samRecordRDD)
    bundle.setFlag = true
    bundle
  }
}

class SAMBundle(key: String,
                val samHeaderInfo: SamHeaderInfo,
                var samRecordRDD: RDD[BasicSamRecord]) extends AbstractResource(key) {

}
