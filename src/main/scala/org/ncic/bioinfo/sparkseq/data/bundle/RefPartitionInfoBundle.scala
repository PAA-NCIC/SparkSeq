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

import org.ncic.bioinfo.sparkseq.data.common.RefPartitionInfo
import org.ncic.bioinfo.sparkseq.resource.AbstractResource

/**
  * @author wbc
  */
object RefPartitionInfoBundle {
  /**
    * 未实现内容
    *
    * @param key
    * @return
    */
  def undefined(key: String): RefPartitionInfoBundle = {
    val bundle = new RefPartitionInfoBundle(key, null)
    bundle.setFlag = false
    bundle
  }

  /**
    * 已实现内容
    *
    * @param key
    * @param refPartitionInfo
    * @return
    */
  def defined(key: String, refPartitionInfo: RefPartitionInfo): RefPartitionInfoBundle = {
    val bundle = new RefPartitionInfoBundle(key, refPartitionInfo)
    bundle.setFlag = true
    bundle
  }
}

class RefPartitionInfoBundle(key: String,
                             var refPartitionInfo: RefPartitionInfo) extends AbstractResource(key) {
}
