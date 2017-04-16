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

/**
  * Author: wbc
  */
object FastqRecord extends Serializable {

  def apply(descriptionLine: String, sequence: Array[Byte], quality: Array[Byte], compressFlag: Boolean): FastqRecord = {
    new FastqRecord(compressFlag, descriptionLine, sequence, quality)
  }
}

class FastqRecord(val compressFlag: Boolean,
                  val descriptionLine: String,
                  val sequence: Array[Byte],
                  val quality: Array[Byte]
                 ) extends Serializable {

  override def toString(): String = {
    "%s\n%s\n+\n%s".format(descriptionLine, sequence, quality)
  }

}
