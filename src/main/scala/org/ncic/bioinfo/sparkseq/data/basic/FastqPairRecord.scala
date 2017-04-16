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

import org.ncic.bioinfo.sparkseq.exceptions.PipelineException

/**
  * 为了进行pair-end的mapping
  * Author: wbc
  */
object FastqPairRecord extends Serializable {

  def apply(fastqRecord1: FastqRecord, fastqRecord2: FastqRecord): FastqPairRecord = {
    if (fastqRecord1.compressFlag != fastqRecord2.compressFlag) {
      throw new PipelineException("Can't zip a compressed fastq with a none-compressed fastq")
    }
    new FastqPairRecord(fastqRecord1.compressFlag, fastqRecord1.descriptionLine,
      fastqRecord1.sequence, fastqRecord2.sequence,
      fastqRecord1.quality, fastqRecord2.quality)
  }
}

class FastqPairRecord(val compressFlag: Boolean,
                      val descriptionLine: String,
                      val sequence1: Array[Byte],
                      val sequence2: Array[Byte],
                      val quality1: Array[Byte],
                      val quality2: Array[Byte]) extends Serializable {

}
