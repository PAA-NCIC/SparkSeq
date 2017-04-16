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
package org.ncic.bioinfo.sparkseq.utils

import org.ncic.bioinfo.sparkseq.data.common.Flags


/**
  * Author: wbc
  */
object FlagUtils {

  val MULTI_SEGMENT = 1
  val EACH_SEGMENT_ALIGNED = 2
  val SEGMENT_UNMAPPED = 4
  val NEXT_SEGMENT_UNMAPPED = 8
  val SEQUENCE_REVERSED = 16
  val NEXT_SEQUENCE_REVERSED = 32
  val FIRST_SEGMENT_IN_MAPPED = 64
  val LAST_SEGMENT_IN_MAPPED = 128
  val SECONDARY_ALIGNMENT = 256
  val NO_PASSING_FILTERS = 512
  val PCR_DUPLICATE = 1024
  val SUPPLEMENTARY_ALIGNMENT = 2048

  def initialFlag = new Flags
}