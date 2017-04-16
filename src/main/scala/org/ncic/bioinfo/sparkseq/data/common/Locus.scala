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

import org.ncic.bioinfo.sparkseq.utils.StringUtils

/**
  * Author: wbc
  */
object Locus {
  def apply(line: String, refContigInfo: RefContigInfo): Locus = {
    val iter = StringUtils.split(line, "\t")
    val contigName = iter.next()
    val contigId = refContigInfo.getId(contigName)
    val start = iter.next().toInt
    val stop = iter.next().toInt
    /*val pos1 = line.indexOf(':')
    val contigName = line.substring(0, pos1)
    val contigId = refContigInfo.getId(contigName)
    val pos2 = line.indexOf('-')
    val start = line.substring(pos1 + 1, pos2).toInt
    val stop = line.substring(pos2 + 1).toInt*/
    new Locus(contigId, contigName, start, stop)
  }
}

class Locus(val contigId: Int, val contigName: String, val start: Int, val stop: Int) extends Serializable {

}
