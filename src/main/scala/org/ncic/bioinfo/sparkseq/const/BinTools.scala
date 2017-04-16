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
package org.ncic.bioinfo.sparkseq.const

import java.io.FileInputStream
import java.util.Properties

import org.ncic.bioinfo.sparkseq.utils.FileUtils

/**
  * Author: wbc
  */
object BinTools {
  val binDirPath = FileUtils.join(FileUtils.getDirPathNoEndSeparator(
    FileUtils.getDirPathNoEndSeparator(
      FileUtils.getDirPathNoEndSeparator(
        FileUtils.getDirPathNoEndSeparator(
          FileUtils.getDirPathNoEndSeparator(
            FileUtils.getDirPathNoEndSeparator(
              FileUtils.getDirPathNoEndSeparator(
                FileUtils.getDirPathNoEndSeparator(
                  this.getClass().getResource("").getPath())))))))), "bin")

  val bwaPath = {
    val tmpPath = FileUtils.join(binDirPath, "bwa")
    if (tmpPath.startsWith("file:")) tmpPath.substring(5) else tmpPath
  }

  val bwaLibPath = {
    val tmpPath = FileUtils.join(binDirPath, "libbwajni.so")
    if (tmpPath.startsWith("file:")) tmpPath.substring(5) else tmpPath
  }

  val confPath = {
    val tmpPath = FileUtils.join(binDirPath, "config.properties")
    if (tmpPath.startsWith("file:")) tmpPath.substring(5) else tmpPath
  }

  val localTmpFilePath = {
    val properties = new Properties()
    properties.load(new FileInputStream(confPath))
    properties.getProperty("localWorkDir")
  }

  val publicTmpFilePath = {
    val properties = new Properties()
    properties.load(new FileInputStream(confPath))
    properties.getProperty("globalWorkDir")
  }

  val processOptimize = {
    val properties = new Properties()
    properties.load(new FileInputStream(confPath))
    properties.getProperty("processOptimize").toBoolean
  }

  val shuffleCompress = {
    val properties = new Properties()
    properties.load(new FileInputStream(confPath))
    properties.getProperty("shuffleCompress").toBoolean
  }

  val DEFAULT_PARTITION_LENGTH = {
    val properties = new Properties()
    properties.load(new FileInputStream(confPath))
    properties.getProperty("partitonLength").toInt
  }

  val splitPartitionThres = {
    val properties = new Properties()
    properties.load(new FileInputStream(confPath))
    properties.getProperty("splitPartitionThres").toInt
  }

  val bqsrGatherThreads = {
    val properties = new Properties()
    properties.load(new FileInputStream(confPath))
    properties.getProperty("bqsrGatherThreads").toInt
  }

  val repartitionCount = {
    val properties = new Properties()
    properties.load(new FileInputStream(confPath))
    properties.getProperty("activeRegionRepartitionCount").toInt
  }
}
