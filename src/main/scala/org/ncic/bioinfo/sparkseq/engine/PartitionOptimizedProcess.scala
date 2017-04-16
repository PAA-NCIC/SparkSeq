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
package org.ncic.bioinfo.sparkseq.engine

import org.apache.spark.rdd.RDD
import org.ncic.bioinfo.sparkseq.data.partition.BundlePartition
import org.ncic.bioinfo.sparkseq.exceptions.PipelineException

/**
  * Author: wbc
  */
abstract class PartitionOptimizedProcess(name: String) extends AbstractProcess(name) {


  def runProcess(): Unit = {

    var inputBundle = if (this.isInstanceOf[PartitionConsumer]
      && this.asInstanceOf[PartitionConsumer].setAsChainConsumer) {
      this.asInstanceOf[PartitionConsumer].parentProcess.resultBundleRDD
    } else {
      getBundlePartition()
    }

    if (this.isInstanceOf[PartitionGenerator]) {
      //如果是generator，则根据是否设置为generator决定要把partition存起来还是消费掉
      val generator = this.asInstanceOf[PartitionGenerator]
      if (generator.setAsChainGenerator) {
        generator.resultBundleRDD = generator.generateBundlePartition(inputBundle)
      } else {
        generator.consumeBundlePartition(inputBundle)
      }
    } else if (this.isInstanceOf[PartitionConsumer]) {
      //如果是纯的consumer，就直接消费掉就好了
      val consumer = this.asInstanceOf[PartitionConsumer]
      consumer.consumeBundlePartition(inputBundle)
    }

  }

  protected def getBundlePartition(): RDD[BundlePartition]

  override def run(): Unit = {
    if (resourcePool == null) {
      throw new PipelineException("Resource pool is not set yet <process: " + name + ">")
    }

    // 运行process
    runProcess()

    done = true
  }
}
