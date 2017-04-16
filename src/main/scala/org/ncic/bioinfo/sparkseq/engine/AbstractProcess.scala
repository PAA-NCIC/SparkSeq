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

import org.apache.spark.SparkContext
import org.ncic.bioinfo.sparkseq.exceptions.PipelineException
import org.ncic.bioinfo.sparkseq.resource.{Resource, ResourcePool}

/**
  * Author: wbc
  */
abstract class AbstractProcess() extends Process {

  // Process接口中resourcePool和sc的实现
  protected var pipeline: Pipeline = null
  protected var resourcePool: ResourcePool = null
  protected var sc: SparkContext = null
  var name: String = null

  def this(name: String) {
    this()
    this.name = name
    this.inputResources ++= getInputResourceList()
    this.outputResources ++= getOutputResourceList()
  }

  /**
    * 负责将output的资源set掉。
    */
  def runProcess(): Unit

  def getInputResourceList(): List[Resource]

  def getOutputResourceList(): List[Resource]

  override def run(): Unit = {
    if (resourcePool == null) {
      throw new PipelineException("Resource pool is not set yet <process: " + name + ">")
    }

    // 检查存在未set的输出resource
    if (outputResources == null || !outputResources.exists(resource => !resource.isSet)) {
      return
    }

    // 检查所有的输入resource是否已经存在
    if (inputResources != null && inputResources.exists(
      resource => !(resourcePool.containsResource(resource) && resource.isSet))) {
      throw new PipelineException("Missing resource for <process" + name + ">")
    }

    // 运行process
    runProcess()

    // 将输出的resource加入资源池
    if (outputResources != null) {
      outputResources.foreach(resource => {
        if (!resource.isSet) {
          throw new PipelineException("Output resource is not set in <process" + name + ">")
        }
        resourcePool.addResource(resource)
      })
    }

    done = true
  }
}
