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
package org.ncic.bioinfo.sparkseq.resource

import org.ncic.bioinfo.sparkseq.exceptions.ResourceException

import scala.collection.mutable

/**
  * Author: wbc
  */
object ResourcePool {
  def apply(): ResourcePool = new ResourcePool()
}

class ResourcePool {

  // 维护两套hash，以hashMap为主，通过resource自带的key索引，hashSet用于快速值索引
  private val resourceMap: mutable.HashMap[String, Resource] = mutable.HashMap()
  private val resourceSet: mutable.HashSet[Resource] = mutable.HashSet()

  /**
    * 将资源添加进资源池
    * 如果有重复的key，则抛出异常
    *
    * @param resource 待添加的资源
    * @throws ResourceException
    */
  def addResource(resource: Resource): Unit = {
    if (resourceMap.contains(resource.key)) {
      throw new ResourceException("Same resource key: " + resource.key)
    }
    resourceMap.put(resource.key, resource)
    resourceSet.add(resource)
  }

  /**
    * 将资源添加进资源池，如果已经有重复的key，则替换原有的resource
    *
    * @param resource
    */
  def replaceResource(resource: Resource): Unit = {
    if (resourceMap.contains(resource.key)) {
      val resourceOld = resourceMap.get(resource.key).get
      //删除原有的resource，即使他们可能只是key相同
      resourceSet.remove(resourceOld)
    }
    resourceMap.put(resource.key, resource)
    resourceSet.add(resource)
  }

  def containsResourceKey(key: String): Boolean = {
    resourceMap.contains(key)
  }

  def containsResource(resource: Resource): Boolean = {
    resourceSet.contains(resource)
  }

  def getResourceByKey(key: String): Option[Resource] = {
    resourceMap.get(key)
  }

  def copy(): ResourcePool = {
    val newPool = ResourcePool()
    resourceSet.foreach(resource => newPool.addResource(resource))
    newPool
  }

}
