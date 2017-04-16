package org.ncic.bioinfo.sparkseq.engine

import org.apache.spark.rdd.RDD
import org.junit.runner.RunWith
import org.ncic.bioinfo.sparkseq.data.partition.BundlePartition
import org.ncic.bioinfo.sparkseq.resource.Resource
import org.scalatest.FunSuite
import org.scalatest.junit.JUnitRunner

/**
  * Author: wbc
  */
@RunWith(classOf[JUnitRunner])
class TestOptimizeEngine extends FunSuite {

  test("Test optimize pipeline engine") {
    val pipeline = Pipeline("pipeName", null)
    val res1 = new MockResource("res 1")
    val res2 = new MockResource("res 2")
    res1.setFlag = true
    res2.setFlag = true
    pipeline.addResource(res1)
    pipeline.addResource(res2)

    val res3: Resource = new MockResource("res 3")
    val res4: Resource = new MockResource("res 4")

    val normal1 = new Normal("normal1", List(res1), List(res3))
    val normal2 = new Normal("normal2", List(res1, res2), List(res4))

    val res5 = new MockResource("res 5")
    val res6 = new MockResource("res 6")
    val generator1 = new Generator("generator1", List(res4, res3), List(res5, res6))

    val res7 = new MockResource("res 7")
    val res8 = new MockResource("res 8")
    val generator2 = new Generator("generator2", List(res6, res5), List(res7, res8))

    val res9 = new MockResource("res 9")
    val res10 = new MockResource("res 10")
    val consumer1 = new Consumer("consumer1", List(res7, res8), List(res9, res10))

    val res11 = new MockResource("res 11")
    val generator3 = new Generator("generator3", List(res9, res10), List(res11))

    val res12 = new MockResource("res 12")
    val res13 = new MockResource("res 13")
    val consumer2 = new Consumer("consumer2", List(res11), List(res12, res13))

    pipeline.addProcess(generator2)
    pipeline.addProcess(consumer2)
    pipeline.addProcess(normal1)
    pipeline.addProcess(normal2)
    pipeline.addProcess(consumer1)
    pipeline.addProcess(generator3)
    pipeline.addProcess(generator1)

    pipeline.run()

    assert(res12.isSet)
    assert(res13.isSet)
  }

  test("Test optimize pipeline engine2") {
    val pipeline = Pipeline("pipeName", null)
    val res1 = new MockResource("res 1")
    val res2 = new MockResource("res 2")
    res1.setFlag = true
    res2.setFlag = true
    pipeline.addResource(res1)
    pipeline.addResource(res2)

    val res3: Resource = new MockResource("res 3")
    val res4: Resource = new MockResource("res 4")

    val normal1 = new Normal("normal1", List(res1), List(res3))
    val normal2 = new Normal("normal2", List(res2), List(res4))

    val res5 = new MockResource("res 5")
    val generator1 = new Generator("generator1", List(res3), List(res5))

    val res6 = new MockResource("res 6")
    val generator2 = new Generator("generator2", List(res4, res5), List(res6))

    val res7 = new MockResource("res 7")
    val consumer1 = new Consumer("consumer1", List(res6), List(res7))


    pipeline.addProcess(generator2)
    pipeline.addProcess(normal1)
    pipeline.addProcess(normal2)
    pipeline.addProcess(consumer1)
    pipeline.addProcess(generator1)

    pipeline.run()

    assert(res7.isSet)
  }

  class Generator(name: String,
                  requiredResources: List[Resource],
                  generateResources: List[Resource])
    extends PartitionOptimizedProcess(name) with PartitionGenerator {

    def getInputResourceList(): List[Resource] = requiredResources

    def getOutputResourceList(): List[Resource] = generateResources

    protected def getBundlePartition(): RDD[BundlePartition] = {
      null
    }

    def generateBundlePartition(bundle: RDD[BundlePartition]): RDD[BundlePartition] = {
      bundle
    }

    def consumeBundlePartition(bundle: RDD[BundlePartition]): Unit = {
      outputResources.foreach(resource => resource.asInstanceOf[MockResource].setFlag = true)
    }

    override def toString():String = name
  }

  class Consumer(name: String,
                 requiredResources: List[Resource],
                 generateResources: List[Resource])
    extends PartitionOptimizedProcess(name) with PartitionConsumer {

    def getInputResourceList(): List[Resource] = requiredResources

    def getOutputResourceList(): List[Resource] = generateResources

    protected def getBundlePartition(): RDD[BundlePartition] = {
      null
    }

    def consumeBundlePartition(bundle: RDD[BundlePartition]): Unit = {
      outputResources.foreach(resource => resource.asInstanceOf[MockResource].setFlag = true)
    }

    override def toString():String = name
  }

  class Normal(name: String,
               requiredResources: List[Resource],
               generateResources: List[Resource])
    extends AbstractProcess(name) {

    def getInputResourceList(): List[Resource] = requiredResources

    def getOutputResourceList(): List[Resource] = generateResources

    def runProcess(): Unit = {
      outputResources.foreach(resource => resource.asInstanceOf[MockResource].setFlag = true)
    }

    override def toString():String = name
  }

  class MockResource(resourceKey: String) extends Resource() {
    val key: String = resourceKey
    var setFlag: Boolean = false

    def isSet: Boolean = {
      setFlag
    }

    def getValue: Object = {
      null
    }

    def setValue(value: Object) = {}
  }

}