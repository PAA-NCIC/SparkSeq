package org.ncic.bioinfo.sparkseq.engine

import org.junit.runner.RunWith
import org.ncic.bioinfo.sparkseq.resource.Resource
import org.scalatest.FunSuite
import org.scalatest.junit.JUnitRunner

/**
  * Author: wbc
  */
@RunWith(classOf[JUnitRunner])
class TestEngine extends FunSuite {

  test("Test pipeline engine") {
    val dictFilePath = getClass().getResource("/littleFasta.dict").getFile()
    val referenceFilePath = getClass().getResource("/littleFasta.fasta").getFile()
    val pipeline = Pipeline("pipeName", null)
    val res1 = new MockResource("res 1")
    val res2 = new MockResource("res 2")
    res1.setFlag = true
    res2.setFlag = true
    pipeline.addResource(res1)
    pipeline.addResource(res2)

    val res3: Resource = new MockResource("res 3")
    val res4: Resource = new MockResource("res 4")
    val proc1 = new MockProcess("proc1", List(res1), List(res3))
    val proc2 = new MockProcess("proc2", List(res1, res2), List(res4))

    val res5 = new MockResource("res 5")
    val proc3 = new MockProcess("proc3", List(res1, res3), List(res5))

    val res6 = new MockResource("res 6")
    val proc4 = new MockProcess("proc4", List(res1, res3, res5), List(res6))

    val res7 = new MockResource("res 7")
    val proc5 = new MockProcess("proc5", List(res5, res6), List(res7))

    pipeline.addProcess(proc5)
    pipeline.addProcess(proc4)
    pipeline.addProcess(proc3)
    pipeline.addProcess(proc2)
    pipeline.addProcess(proc1)

    pipeline.run()

    assert(res7.isSet)
  }

  class MockProcess(name: String,
                    requiredResources:List[Resource],
                    generateResources:List[Resource]) extends AbstractProcess(name) {

    def getInputResourceList(): List[Resource] = requiredResources

    def getOutputResourceList(): List[Resource] = generateResources

    override def runProcess(): Unit = {
      print("Run " + name)
      outputResources.foreach(resource => resource.asInstanceOf[MockResource].setFlag = true)
    }

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