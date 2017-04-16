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
package org.ncic.bioinfo.sparkseq

import java.util

import org.apache.spark.{SparkConf, SparkContext}
import org.kohsuke.args4j.{Argument, CmdLineParser, Option}
import org.ncic.bioinfo.sparkseq.algorithms.data.vcf.RODNames
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLoc
import org.ncic.bioinfo.sparkseq.algorithms.walker.SerializableActiveRegionMapData
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.ActiveRegion
import org.ncic.bioinfo.sparkseq.data.basic.{BasicSamRecord, FastqPairRecord, FastqRecord, VcfRecord}
import org.ncic.bioinfo.sparkseq.data.bundle.{FASTQPairBundle, RefPartitionInfoBundle, SAMBundle, VCFBundle}
import org.ncic.bioinfo.sparkseq.data.common._
import org.ncic.bioinfo.sparkseq.engine.Pipeline
import org.ncic.bioinfo.sparkseq.fileio.{NormalFileLoader, NormalFileWriter}
import org.ncic.bioinfo.sparkseq.processes.ReadRepartitioner
import org.ncic.bioinfo.sparkseq.processes.cleaning.{BaseRecalibrationProcess, IndelRealignProcess, MarkDuplicateProcess, PartitionMarkDuplicateProcess}
import org.ncic.bioinfo.sparkseq.processes.mapping.{BwaMappingProcess, JNIBwaMemProcess}
import org.ncic.bioinfo.sparkseq.processes.variantcalling.HaplotypeCallerProcess
import org.ncic.bioinfo.sparkseq.utils.StringUtils

/**
  * Author: wbc
  */
object WGSPipeline {
  @Option(required = true, name = "-ref", usage = "reference")
  val reference: String = StringUtils.EMPTY

  @Option(required = true, name = "-dict", usage = "dictionary of reference")
  val dictFilePath: String = StringUtils.EMPTY

  @Option(required = true, name = "-fq1", usage = "fastq r1 for sample")
  val fastq1: String = StringUtils.EMPTY

  @Option(required = false, name = "-fq2", usage = "fastq r2 for sample")
  val fastq2: String = StringUtils.EMPTY

  @Option(required = true, name = "-output", usage = "result path")
  val resultVcfPath: String = StringUtils.EMPTY

  @Option(required = false, name = "-1000gindel", usage = "1000G known indels")
  val indels1000G: String = StringUtils.EMPTY

  @Option(required = false, name = "-millsindel", usage = "mills and 1000G indels")
  val indelsMills: String = StringUtils.EMPTY

  @Option(required = false, name = "-dbsnp", usage = "known vcfs")
  val dbsnp: String = StringUtils.EMPTY

  @Argument
  val arguments: util.ArrayList[String] = new util.ArrayList[String]()

  def main(args: Array[String]): Unit = {
    // parse argument
    val parser: CmdLineParser = new CmdLineParser(this)
    parser.setUsageWidth(300)
    // parse the arguments.
    val argList = new util.ArrayList[String]()
    args.foreach(arg => argList.add(arg))
    parser.parseArgument(argList)

    val conf = new SparkConf()
      .setAppName("SparkSeq")
      .set("spark.driver.maxResultSize","24g")
      .set("spark.serializer","org.apache.spark.serializer.KryoSerializer")
      .registerKryoClasses(Array(classOf[BasicSamRecord], classOf[FastqPairRecord], classOf[FastqRecord], classOf[VcfRecord],
        classOf[SerializableActiveRegionMapData], classOf[ActiveRegion], classOf[GenomeLoc]))

    if (conf.getOption("spark.master").isEmpty) {
      conf.setMaster("local[%d]".format(Runtime.getRuntime.availableProcessors()))
    }

    val sc = new SparkContext(conf)
    val pipelineName = "myPipeline"
    val pipeline = Pipeline(pipelineName, sc)

    val refContigInfo = RefContigInfo(dictFilePath)
    val refPartitionInfo = RefPartitionInfo(refContigInfo)

    val samHeaderInfo = SamHeaderInfo.unsortedHeader(refContigInfo, null)
    val fastqPairRdd = NormalFileLoader.loadFastqPairToRdd(sc, fastq1, fastq2)
    val fastqPairBundle = FASTQPairBundle.defined("fastqPair", fastqPairRdd)
    val alignedSAMBundle = SAMBundle.undefined("alignedSam", samHeaderInfo)

    val mappingProcess = JNIBwaMemProcess.pairEnd("MyBwaMapping",
      refContigInfo, reference, fastqPairBundle, alignedSAMBundle)
    pipeline.addProcess(mappingProcess)

    val repartitionInfoBundle = RefPartitionInfoBundle.undefined("repartition")
    val repartitionProcess = ReadRepartitioner("MyRepartition", alignedSAMBundle, repartitionInfoBundle, refPartitionInfo)
    pipeline.addProcess(repartitionProcess)

    val dedupedSamBundle = SAMBundle.undefinedFromParent("dedupedSam", alignedSAMBundle)
    val markDuplicateProcess = PartitionMarkDuplicateProcess("MyMarkDuplicate", reference, repartitionInfoBundle, List(alignedSAMBundle), List(dedupedSamBundle))
    pipeline.addProcess(markDuplicateProcess)

    val realignedSamBundle = SAMBundle.undefinedFromParent("realignedSam", dedupedSamBundle)
    val rodMapRealign = Map(RODNames.INDELS_1000G_PHASE1 -> indels1000G, RODNames.INDELS_MILLS_AND_1000G -> indelsMills)
    val indelRealignerProcess = IndelRealignProcess("MyIndelRealigner", reference, rodMapRealign, repartitionInfoBundle, List(dedupedSamBundle), List(realignedSamBundle))
    pipeline.addProcess(indelRealignerProcess)

    val recaledSamBundle = SAMBundle.undefinedFromParent("recaledSam", realignedSamBundle)
    val rodMapRecal = Map(RODNames.INDELS_1000G_PHASE1 -> indels1000G)
    val baseRecalibrationProcess = BaseRecalibrationProcess("MyBecalibrationProcess", reference, rodMapRecal, repartitionInfoBundle, List(realignedSamBundle), List(recaledSamBundle))
    pipeline.addProcess(baseRecalibrationProcess)

    val vcfHeader = VcfHeaderInfo.newHeader(refContigInfo, List())
    val vcfBundle = VCFBundle.undefined("ResultVCF", vcfHeader)
    val rodMapCall = Map(RODNames.DBSNP -> dbsnp)
    val samHeaderInfoBD = recaledSamBundle.samHeaderInfo
    val vcfHeaderInfoBD = NormalFileLoader.loadVcfHeader(dbsnp)
    val haplotypeCallerProcess = HaplotypeCallerProcess("MyHaplotypeCaller", reference, rodMapCall, repartitionInfoBundle, List(recaledSamBundle), vcfBundle, true,
      samHeaderInfoBD, vcfHeaderInfoBD)
    pipeline.addProcess(haplotypeCallerProcess)

    pipeline.run()

    val vcfRecords = vcfBundle.vcfRecordRDD.collect()

    //将结果RDD中的数据写入文件
    NormalFileWriter.writeVcf(VcfHeaderInfo.newHeader(refContigInfo, null), vcfRecords, resultVcfPath)
  }
}
