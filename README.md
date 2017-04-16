SparkSeq
====

# Introduction

SparkSeq is a programming framework for big genome data analysis process based on [Apache Spark](http://spark.incubator.apache.org/).

In the latest version, WGS pipline is implemented based on this framwork.

# Get Started

## Installing Spark
A Spark release must be on your system before running SparkSeq.
Our work default to Spark 2.1.0 and Hadoop 2.7.0. We also test the program on Spark 1.6.2 and Hadoop 2.6.4.
The latest release of Spark can be download from [Spark website](http://spark.apache.org/downloads.html). More information about installing Spark refers to [Spark Installing Document](https://github.com/apache/spark).

## Build SparkSeq
SparkSeq is built using Apache Maven. To build SparkSeq and its example program, run the following command in shell:
```
$ git clone https://github.com/PAA-NCIC/SparkSeq.git
$ cd SparkSeq
$ mvn clean package
...
[INFO] ------------------------------------------------------------------------
[INFO] BUILD SUCCESS
[INFO] ------------------------------------------------------------------------
[INFO] Total time: 01:00 min
[INFO] Finished at: 2017-04-03T19:43:34+08:00
[INFO] Final Memory: 140M/2083M
[INFO] ------------------------------------------------------------------------
```

## Run a example of WGS
A pipline to call raw variants in VCF file format from raw reads in FASTQ file format is implemented as an example program in GPF.

Before running a WGS pipline, a human b37 reference and its index file must be provided in a storage which can be seen by each node in Spark cluster, as the BWA mem task need to load the index into memory when run .
The known indel/snp vcf files and input FASTQ files must be provided in HDFS or NFS system. 
We provide an example running script named "runSparkWGS.sh" in directory "SparkSeq/bin" like following:
```
spark_master=spark://master:7077
driver_memory=30G
executor_memory=30G
total_executor_cores=1024

spark-submit --class org.ncic.bioinfo.sparkseq.WGSPipeline \
 --master ${spark_master} \
 --driver-memory ${driver_memory} \
 --executor-memory ${executor_memory} \
 --total-executor-cores ${total_executor_cores} \
 /PATH/TO/SparkSeq/target/spark-seq-0.9.0-jar-with-dependencies.jar \
 -ref /PATH/TO/human_g1k_v37.fasta \
 -dict /PATH/TO/human_g1k_v37.dict \
 -fq1 /PATH/TO/DATA/1.fastq \
 -fq2 /PATH/TO/DATA/2.fastq \
 -output /PATH/TO/OUTPUT/result.vcf \
 -1000gindel /PATH/TO/1000G_phase1.indels.b37.vcf \
 -millsindel /PATH/TO/Mills_and_1000G_gold_standard.indels.b37.vcf \
 -dbsnp /PATH/TO/dbsnp_138.b37.vcf
```

You can run the script by run the following command in shell:
```
sh runSparkWGS.sh
```

The WGS pipeline support all arguments defined by [Apache Spark](http://spark.incubator.apache.org/), and also defines a series of arguments for WGS pipline.
````
Arguments for process defination.
 -fq1				: Path to input fastq file 1.
 -fq2				: Path to input fastq file 2.
 -ref				: Path to b37 reference.
 -dict				: Path to dict file of b37 reference.
 -output			: Path to write the output VCF file.
 -1000gindel			: Path to VCF file: 1000G known indels
 -millsindel			: Path to VCF file: mills and 1000G indels
 -dbsnp				: Path to VCF file: dbsnp
````

# License
SparkSeq is released under a [GNU General Public License](https://github.com/PAA-NCIC/SparkSeq/master/LICENSE).
