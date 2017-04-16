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
