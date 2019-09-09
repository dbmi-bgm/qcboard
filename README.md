# QCBoard


## Support
* GRCh37 (hg19)
* GRCh38 (hg38)

## Prerequisites
* samtools (>1.7)
* PICARD (current version)
* bedtools (current version)
* fastqc
* unzip

<!--## Installation
```
pip install qcboard
```
-->
## Usage of BAMQC
```
python qcboard_v1.py -bam [BAM] -out [OUT TITLE]
```


### Pre-executed commands
```
python bin/bamstat.py data/test1.bam > data/test1.bam.stat
fastqc ./data/test1.bam
java -Xms5G -Xmx5G -jar picard.jar CollectMultipleMetrics -VALIDATION_STRINGENCY LENIENT -INPUT data/test1.bam -OUTPUT data/test1.bam.cmm
python qcboard_v1.py -bam data/test1.bam -out data/test1
bedtools genomecov -ibam data/test1.bam -g data/PUB/b38.iff > data/test1.bam.hist
```

#### Bamstat.py

* input: ordered bam file
* output: read count table based on bam index (bai)
* dependency:
	* samtools (>1.1)
* expected running time : < 1 sec

```
python bin/bamstat.py [BAM] > [OUT.stat]

ex) samtools idxstats data/test1.bam > data/test1.bam.stat
```


<!--#### Samtools idxstat

* input: ordered bam file
* output: read count table based on bam index (bai)
* dependency:
	* samtools (>1.1)
* expected running time : < 1 sec

```
samtools idxstats [BAM] > [OUT.stat]

ex) samtools idxstats data/test1.bam > data/test1.bam.samtools.idxstats
```-->


#### FASTQC

* input: ordered bam file
* output: fastqc result file (html and metrics files)
	* `test1_fastqc.html`, `test1_fastqc.zip`
* dependency:
	* fastqc ([https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)) (>= 0.11.8)
	* JAVA
* expected running time for x30 : 10~20min

```
fastqc [BAM]

ex) fastqc ./data/test1.bam
```

#### PICARD CollectMultipleMetrics

* input: ordered bam file
* output: metrics files
	* `test1.bam.cmm.alignment_summary_metrics`
	* `test1.bam.cmm.base_distribution_by_cycle_metrics`
	* `test1.bam.cmm.base_distribution_by_cycle.pdf`
	* `test1.bam.cmm.insert_size_histogram.pdf`
	* `test1.bam.cmm.insert_size_metrics`
	* `test1.bam.cmm.quality_by_cycle_metrics`
	* `test1.bam.cmm.quality_by_cycle.pdf`
	* `test1.bam.cmm.quality_distribution_metrics`
	* `test1.bam.cmm.quality_distribution.pdf`
* dependency:
	* picard (>1.0) : [https://broadinstitute.github.io/picard/](https://broadinstitute.github.io/picard/)
* expected running time for x30 : 10~20min

```
java -Xms5G -Xmx5G -jar picard.jar \
     CollectMultipleMetrics  \
     -VALIDATION_STRINGENCY LENIENT \
     -INPUT [BAM] \
     -OUTPUT [BAM].cmm


ex) java -Xms5G -Xmx5G -jar picard.jar \
     CollectMultipleMetrics  \
     -VALIDATION_STRINGENCY LENIENT \
     -INPUT data/test1.bam \
     -OUTPUT data/test1.bam.cmm
```



#### BEDTOOLS
* input: ordered bam file
* output: coverage file (.hist)
* dependency: bedtools

```
bedtools genomecov -ibam [BAM] -g [IFF file] > [BAM].hist

ex) bedtools genomecov -ibam data/test1.bam -g data/PUB/b38.iff > data/test1.bam.hist
```


<!--

```
qcboard -bam [BAM] -out [OUT]
	-cmm [CMM]
	-fastqc [FASTQC]
```

-->
<!--### For VCF QC
```
qcboard -vcf [VCF] -out [OUT]
```
-->

