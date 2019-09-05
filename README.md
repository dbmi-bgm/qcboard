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
## Usage
```
python qcboard_v1.py -bam [BAM] -out [OUT TITLE]
```

### For BAM QC

#### Single individual

```
python qcboard_v1.py -bam [BAM] -out [OUT TITLE]
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

