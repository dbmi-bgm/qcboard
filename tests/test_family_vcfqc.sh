

python ../qcboard.py vcfstat -vcf data/Ashkenazim/HG002.hs37d5.60x.mnvindelmixed.vqsr.vcf.gz \
	data/Ashkenazim/HG003.hs37d5.60x.mnvindelmixed.vqsr.vcf.gz \
	data/Ashkenazim/HG004.hs37d5.60x.mnvindelmixed.vqsr.vcf.gz \
	-ped C,F,M -out HG002.hs37d5.60x.trio
python ../qcboard.py vcfqc -vcf data/Ashkenazim/HG002.hs37d5.60x.mnvindelmixed.vqsr.vcf.gz \
	data/Ashkenazim/HG003.hs37d5.60x.mnvindelmixed.vqsr.vcf.gz \
	data/Ashkenazim/HG004.hs37d5.60x.mnvindelmixed.vqsr.vcf.gz \
	-ped C,F,M -out HG002.hs37d5.60x.trio

