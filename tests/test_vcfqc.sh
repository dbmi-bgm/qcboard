# python ../qcboard.py vcfqc -vcf data/test.vcf -out test_vcf

### vcfstat 
# python ../qcboard.py vcfstat -vcf data/Ashkenazim/HG002.hs37d5.60x.mnvindelmixed.vqsr.vcf.gz
# python ../qcboard.py vcfstat -vcf data/Ashkenazim/HG003.hs37d5.60x.mnvindelmixed.vqsr.vcf.gz
# python ../qcboard.py vcfstat -vcf data/Ashkenazim/HG004.hs37d5.60x.mnvindelmixed.vqsr.vcf.gz

python ../qcboard.py vcfqc -vcf data/Ashkenazim/HG002.hs37d5.60x.mnvindelmixed.vqsr.vcf.gz -out HG002.hs37d5.60x
python ../qcboard.py vcfqc -vcf data/Ashkenazim/HG003.hs37d5.60x.mnvindelmixed.vqsr.vcf.gz -out HG003.hs37d5.60x
python ../qcboard.py vcfqc -vcf data/Ashkenazim/HG004.hs37d5.60x.mnvindelmixed.vqsr.vcf.gz -out HG004.hs37d5.60x
