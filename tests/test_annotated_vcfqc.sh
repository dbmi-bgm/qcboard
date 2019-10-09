# python ../qcboard.py vcfstat -vcf data/Ashkenazim/HG002.hs37d5.60x.mnvindelmixed.vqsr.annot.vcf.gz
python ../qcboard.py vcfqc -vcf data/Ashkenazim/HG002.hs37d5.60x.mnvindelmixed.vqsr.annot.vcf.gz -out HG002.hs37d5.60x.annot
python ../qcboard.py vcfqc -vcf data/Ashkenazim/HG003.hs37d5.60x.mnvindelmixed.vqsr.annot.vcf.gz -out HG003.hs37d5.60x.annot
python ../qcboard.py vcfqc -vcf data/Ashkenazim/HG004.hs37d5.60x.mnvindelmixed.vqsr.annot.vcf.gz -out HG004.hs37d5.60x.annot
# python ../qcboard.py vcfqc -vcf data/Ashkenazim/HG002.hs37d5.60x.mnvindelmixed.vqsr.annot.vcf.gz data/Ashkenazim/HG003.hs37d5.60x.mnvindelmixed.vqsr.annot.vcf.gz data/Ashkenazim/HG004.hs37d5.60x.mnvindelmixed.vqsr.annot.vcf.gz -ped C,F,M -out Trio.hs37d5.60x.annot

