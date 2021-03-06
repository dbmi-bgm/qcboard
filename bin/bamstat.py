#!/usr/bin/env python
# -*- coding: utf-8 -*-
#### bamstat.py
#### made by Daniel Minseok Kwon
#### 2017-10-29 16:07:31
#########################
import sys
import os

samtools = "samtools"

EXOM_SEQ_COVERED_REGION_NimbleGen_SeqCap=64190747
EXOM_SEQ_COVERED_REGION = EXOM_SEQ_COVERED_REGION_NimbleGen_SeqCap

SEQVER_CHR1LEN={249250621:'GRCh37',248956422:'GRCh38'}

MAPPABLE_LEN = {}
MAPPABLE_LEN['GRCh37'] = {'1':223574127,'2':237180652,'3':194676623,'4':187355871,'5':176782338,'6':167110905,'7':154688534,'8':142361290,'9':117088589,'10':130479461,'11':130841171,'12':130158001,'13':95519122,'14':88220591,'15':81022884,'16':78305341,'17':77413192,'18':74629856,'19':55591264,'20':59388246,'21':34967844,'22':34676104,'X':149110827,'Y':18340719,'MT':16569}
MAPPABLE_LEN['GRCh38'] = {'1':228278438,'2':239247912,'3':197472479,'4':188935388,'5':177502130,'6':169171305,'7':157169456,'8':143388316,'9':118976756,'10':131970072,'11':133841004,'12':132519294,'13':97763523,'14':88369084,'15':82225345,'16':79944795,'17':79417227,'18':78823041,'19':55570863,'20':62882875,'21':35754804,'22':36657544,'X':150912695,'Y':18326499,'M':16569}

def comma(value):
    return "{:,}".format(value)

def run_cmd(scmd, flag=False):
    if flag:
        print (scmd)
    rst = os.popen(scmd)
    rst_cont = rst.read()
    return rst_cont


def bamstat(bam):
    read_len = 0
    cmd_cont = run_cmd(samtools+ " view "+bam+" | head -1000")
    i = 0
    for line in cmd_cont.strip().split('\n'):
        if len(line) > 0 and line[0] != '@':
            arr= line.split('\t')
            if len(arr[9]) > read_len:
                read_len = len(arr[9])
            i += 1
            if i > 10:
                break

    cont = []
    cont.append("CHROM")
    cont.append("LEN")
    cont.append("MAAPED")
    cont.append("UNMAPPED")
    cont.append("TOTAL")
    cont.append("MAPPED_RATIO")
    cont.append("COVERAGE")
    cont.append("ADJUSTED_COVERAGE")
    print('\t'.join(cont))

    total_unmapped = 0
    total_mapped = 0
    total_len = 0
    main_mapped = 0
    main_unmapped = 0
    main_len = 0
    main_len_adj = 0
    cmd = samtools + " idxstats " + bam
    cmd_cont = run_cmd(cmd)
    x_chrom_total_reads = 0
    y_chrom_total_reads = 0
    seqver = ''

    for line in cmd_cont.strip().split('\n'):
        arr = line.split('\t')
        #print(arr)
        mapped = int(arr[2])
        unmapped = int(arr[3])
        total_line = mapped+unmapped 
        total_mapped += mapped
        if arr[0][:3] == "chr":
            chrom = arr[0][3:]
        else:
            chrom = arr[0]

        if chrom == "*":
            total_unmapped += unmapped
        chrom_len = int(arr[1])
        if chrom == '1':
            seqver = SEQVER_CHR1LEN[chrom_len]
        try:
            chrom_len_adj = MAPPABLE_LEN[seqver][chrom]
        except KeyError:
            chrom_len_adj = chrom_len
            pass
        total_len += chrom_len

        if chrom_len == 0:
            cov = 0
            cov_adj = 0
        else:
            cov = (read_len * mapped) / chrom_len
            if chrom_len_adj == 0:
                cov_adj = 0
            else:
                cov_adj = (read_len * mapped) / chrom_len_adj

        if len(chrom) <= 2:
            main_mapped += mapped
            if chrom == "*":
                main_unmapped += unmapped
            main_len += chrom_len
            main_len_adj += chrom_len_adj

        

        cont = []
        cont.append(chrom)
        cont.append(comma(chrom_len))
        cont.append(comma(mapped))
        cont.append(comma(unmapped))
        cont.append(comma(total_line))
        if total_line > 0:
            cont.append(str(round(100*mapped/total_line,2))+'%')
        else:
            cont.append("0%")
        cont.append(str(round(cov,1)))
        cont.append(str(round(cov_adj,1)))

        if chrom == "X":
            x_chrom_total_reads = total_line
        if chrom == "Y":
            y_chrom_total_reads = total_line

        print('\t'.join(cont))

    total_reads = total_unmapped + total_mapped
    cov = (read_len * total_mapped) / total_len
    exomcov = (read_len * total_mapped) / EXOM_SEQ_COVERED_REGION
    print ("================ALL CHROM=================")
    print("REFERENCE VERSION:",seqver)
    print("READ_LEN:",read_len)
    print ("MAPPED:",comma(total_mapped))
    print ("UNMAPPED:",comma(total_unmapped))
    print ("TOTAL READS:",comma(total_reads))
    print ("TOTAL LENGTH:",comma(total_len))
    print ("COVERAGE:",round(cov,1))
    print ("EXOM COVERAGE:",round(exomcov,1))

    main_reads = main_mapped + main_unmapped
    if main_len == 0:
        main_cov = 0.0
        main_cov_adj = 0.0
    else:
        main_cov = (read_len * main_mapped) / main_len
        main_cov_adj = (read_len * main_mapped) / main_len_adj
    
    main_exomcov = (read_len * main_mapped) / EXOM_SEQ_COVERED_REGION
    print ("=======ONLY MAIN CHROM (1~22,X,Y,MT)=======")
    print ("MAPPED:",comma(main_mapped), " (",round(100*main_mapped/total_mapped,1),"% )")
    #if total_unmapped == 0:
    #   print ("UNMAPPED:",comma(main_unmapped), " (",round(100*0,1),"% )")
    #else:
    #   print ("UNMAPPED:",comma(main_unmapped), " (",round(100*main_unmapped/total_unmapped,1),"% )")
    #print ("TOTAL READS:",comma(main_reads) , " (",round(100*main_reads/total_reads,1),"% )")
    print ("TOTAL LENGTH:",comma(main_len) , " (",round(100*main_len/total_len,1),"% )")
    print ("COVERAGE:",round(main_cov,1))
    print ("ADJUSTED COVERAGE:",round(main_cov_adj,1))
    print ("EXOM COVERAGE:",round(main_exomcov,1))
    print ("=======X,Y CHROM=======")
    print ("chrX READS:",comma(x_chrom_total_reads) + ' (' + str(round(100*x_chrom_total_reads/total_reads,3)) + '%)')
    print ("chrY READS:",comma(y_chrom_total_reads) + ' (' + str(round(100*y_chrom_total_reads/total_reads,3)) + '%)')
    if y_chrom_total_reads <= 0:
        print ("X/Y RATIO: 0")
        print ("EXT. GENDER: NA")
    else:
        print ("X/Y RATIO:",round(x_chrom_total_reads/y_chrom_total_reads,2))
        
        if x_chrom_total_reads/y_chrom_total_reads<200:
            gender = "Male"
        else:
            gender = "Female"
        print ("EXT. GENDER:",gender)

if __name__ == "__main__":
    bam = sys.argv[1]
    bamstat(bam)
