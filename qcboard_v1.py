#!/usr/bin/env python
# -*- coding: utf-8 -*-
#### qcboard_v1.py
#### made by Daniel Minseok Kwon
#### 2019-09-05 11:03:02
#########################
import sys
import os
import argparse

VERSION = "0.1"
VERSION_DATE = "19.09.05"

FIELDLIST_BAMQC = {}
FIELDLIST_BAMQC['PICARD.CMM'] = ['NO_PAIR','NO_PAIR_1','NO_PAIR_2']
FIELDLIST_BAMQC['FASTQC'] = ['SEQUENCE_LENGTH','GC_PER','DUPLICATION_PER']
FIELDLIST_BAMQC['SAMTOOLS'] = ['NO_MAPPED_READS','NO_UNMAPPED_READS','MAPPED_RATIO','XY_RATIO','EST_GENDER','CHROM_COVERAGE_TAB','COVERAGE_ALL_CHROM','COVERAGE_MAIN_CHROM']

def run_cmd(scmd, flag=False):
    if flag:
        print (scmd)
    rst = os.popen(scmd)
    rst_cont = rst.read()
    return rst_cont

def fileOpen(path):
    f = open(path, "r")
    return f.read() 

def fileSave (path, cont, opt, gzip_flag = "n"):
    if gzip_flag == "gz":
        import gzip
        f = gzip.open(path, opt)
        f.write(cont)
        f.close()
    else:
        f = open (path, opt)
        f.write(cont)
        f.close

def comma(value):
    return "{:,}".format(value)


def gzopen(fname):
    if fname.endswith(".gz") or fname.endswith(".zip"):
        import gzip
        f1 = gzip.GzipFile(fname, "r")
    else:
        f1 = open(fname)
    return f1

def get_options():
    parser = argparse.ArgumentParser(usage='%(prog)s <sub-command> [options]', description='%(prog)s ver' + VERSION + " (" + VERSION_DATE + ")" + ': convert bam to image')
    parser.add_argument('-v', '--version', action='version', version="%(prog)s ver" + VERSION + " (" + VERSION_DATE + ")")
    
    parser.add_argument('-bam', dest='bam', default='', help='bam file')
    parser.add_argument('-out', dest='out', default='', help='title of output file')
    parser.add_argument('-temp', dest='temp', default='qcboard_bamqc.html', help='template html file')
    parser.add_argument('-silence', dest='silence', action="store_true", default=False, help='don\'t print any log.')
    parser.add_argument('-debug', dest='debug', action="store_true", default=False, help='turn on the debugging mode')

    if len(sys.argv) == 1 or (len(sys.argv) == 2 and sys.argv[1][0] != '-'):
        sys.argv.append('-h')
    opt = vars(parser.parse_args())
    return opt


def cli():
    opt = get_options()
    qcb = QCBoard(opt)
    qcb.run()


class QCBoard():
    has_opt_error = False
    out_html = ""
    out_json = ""
    mode = "BAMQC" ## BAMQC or VCFQC

    def __init__(self, opt):
        self.opt = opt
        self.set_outfilenames()
        self.set_infilenames()
        self.qcstat = {}
        self.qcstat['SAMTOOLS'] = {}
        self.qcstat['PICARD.CMM'] = {}
        self.qcstat['FASTQC'] = {}
        
    def set_infilenames(self):
        self.infile = {}
        # self.infile['samtools.idxstats'] = self.opt['bam'] + '.idxstats'
        self.infile['samtools.idxstats'] = self.opt['bam'] + '.stat'
        self.infile['picard.cmm.alignment_summary_metrics'] = self.opt['bam'] + '.cmm.alignment_summary_metrics'
        self.infile['fastqc'] = self.opt['bam'][:-4] + '_fastqc.zip'


    def set_outfilenames(self):
        if self.opt['out'].endswith('.html'):
            self.out_html = self.opt['out']
            self.out_json = self.opt['out'][:-5] + '.json'
        else:
            self.out_html = self.opt['out'] + '.html'
            self.out_json = self.opt['out'] + '.json'

    def set_refeq_from_ucsc(self, pos1):
        spos = pos1['g_spos']-self.opt['margin'] - 500
        epos = pos1['g_epos']+self.opt['margin']+1 + 500
        url = "http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=chr"+pos1['chrom']+":"+str(spos+1)+","+str(epos+1)
        cont = get_url(url)
        seq = ""
        for line in cont.strip().split('\n'):
            if line[0] != '<':
                seq += line.strip().upper()
        i = 0
        for gpos in range(spos,epos):
            self.refseq[gpos] = seq[i]
            i += 1

    def set_refseq_from_localfasta(self, pos1):
        spos = pos1['g_spos']-self.opt['margin'] - 500
        epos = pos1['g_epos']+self.opt['margin']+1 + 500
        seq = self.get_refseq(self.opt['ref'], pos1['chrom'], spos,epos)
        i = 0
        for gpos in range(spos, epos):
            self.refseq[gpos] = seq[i]
            i += 1
    
    def get_refseq(self, ref, chrom, spos,epos):
        f = Fasta(ref)
        fastachrommap = {}
        for c1 in list(f.keys()):
            arr = c1.split(' ')
            tchrom = arr[0]
            fastachrommap[tchrom] = c1
        # fasta_chrom = chrom + ' dna:chromosome chromosome:GRCh37:'+chrom+':1:'+str(CHROM_LEN['b37d5'][chrom])+':1'
        refseq = f[fastachrommap[chrom]][spos:epos+1]
        return refseq

    
    def save_html(self):
        cont = fileOpen(self.opt['temp'])
        cont = cont.replace('##BAMFILE##',self.opt['bam'])

        for k1 in FIELDLIST_BAMQC.keys():
            if len(self.qcstat[k1].keys()) > 0:
                for k2 in FIELDLIST_BAMQC[k1]:
                    if type(self.qcstat[k1][k2]) == type(1) or type(self.qcstat[k1][k2]) == type(1.1):
                        v1 = comma(self.qcstat[k1][k2])
                    else:
                        v1 = self.qcstat[k1][k2]
                    cont = cont.replace('##'+k1+'.'+k2+'##', v1)

        fileSave(self.out_html, cont, 'w')
        print ('Saved '+self.out_html)


    def mk_run_script(self):
        pass

    def mk_output(self):
        pass

    def get_samtools_idxstats(self):
        flag = ""
        chromcovtab = ""

        for line in open(self.infile['samtools.idxstats']):
            arr = line.split('\t')
            arr[-1] = arr[-1].strip()
            if arr[0] == "CHROM":
                flag = "chrcov"

            if flag == "chrcov":
                if arr[0] == "CHROM":
                    chromcovtab += '<tr><th>'+'</th><th>'.join(arr)+'</th></tr>'
                else:
                    chromcovtab += '<tr><td>'+'</td><td>'.join(arr)+'</td></tr>'
                if arr[0] == "MT":
                    flag = ""
            
            if flag=="" and "=ALL CHROM=" in line:
                flag = "ALL"
            if flag=="ALL" and "=ONLY MAIN CHROM" in line:
                flag = "MAIN"
            if flag=="MAIN" and "=X,Y CHROM=" in line:
                flag = "XY"

            if flag == "ALL":
                if line[:len("COVERAGE:")] == "COVERAGE:":
                    self.qcstat['SAMTOOLS']['COVERAGE_ALL_CHROM'] = line.replace("COVERAGE:","").strip()
                if line[:len("MAPPED:")] == "MAPPED:":
                    self.qcstat['SAMTOOLS']['NO_MAPPED_READS'] = line.replace("MAPPED:","").strip()
                if line[:len("UNMAPPED:")] == "UNMAPPED:":
                    self.qcstat['SAMTOOLS']['NO_UNMAPPED_READS'] = line.replace("UNMAPPED:","").strip()

            if flag == "MAIN":
                if line[:len("COVERAGE:")] == "COVERAGE:":
                    self.qcstat['SAMTOOLS']['COVERAGE_MAIN_CHROM'] = line.replace("COVERAGE:","").strip()
            if flag == "XY":
                if line[:len("X/Y RATIO:")] == "X/Y RATIO:":
                    self.qcstat['SAMTOOLS']['XY_RATIO'] = line.replace("X/Y RATIO:","").strip()
                if line[:len("EXT. GENDER:")] == "EXT. GENDER:":
                    self.qcstat['SAMTOOLS']['EST_GENDER'] = line.replace("EXT. GENDER:","").strip()
            print (arr)

        self.qcstat['SAMTOOLS']['CHROM_COVERAGE_TAB'] = chromcovtab
        self.qcstat['SAMTOOLS']['MAPPED_RATIO'] = round(100*int(self.qcstat['SAMTOOLS']['NO_MAPPED_READS'].replace(',','')) / (int(self.qcstat['SAMTOOLS']['NO_MAPPED_READS'].replace(',','')) + int(self.qcstat['SAMTOOLS']['NO_UNMAPPED_READS'].replace(',',''))), 3)

    
    def get_picard_cmm(self):
        for line in open(self.infile['picard.cmm.alignment_summary_metrics']):
            if line[0] != '#':
                arr = line.split('\t')
                if arr[0] == "FIRST_OF_PAIR":
                    self.qcstat['PICARD.CMM']['NO_PAIR_1'] = int(arr[1])
                if arr[0] == "SECOND_OF_PAIR":
                    self.qcstat['PICARD.CMM']['NO_PAIR_2'] = int(arr[1])
                if arr[0] == "PAIR":
                    self.qcstat['PICARD.CMM']['NO_PAIR'] = int(arr[1])

    def get_fastqc(self):
        cmd = "unzip -f " + self.infile['fastqc'] + " -d " + '/'.join(self.infile['fastqc'].split('/')[:-1])
        run_cmd(cmd)
        print (self.infile['fastqc'][:-4] + '/fastqc_data.txt')
        for line in open(self.infile['fastqc'][:-4] + '/fastqc_data.txt'):
            if line[:len('Sequence length')] == 'Sequence length':
                self.qcstat['FASTQC']['SEQUENCE_LENGTH'] = line.split('\t')[1].strip()
            if line[:len('%GC')] == '%GC':
                self.qcstat['FASTQC']['GC_PER'] = line.split('\t')[1].strip()
            if line[:len('#Total Deduplicated Percentage')] == '#Total Deduplicated Percentage':
                self.qcstat['FASTQC']['DUPLICATION_PER'] = str(round(float(line.split('\t')[1].strip()),3))
            
    def run(self):
        self.get_samtools_idxstats()
        self.get_picard_cmm()
        self.get_fastqc()
        self.save_html()
        


if __name__ == "__main__":
    print ("#USAGE: python qcboard_v1.py -bam [BAM FILE] -out [OUTPUT TITLE]")
    cli()
