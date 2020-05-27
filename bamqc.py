#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os.path
import json
import qcutil
import conf

FIELDLIST_BAMQC = {}
FIELDLIST_BAMQC['PICARD.CMM'] = ['NO_PAIR','NO_PAIR_1','NO_PAIR_2','MEDIAN_INSERT_SIZE','MIN_INSERT_SIZE','MAX_INSERT_SIZE','MEAN_INSERT_SIZE','STANDARD_DEVIATION']
FIELDLIST_BAMQC['FASTQC'] = ['GC_PER','DUPLICATION_PER']

# FIELDLIST_BAMQC['FASTQC'].append('Basic_Statistics')
# FIELDLIST_BAMQC['FASTQC'].append('Per_base_sequence_quality')
# FIELDLIST_BAMQC['FASTQC'].append('Per_tile_sequence_quality')
# FIELDLIST_BAMQC['FASTQC'].append('Per_sequence_quality_scores')
# FIELDLIST_BAMQC['FASTQC'].append('Per_base_sequence_content')
# FIELDLIST_BAMQC['FASTQC'].append('Per_sequence_GC_content')
# FIELDLIST_BAMQC['FASTQC'].append('Per_base_N_content')
# FIELDLIST_BAMQC['FASTQC'].append('Sequence_Length_Distribution')
# FIELDLIST_BAMQC['FASTQC'].append('Sequence_Duplication_Levels')
# FIELDLIST_BAMQC['FASTQC'].append('Overrepresented_sequences')
# FIELDLIST_BAMQC['FASTQC'].append('Adapter_Content')
# FIELDLIST_BAMQC['FASTQC'].append('Kmer_Content')

FIELDLIST_BAMQC['SAMTOOLS'] = ['SEQUENCE_LENGTH','NO_MAPPED_READS','NO_UNMAPPED_READS','MAPPED_RATIO','XY_RATIO','EST_GENDER','CHROM_COVERAGE_TAB','COVERAGE_ALL_CHROM','COVERAGE_MAIN_CHROM','SEQVERSION','ADJUSTED_COVERAGE_MAIN_CHROM']
CHROMCOVTAB_HEADERMAP = {'CHROM':'Chromosome','LEN':'Length','MAAPED':'Num Mapped','UNMAPPED':'Num Unmapped','TOTAL':'Total','MAPPED_RATIO':'Mapped Ratio','COVERAGE':'Coverage','ADJUSTED_COVERAGE':'Adjusted Coverage'}


class QCBoardBAM():
    has_opt_error = False
    out_html = ""
    out_json = ""

    def __init__(self, opt):
        self.opt = opt
        self.set_outfilenames()
        self.set_infilenames()
        self.qcstat = {}
        self.qcstat['SAMTOOLS'] = {}
        self.qcstat['PICARD.CMM'] = {}
        self.qcstat['FASTQC'] = {}
        self.qcstat['FASTQC']['GC_PER'] = self.opt['per_gc']
        self.qcstat['FASTQC']['DUPLICATION_PER'] = self.opt['per_dup']

    def set_infilenames(self):
        self.infile = {}
        # self.infile['samtools.idxstats'] = self.opt['bam'] + '.idxstats'
        self.infile['samtools.idxstats'] = self.opt['bam'] + '.stat'
        self.infile['picard.cmm.alignment_summary_metrics'] = self.opt['bam'] + '.cmm.alignment_summary_metrics'
        self.infile['picard.cmm.quality_by_cycle_metrics'] = self.opt['bam'] + '.cmm.quality_by_cycle_metrics'
        self.infile['picard.cmm.insert_size_metrics'] = self.opt['bam'] + '.cmm.insert_size_metrics'
        self.infile['fastqc'] = self.opt['bam'][:-4] + '_fastqc.zip'


    def set_outfilenames(self):
        if self.opt['out'].endswith('.html'):
            self.out_html = self.opt['out']
            self.out_json = self.opt['out'][:-5] + '.json'
            self.out_tsv = self.opt['out'][:-5] + '.tsv'
        else:
            self.out_html = self.opt['out'] + '.html'
            self.out_json = self.opt['out'] + '.json'
            self.out_tsv = self.opt['out'] + '.tsv'

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


    def str_num(self, x):
        subst_x = str(x).replace(',', '').replace('%', '')
        try:
            return int(subst_x)
        except:
            try:
                return float(subst_x)
            except:
                return subst_x

    def save_html(self):
        cont = qcutil.fileOpen(os.path.join(conf.DIRPATH, self.opt['temp']) )
        cont = cont.replace('##BAMFILE##',self.opt['bam'])

        for k1 in FIELDLIST_BAMQC.keys():
            if len(self.qcstat[k1].keys()) > 0:
                for k2 in FIELDLIST_BAMQC[k1]:
                    if k2 == 'CHROM_COVERAGE_TAB':
                        k2 = 'CHROM_COVERAGE_TAB_HTML'
                    if type(self.qcstat[k1][k2]) == type(1) or type(self.qcstat[k1][k2]) == type(1.1):
                        v1 = qcutil.comma(self.qcstat[k1][k2])
                    else:
                        v1 = self.qcstat[k1][k2]

                    if k1 == "FASTQC":
                        v1 = v1.replace('PASS','<span class="badge badge-success">PASS</span>')
                        v1 = v1.replace('WARN','<span class="badge badge-warning">WARN</span>')
                        v1 = v1.replace('FAIL','<span class="badge badge-danger">FAIL</span>')

                    cont = cont.replace('##'+k1+'.'+k2+'##', v1)

        qcutil.fileSave(self.out_html, cont, 'w')
        print ('Saved '+self.out_html)

    def save_json(self):
        conv_dict = {
            'SEQVERSION': 'Reference version',
            'EST_GENDER': 'Estimated gender',
            'COVERAGE_ALL_CHROM': 'Coverage all chromosomes',
            'COVERAGE_MAIN_CHROM': 'Coverage main chromosomes',
            'ADJUSTED_COVERAGE_MAIN_CHROM': 'Adjusted coverage main chromosomes',
            'SEQUENCE_LENGTH': 'Sequence length',
            'XY_RATIO': 'XY ratio',
            'NO_UNMAPPED_READS':'Num unmapped reads',
            'NO_MAPPED_READS': 'Num mapped reads',
            'NO_PAIR': 'Num pairs',
            'NO_PAIR_1': 'Num 1st of pair',
            'NO_PAIR_2': 'Num 2nd of pair',
            'MEDIAN_INSERT_SIZE': 'Median insertion size',
            'MEDIAN_ABSOLUTE_DEVIATION': 'Median absolute deviation of insertion size distribution',
            'MEAN_INSERT_SIZE': 'Mean insertion size',
            'STANDARD_DEVIATION': 'Standard deviation of insertion size',
            'MIN_INSERT_SIZE': 'Minimum insertion size',
            'MAX_INSERT_SIZE': 'Maximum insertion size',
            'MAPPED_RATIO': 'Mapped ratio'
        }
        d = {}
        for key in self.qcstat:
            if key == 'SAMTOOLS':
                # parsing samtools info
                for k, v in self.qcstat['SAMTOOLS'].items():
                    if k == 'CHROM_COVERAGE_TAB_HTML':
                        pass
                    elif k == 'CHROM_COVERAGE_TAB':
                        d.setdefault('Per chromosome coverage', [])
                        for chr_str in v.split('|')[1:-1]:
                            chr, len, map, unmap, tot, map_ratio, coverage, adj_coverage = chr_str.split(':')
                            d_chr = {'Chromosome': self.str_num(chr), 'Chromosome length': self.str_num(len), 'Num mapped reads': self.str_num(map),
                                     'Num unmapped reads': self.str_num(unmap), 'Total reads': self.str_num(tot),
                                     'Mapped ratio': self.str_num(map_ratio), 'Coverage': self.str_num(coverage),'Adjusted coverage': self.str_num(adj_coverage)
                                    }
                            d['Per chromosome coverage'].append(d_chr)
                    else:
                        d.setdefault(conv_dict[k], self.str_num(v))
            elif key == 'PICARD.CMM':
                for k, v in self.qcstat['PICARD.CMM'].items():
                    d.setdefault(conv_dict[k], self.str_num(v))

        # writing json
        with open(self.out_json, 'w') as outfile:
            json.dump(d, outfile)
        print ('Saved '+self.out_json)

    def mk_run_script(self):
        pass

    def mk_output(self):
        pass

    def get_samtools_idxstats(self):
        global CHROMCOVTAB_HEADERMAP
        flag = ""
        chromcovtab_html = ""
        chromcovtab = ""

        for line in open(self.infile['samtools.idxstats']):
            arr = line.split('\t')
            arr[-1] = arr[-1].strip()
            if arr[0] == "CHROM":
                flag = "chrcov"

            if flag == "chrcov":
                if arr[0] == "CHROM":
                    chromcovtab_html += '<tr>'
                    for a1 in arr:
                        chromcovtab_html += '<th>'+CHROMCOVTAB_HEADERMAP[a1]+'</th>'
                        if a1 != 'CHROM':
                            chromcovtab += ':'
                        chromcovtab += CHROMCOVTAB_HEADERMAP[a1]
                    chromcovtab_html += '</tr>'
                    chromcovtab += '|'
                else:
                    chromcovtab_html += '<tr><td>'+'</td><td>'.join(arr)+'</td></tr>'
                    chromcovtab += ':'.join(arr) + '|'
                if arr[0] == "MT" or arr[0] == "M":
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
                if line[:len("READ_LEN:")] == "READ_LEN:":
                    self.qcstat['SAMTOOLS']['SEQUENCE_LENGTH'] = line.replace("READ_LEN:","").strip()
                if line[:len("REFERENCE VERSION:")] == "REFERENCE VERSION:":
                    self.qcstat['SAMTOOLS']['SEQVERSION'] = line.replace("REFERENCE VERSION:","").strip()

            if flag == "MAIN":
                if line[:len("COVERAGE:")] == "COVERAGE:":
                    self.qcstat['SAMTOOLS']['COVERAGE_MAIN_CHROM'] = line.replace("COVERAGE:","").strip()
                if line[:len("ADJUSTED COVERAGE:")] == "ADJUSTED COVERAGE:":
                    self.qcstat['SAMTOOLS']['ADJUSTED_COVERAGE_MAIN_CHROM'] = line.replace("ADJUSTED COVERAGE:","").strip()
            if flag == "XY":
                if line[:len("X/Y RATIO:")] == "X/Y RATIO:":
                    self.qcstat['SAMTOOLS']['XY_RATIO'] = line.replace("X/Y RATIO:","").strip()
                if line[:len("EXT. GENDER:")] == "EXT. GENDER:":
                    self.qcstat['SAMTOOLS']['EST_GENDER'] = line.replace("EXT. GENDER:","").strip()
            # print (arr)

        self.qcstat['SAMTOOLS']['CHROM_COVERAGE_TAB_HTML'] = chromcovtab_html
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

        colnames = []
        for line in open(self.infile['picard.cmm.insert_size_metrics']):
            if line[0] != '#':
                arr = line.split('\t')
                arr[-1] = arr[-1].strip()
                if len(colnames) > 3:
                    for k in range(len(colnames)):
                        if colnames[k] in FIELDLIST_BAMQC['PICARD.CMM']:
                            if '.' in arr[k]:
                                try:
                                    self.qcstat['PICARD.CMM'][colnames[k]] = float(arr[k].strip())
                                except ValueError:
                                    self.qcstat['PICARD.CMM'][colnames[k]] = arr[k].strip()
                            else:
                                try:
                                    self.qcstat['PICARD.CMM'][colnames[k]] = int(arr[k])
                                except ValueError:
                                    self.qcstat['PICARD.CMM'][colnames[k]] = arr[k].strip()
                    break
                if arr[0] == "MEDIAN_INSERT_SIZE":
                    colnames = arr

    def get_fastqc(self):
        # cmd = "unzip -f " + self.infile['fastqc'] + " -d " + '/'.join(self.infile['fastqc'].split('/')[:-1])
        cmd = "unzip " + self.infile['fastqc'] + " -d " + '/'.join(self.infile['fastqc'].split('/')[:-1])
        qcutil.run_cmd(cmd)
        print (self.infile['fastqc'][:-4] + '/fastqc_data.txt')
        for line in open(self.infile['fastqc'][:-4] + '/fastqc_data.txt'):
            if line[:len('Sequence length')] == 'Sequence length':
                self.qcstat['FASTQC']['SEQUENCE_LENGTH'] = line.split('\t')[1].strip()
            if line[:len('%GC')] == '%GC':
                self.qcstat['FASTQC']['GC_PER'] = line.split('\t')[1].strip()
            if line[:len('#Total Deduplicated Percentage')] == '#Total Deduplicated Percentage':
                self.qcstat['FASTQC']['DUPLICATION_PER'] = str(round(float(line.split('\t')[1].strip()),3))

        print (self.infile['fastqc'][:-4] + '/summary.txt')
        for line in open(self.infile['fastqc'][:-4] + '/summary.txt'):
            line = line.strip()
            arr = line.split('\t')
            self.qcstat['FASTQC'][arr[1].strip().replace(' ','_')] = arr[0].strip()

    def run(self):
        self.get_samtools_idxstats()
        self.get_picard_cmm()
        # self.get_fastqc()
        self.save_html()
        self.save_json()