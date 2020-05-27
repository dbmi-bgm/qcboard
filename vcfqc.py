#!/usr/bin/env python
# -*- coding: utf-8 -*-
import qcutil
import os.path
import conf
import vcfstat
import famvcfstat

class QCBoardVCF():
    has_opt_error = False
    out_html = ""
    out_json = ""

    def __init__(self, opt):
        self.opt = opt
        self.vcflist = self.opt['vcf']
        self.set_outfilenames()
        self.set_infilenames()

    def set_infilenames(self):
        self.infile = {}
        for vcf in self.vcflist:
            self.infile[vcf] = {}
            self.infile[vcf]['vcfstat'] = vcf + '.stat.json'
            # self.infile[vcf]['ethnicpred'] = vcf + '.ethnic'

    def set_outfilenames(self):
        if self.opt['out'].endswith('.html'):
            self.out_html = self.opt['out']
            self.out_json = self.opt['out'][:-5] + '.json'
            self.out_tsv = self.opt['out'][:-5] + '.tsv'
            self.out_famstat = self.opt['out'][:-5] + '.famstat.json'
        else:
            self.out_html = self.opt['out'] + '_vcfqc.html'
            self.out_json = self.opt['out'] + '.json'
            self.out_tsv = self.opt['out'] + '.tsv'
            self.out_famstat = self.opt['out'] + '.famstat.json'

    def save_html(self):
        for vcf in self.vcflist:
            self.qcstat[vcf]['MAIN_CHROM_LIST'] = conf.MAIN_CHROM_LIST
            self.qcstat[vcf]['VCFFILE'] = vcf
            self.qcstat[vcf] = qcutil.add_comma_with_dict(self.qcstat[vcf])
            qcutil.render(self.out_html,self.opt['temp'],self.qcstat[vcf])
            print ('Saved '+self.out_html)

    def save_json(self):
        conv_dict = {
            'EST_GENDER': 'Estimated gender',
            'COVERAGE_ALL_CHROM': 'Coverage all chromosomes',
            'COVERAGE_MAIN_CHROM': 'Coverage main chromosomes',
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
        # for key in self.qcstat:
        #     if key == 'SAMTOOLS':
        #         # parsing samtools info
        #         for k, v in self.qcstat['SAMTOOLS'].items():
        #             if k == 'CHROM_COVERAGE_TAB_HTML':
        #                 pass
        #             elif k == 'CHROM_COVERAGE_TAB':
        #                 d.setdefault('Per chromosome coverage', [])
        #                 for chr_str in v.split('|')[1:-1]:
        #                     chr, len, map, unmap, tot, map_ratio, coverage = chr_str.split(':')
        #                     d_chr = {'Chromosome': self.str_num(chr), 'Chromosome length': self.str_num(len), 'Num mapped reads': self.str_num(map),
        #                                 'Num unmapped reads': self.str_num(unmap), 'Total reads': self.str_num(tot),
        #                                 'Mapped ratio': self.str_num(map_ratio), 'Coverage': self.str_num(coverage)
        #                             }
        #                     d['Per chromosome coverage'].append(d_chr)
        #             else:
        #                 d.setdefault(conv_dict[k], self.str_num(v))
        #     elif key == 'PICARD.CMM':
        #         for k, v in self.qcstat['PICARD.CMM'].items():
        #             d.setdefault(conv_dict[k], self.str_num(v))

        d = self.qcstat

        qcutil.savejson(self.out_json, d)
        print ('Saved '+self.out_json)


    def get_vcfstat(self):
        self.qcstat = {}
        for vcf in self.vcflist:
            if not qcutil.is_exist(self.infile[vcf]['vcfstat']):
                tmpopt = {}
                tmpopt['vcf'] = vcf
                qs = vcfstat.QCBoardVCFSTAT(tmpopt)
                qs.run()
            self.qcstat[vcf] = qcutil.loadjson(self.infile[vcf]['vcfstat'])

    def get_family_stat(self):
        qs = famvcfstat.QCBoardFamilyVCFSTAT(self.opt)
        qs.run()
        self.famqcstat = qcutil.loadjson(self.out_famstat)

    def run(self):
        self.get_vcfstat()
        if self.opt['relation'] != '':
            self.get_family_stat()
            
        self.save_html()
        self.save_json()
        