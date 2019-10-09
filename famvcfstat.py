#!/usr/bin/env python
# -*- coding: utf-8 -*-
import qcutil
import os.path
import conf

UNITREADLINE = 20000

class QCBoardFamilyVCFSTAT():
    def __init__(self, opt):
        self.opt = opt
        self.vcflist = self.opt['vcflist']
        self.set_outfilenames()
        self.qcstat = self.init_qcstat()

    def set_outfilenames(self):
        if 'out' in self.opt.keys() and self.opt['out'].endswith('.json'):
            self.out_json = self.opt['out'][:-5] + '.famstat.json'
        elif 'out' in self.opt.keys() and self.opt['out'] != '':
            self.out_json = self.opt['out'] + '.famstat.json'
        else:
            self.out_json = self.opt['vcf'] + '.famstat.json'

    def init_qcstat(self):
        stat = {}
        for k1 in list(conf.FIELDLIST_VCFQC.keys()):
            stat[k1] = {}
            for k2 in conf.FIELDLIST_VCFQC[k1]:
                stat[k1][k2] = 0
                stat[k1][k2+'_PASS'] = 0
        return stat

    def get_read_block_vcf(self, fp, tchrom, tpos, max_read_line=0, flag_eof=False):
        blockvar = {}
        i = 0
        chrom = ''
        pos = 0
        while True:
            line = fp.readline().decode('UTF-8')
            try:
                tmp = line[0]
            except IndexError:
                break
                flag_eof = True
            if line[0] != '#':
                i += 1
                arr = line.split('\t')
                arr[-1] = arr[-1].strip()
                chrom = arr[0].replace('chr','').replace('MT','M')
                pos = int(arr[1])
                vqsr_filter = arr[conf.VCF_HEADCOL.index('FILTER')].strip()
                ref = arr[conf.VCF_HEADCOL.index('REF')].strip()
                alt = arr[conf.VCF_HEADCOL.index('ALT')].strip()
                
                k = 9
                callinfo = arr[k].strip().split(':')
                gt = callinfo[0]
                
                homhet = qcutil.get_homhet(gt)
                (adj_ref, adj_alt, adj_gt) = qcutil.get_adjusted_refalt(ref,alt,gt)
                vartype = qcutil.get_vartype(adj_ref,adj_alt)

                if vartype == "SNV":
                    d = {}
                    d['geno'] = qcutil.get_genotype(adj_gt,adj_ref, adj_alt)
                    # key1 = '_'.join([pos,adj_ref,adj_alt])
                    # print (pos, adj_ref, adj_alt,ref, alt, geno)
                    if vqsr_filter == "PASS":
                        d['vqsr'] = 'p'
                    else:
                        d['vqsr'] = ''
                    d['ref'] = ref
                    blockvar[pos] = d
                if not flag_eof:
                    if max_read_line>0:
                        if i > max_read_line:
                            break
                    else:
                        if pos >= tpos or chrom != tchrom:
                            break
        if max_read_line>0:
            tpos = pos
            tchrom = chrom
        return blockvar, tchrom, tpos, flag_eof

    def cal_mendelstat_from_block(self, mstat, blockC, blockF, blockM):
        for p1 in list(blockC.keys()):
            bC = blockC[p1]
            gC = bC['geno']
            rr = bC['ref']+bC['ref']
            flag = True
            try:
                gF = blockF[p1]['geno']
            except KeyError:
                gF = rr
                # flag = False

            try:
                gM = blockM[p1]['geno']
            except KeyError:
                gM = rr
                # flag = False

            if flag:
                vt = ''
                if bC['vqsr'] == "p":
                    vt = '_PASS'

                mstat['N'+vt] += 1
                if gC != gF and gC != gM:
                    if (gC[0] in gF and gC[1] in gM) or (gC[1] in gF and gC[0] in gM):
                        pass
                    else:
                        mstat['MENDEL_ERROR'+vt] += 1
                        if gF == rr and gM == rr:
                            mstat['DENOVO'+vt] += 1
                        # print (gC, gF, gM)
        return mstat

    def load_vcflist(self):
        global UNITREADLINE
        stat = self.qcstat
        # print (stat)
        
        flist = []
        for vcf in self.vcflist:
            flist.append(qcutil.gzopen(vcf))

        mendelstat = {}
        blockvar = {}
        blockvarchild = {}
        tchrom = ''
        tpos = 0

        iblock = 0
        while True:
            iblock += 1
            blockvar[self.pedmap['C']], tchrom, tpos, flag_eof = self.get_read_block_vcf(flist[self.pedmap['C']], tchrom,tpos, UNITREADLINE)
            blockvar[self.pedmap['F']], tchrom, tpos, flag_eof = self.get_read_block_vcf(flist[self.pedmap['F']], tchrom,tpos, UNITREADLINE, flag_eof)
            blockvar[self.pedmap['M']], tchrom, tpos, flag_eof = self.get_read_block_vcf(flist[self.pedmap['M']], tchrom,tpos, UNITREADLINE, flag_eof)
            
            mendelstat = qcutil.check_key(mendelstat, 'MENDEL',{'N':0,'MENDEL_ERROR':0,'DENOVO':0,'MENDEL_ERROR_RATIO':0.0,'DENOVO_RATIO':0.0})
            mendelstat['MENDEL'] = self.cal_mendelstat_from_block(mendelstat['MENDEL'],blockvar[self.pedmap['C']], blockvar[self.pedmap['F']], blockvar[self.pedmap['M']])
            print (iblock, len(blockvar[self.pedmap['C']].keys()), mendelstat, tchrom, tpos)
            if iblock > 2000:
                pass
                break

        self.qcstat = mendelstat
        self.cal_qcstat()

    def cal_qcstat(self):
        stat = self.qcstat
        print (stat)
        # stat['QCBOARD']['NOSAMPLE'] = len(self.vcflist)
        for vt in ['','_PASS']:
            if stat['MENDEL']['N'+vt] == 0:
                stat['MENDEL']['MENDEL_ERROR_RATIO'+vt] = 0
                stat['MENDEL']['DENOVO_RATIO'+vt] = 0
            else:
                stat['MENDEL']['MENDEL_ERROR_RATIO'+vt] = round(stat['MENDEL']['MENDEL_ERROR'+vt] / stat['MENDEL']['N'+vt],4)
                stat['MENDEL']['DENOVO_RATIO'+vt] = round(stat['MENDEL']['DENOVO'+vt] / stat['MENDEL']['N'+vt],4)
        self.qcstat = stat

    def save_stat(self):
        qcutil.savejson(self.out_json, self.qcstat)
        print ('Saved '+self.out_json)

    def load_ped(self):
        pedlist = self.opt['ped'].strip().split(',')
        self.pedmap = {'F':'','M':'','C':''}
        for i in range(len(pedlist)):
            ped = pedlist[i].strip()
            if ped == 'F' or ped == 'M' or ped == 'C':
                self.pedmap[ped]=i

    def run(self):
        self.load_ped()
        self.load_vcflist()
        self.save_stat()