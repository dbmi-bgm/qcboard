#!/usr/bin/env python
# -*- coding: utf-8 -*-
import qcutil
import os.path
import conf

UNITREADLINE = 20000

MENDELSTAT_DEFAULT_FIELD = ['N','MENDEL_ERROR','DENOVO','MENDEL_INHERITANCE']
VARTYPE_LIST = ['ALL','SNV', 'INS','DEL', 'MNV', 'MIXED']

class QCBoardFamilyVCFSTAT():
    def __init__(self, opt):
        self.relation = []
        self.sample_relation_map = {}
        self.relation_sample_map = {}
        self.opt = opt
        self.vcflist = self.opt['vcf']
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

    def get_read_block_from_one_vcf(self, fp, tchrom, tpos, max_read_line=0, flag_eof=False):
        blockvar = {}
        i = 0
        chrom = ''
        pos = 0
        samplelist = []

        while True:
            line = fp.readline().decode('UTF-8')

            try:
                line[0]
            except IndexError:
                flag_eof = True
                break

            if line[0] == '#':
                if line[:len('#CHROM')] == '#CHROM':
                    arr = line.strip().split('\t')
                    samplelist = arr[9:]
            else:
                i += 1
                arr = line.split('\t')
                arr[-1] = arr[-1].strip()
                chrom = qcutil.convert_chrom(arr[0])
                pos = int(arr[1])
                chrompos = chrom + ':' + str(pos)
                vqsr_filter = arr[conf.VCF_HEADCOL.index('FILTER')].strip()
                ref = arr[conf.VCF_HEADCOL.index('REF')].strip()
                alt = arr[conf.VCF_HEADCOL.index('ALT')].strip()
                
                for sidx in range(len(samplelist)):
                    callinfo = arr[sidx+9].strip().split(':')
                    gt = callinfo[0]
                    
                    homhet = qcutil.get_homhet(gt)
                    (adj_ref, adj_alt, adj_gt) = qcutil.get_adjusted_refalt(ref,alt,gt)
                    vartype = qcutil.get_vartype(adj_ref,adj_alt)


                    d = {}
                    d['vartype'] = vartype
                    d['numgt'] = gt.replace('|','/')
                    d['geno'] = qcutil.get_genotype(adj_gt, adj_ref, adj_alt)  # useless
                    
                    if vqsr_filter == "PASS":
                        d['vqsr'] = 'p'
                    else:
                        d['vqsr'] = ''
                    d['ref'] = ref
                    try:
                        blockvar[samplelist[sidx]]
                    except KeyError:
                        blockvar[samplelist[sidx]] = {}
                    blockvar[samplelist[sidx]][chrompos] = d

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
        

    def cal_mendelstat_from_block(self, mstat, blockvar):
        probandID = self.relation_sample_map['Proband']
        fatherID = self.relation_sample_map['Father']
        motherID = self.relation_sample_map['Mother']
        for p1 in blockvar[probandID].keys():
            bP = blockvar[probandID][p1]
            vartype = bP['vartype']
            gP = bP['numgt']
            rr = "0/0"
            flag = True
            try:
                gF = blockvar[fatherID][p1]['numgt']
            except KeyError:
                gF = rr

            try:
                gM = blockvar[motherID][p1]['numgt']
            except KeyError:
                gM = rr

            gPa = gP.split('/')
            gFa = gF.split('/')
            gMa = gM.split('/')

            if flag:
                vfilter = ''
                if bP['vqsr'] == "p":
                    vfilter = '_PASS'

                for vt in ['ALL_', vartype + '_']:
                    mstat[vt + 'N'+vfilter] += 1
                    flag_mendel_inheritance = True
                    if gP != gF and gP != gM:
                        if (gPa[0] in gFa and gPa[1] in gMa) or (gPa[1] in gFa and gPa[0] in gMa):
                            pass
                        else:
                            flag_mendel_inheritance = False
                            # print(p1, gP, gF, gM)
                            mstat[vt + 'MENDEL_ERROR'+vfilter] += 1
                            if gF == rr and gM == rr:
                                mstat[vt + 'DENOVO'+vfilter] += 1
                            # print (gP, gF, gM)
                    if flag_mendel_inheritance:
                        mstat[vt + 'MENDEL_INHERITANCE'+vfilter] += 1
        return mstat

    def get_mendelstat_default(self):
        global MENDELSTAT_DEFAULT_FIELD, VARTYPE_LIST
        d = {}
        for vt in VARTYPE_LIST:
            for f1 in MENDELSTAT_DEFAULT_FIELD:
                d[vt + '_' + f1] = 0
        return d

    def load_vcflist(self):
        global UNITREADLINE
        stat = self.qcstat
        
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
            if len(flist) > 1:
                # TODO: 
                # blockvar, tchrom, tpos, flag_eof = self.get_read_block_from_multiple_vcf(flist, tchrom,tpos, UNITREADLINE)
                pass
            else:
                blockvar, tchrom, tpos, flag_eof = self.get_read_block_from_one_vcf(flist[0], tchrom,tpos, UNITREADLINE, False)
            mendelstat = qcutil.check_key(mendelstat, 'MENDEL',self.get_mendelstat_default())
            mendelstat['MENDEL'] = self.cal_mendelstat_from_block(mendelstat['MENDEL'],blockvar)
            if iblock > 2000 or flag_eof:
                break

        self.qcstat = mendelstat
        self.cal_qcstat()

    def cal_qcstat(self):
        global MENDELSTAT_DEFAULT_FIELD, VARTYPE_LIST
        stat = self.qcstat
        for vt in VARTYPE_LIST:
            vt = vt + '_'
            for f1 in MENDELSTAT_DEFAULT_FIELD:
                for vfilter in ['','_PASS']:
                    if stat['MENDEL'][vt + f1+vfilter] == 0:
                        stat['MENDEL'][vt + f1 + vfilter + '_RATIO'] = 0.0
                    else:
                        stat['MENDEL'][vt + f1 + vfilter + '_RATIO'] = round(stat['MENDEL'][vt + f1+vfilter] / stat['MENDEL'][vt + 'N'+vfilter],4)
        self.qcstat = stat

    def save_stat(self):
        qcutil.savejson(self.out_json, self.qcstat, sort=True)
        print ('Saved '+self.out_json)

    def load_relation(self):
        self.relation = qcutil.loadjson(self.opt['relation'])
        for r1 in self.relation:
            self.sample_relation_map[r1['sample']] = r1
            try:
                tmp = self.relation_sample_map[r1['relation']]
                if type(tmp) == type([]):
                    self.relation_sample_map[r1['relation']].append(r1['sample'])
                else:
                    self.relation_sample_map[r1['relation']]= [tmp, r1['sample']]
            except KeyError:
                self.relation_sample_map[r1['relation']] = r1['sample']
        
    def run(self, flag_save=True):
        self.load_relation()
        self.load_vcflist()
        if flag_save:
            self.save_stat()