#!/usr/bin/env python
# -*- coding: utf-8 -*-
import qcutil
import os.path
import conf

class QCBoardVCFSTAT():
    def __init__(self, opt):
        self.opt = opt
        self.set_outfilenames()
        self.qcstat = self.init_qcstat()

    def set_outfilenames(self):
        if 'out' in self.opt.keys() and self.opt['out'].endswith('.json'):
            self.out_json = self.opt['out'][:-5] + '.stat.json'
        elif 'out' in self.opt.keys() and self.opt['out'] != '':
            self.out_json = self.opt['out'] + '.stat.json'
        else:
            self.out_json = self.opt['vcf'] + '.stat.json'

    def init_qcstat(self):
        stat = {}
        for k1 in list(conf.FIELDLIST_VCFQC.keys()):
            stat[k1] = {}
            for k2 in conf.FIELDLIST_VCFQC[k1]:
                stat[k1][k2] = 0
                stat[k1][k2+'_PASS'] = 0
                stat[k1][k2+'_NONPASS'] = 0
        return stat

    def get_annot_stat(self, info, annovar_headeridx, annotadd_headeridx):
        astat = {}
        for infofield in info.split(';'):
            if 'ANNOVAR=' in infofield:
                arr = infofield.replace('ANNOVAR=','').split('|')
                exonicfunc = arr[annovar_headeridx['ExonicFunc.ensGene']].strip()
                func = arr[annovar_headeridx['Func.ensGene']].strip()
                if exonicfunc != '':
                    func = exonicfunc
                astat['VARFUNC'] = func
                astat['REPEAT'] = arr[annovar_headeridx['all_repeats.b37']].strip()

            if 'ANNOTADD=' in infofield:
                arr = infofield.replace('ANNOTADD=','').split('|')
                try:
                    popaf = arr[annotadd_headeridx['gnomAD211.AF_popmax']].strip()
                    if popaf == '':
                        astat['POPAF'] = 0.0
                    else:
                        astat['POPAF'] = float(popaf)
                except IndexError:
                    astat['POPAF'] = 0.0

                # cpg = arr[annotadd_headeridx['CADD14.CpG']].strip()
                # print (cpg)
                
        return astat

    def get_headeridx(self, annotype, line, headeridx):
        if line[:len("##INFO=<ID="+annotype)] == "##INFO=<ID="+annotype:
            i = 0
            for field in line.split('Functional annotations:')[-1].replace(' ">','').replace("'","").strip().split(' | '):
                headeridx[field] = i
                i += 1
        return headeridx

    def load_vcf(self):
        stat = self.qcstat
        # print (stat)
        annovar_headeridx = {}
        annotadd_headeridx = {}
        ii = 0
        for line in qcutil.gzopen(self.opt['vcf']):
            if self.opt['vcf'].endswith('.gz'):
                line = line.decode('UTF-8')
            if line[:6] == "#CHROM":
                arr = line.split('\t')
                arr[-1] = arr[-1].strip()
                self.samplenames = arr[9:]
            elif line[0] == '#':
                annovar_headeridx = self.get_headeridx('ANNOVAR', line, annovar_headeridx)
                annotadd_headeridx = self.get_headeridx('ANNOTADD', line, annotadd_headeridx)
            elif line[0] != '#':
                ii += 1
                arr = line.split('\t')
                arr[-1] = arr[-1].strip()
                chrom = arr[0].replace('chr','').replace('MT','M')
                pos = arr[1]
                vqsr_filter = arr[conf.VCF_HEADCOL.index('FILTER')].strip()
                ref = arr[conf.VCF_HEADCOL.index('REF')].strip()
                alt = arr[conf.VCF_HEADCOL.index('ALT')].strip()

                info = arr[conf.VCF_HEADCOL.index('INFO')].strip()
                annotstat = self.get_annot_stat(info, annovar_headeridx, annotadd_headeridx)

                k = 9
                callinfo = arr[k].strip().split(':')
                gt = callinfo[0]
                
                homhet = qcutil.get_homhet(gt)
                (adj_ref, adj_alt, adj_gt) = qcutil.get_adjusted_refalt(ref,alt,gt)
                vartype = qcutil.get_vartype(adj_ref,adj_alt)

                if vqsr_filter == "PASS":
                    vt = '_PASS'
                else:
                    # vt = ''
                    vt = '_NONPASS'

                stat['VARTYPE'][vartype+vt] += 1

                # print (vartype, arr[:5])
                if vartype=='SNV':
                    titv = qcutil.get_titv(adj_ref,adj_alt)
                    subst = qcutil.get_substitution(adj_ref,adj_alt)
                    stat['TITV'][titv + vt] += 1
                    stat['HETHOM'][homhet + vt] += 1 
                    stat['SUBSTITUTION'][subst + vt] += 1
                    stat['CHROM'][chrom+ vt] += 1

                    if 'VARFUNC' in annotstat.keys():
                        stat = qcutil.check_key(stat, 'DNDS',{'DN':0,'DS':0,'DNDS':0.0})
                        if annotstat['VARFUNC']=='synonymous_SNV':
                            stat['DNDS']['DS'+vt] += 1
                        if annotstat['VARFUNC']=='nonsynonymous_SNV':
                            stat['DNDS']['DN'+vt] += 1

                    if 'POPAF' in annotstat.keys():
                        stat = qcutil.check_key(stat, 'RARE',{'RARE1P':0,'RARE01P':0,'SINGLETON':0,'NPOPAF':0,'RARE1P_PER':0.0,'RARE01P_PER':0.0,'SINGLETON_PER':0.0})
                        stat['RARE']['NPOPAF'+vt] += 1
                        if annotstat['POPAF']<0.01:
                            stat['RARE']['RARE1P'+vt] += 1
                        if annotstat['POPAF']<0.001:
                            stat['RARE']['RARE01P'+vt] += 1
                        if annotstat['POPAF']==0:
                            stat['RARE']['SINGLETON'+vt] += 1
                    
                    if 'REPEAT' in annotstat.keys():
                        stat = qcutil.check_key(stat, 'REPEAT',{'REPEAT':0,'NONREPEAT':0,'REPEAT_PER':0.0})
                        if annotstat['REPEAT']!="":
                            stat['REPEAT']['REPEAT'+vt] += 1
                        else:
                            stat['REPEAT']['NONREPEAT'+vt] += 1

                elif vartype=='INS' or vartype == 'DEL':
                    stat['HETHOM_INDEL'][homhet+vt] += 1
                    stat['CHROM_INDEL'][chrom+vt] += 1

                if ii % 500000 == 0:
                    print (ii, 'chr'+chrom, pos)
                    # break
                # break
        
        self.qcstat = stat
        self.cal_qcstat()
        # print (self.qcstat)

    def cal_qcstat(self):
        stat = self.qcstat
        print (stat)
        stat['QCBOARD']['NOSAMPLE'] = len(self.samplenames)
        print (stat['VARTYPE'])

        for k1 in stat.keys():
            for k2 in stat[k1].keys():
                if not '_PASS' in k2 and not '_NONPASS' in k2:
                    stat[k1][k2] = stat[k1][k2+'_PASS'] + stat[k1][k2+'_NONPASS']

        for vt in ['','_PASS','_NONPASS']:
            
            if stat['TITV']['TV'+vt] == 0:
                stat['TITV']['TITV'+vt] = 0
            else:
                stat['TITV']['TITV'+vt] = round(stat['TITV']['TI'+vt] / stat['TITV']['TV'+vt],4)

            for vartype in ['','_INDEL']:

                for chrom in list(stat['CHROM'+vartype].keys()):
                    if not 'TOTAL' in chrom:
                        if (vt == '' and not '_PASS' in chrom) or (vt == '_PASS' and '_PASS' in chrom):
                            # stat['CHROM'+vartype]['TOTAL'+vt] += stat['CHROM'+vartype][chrom]
                            pass

                if stat['HETHOM'+vartype]['HOMALT'+vt] == 0:
                    stat['HETHOM'+vartype]['HETHOM'+vt] = 0
                else:
                    stat['HETHOM'+vartype]['HETHOM'+vt] = round(stat['HETHOM'+vartype]['HET'+vt] / stat['HETHOM'+vartype]['HOMALT'+vt],4)
        
            stat['VARTYPE']['INDEL'+vt] = stat['VARTYPE']['INS'+vt] + stat['VARTYPE']['DEL'+vt]    
            stat['VARTYPE']['ALL'+vt] = stat['VARTYPE']['SNV'+vt] + stat['VARTYPE']['INDEL'+vt] + stat['VARTYPE']['MNV'+vt]

            sumsubst = 0
            for subst in list(stat['SUBSTITUTION'].keys()):
                if (vt == '' and not '_PASS' in subst) or (vt == '_PASS' and '_PASS' in subst):
                    sumsubst += stat['SUBSTITUTION'][subst]

            for subst in list(stat['SUBSTITUTION'].keys()):
                if (vt == '' and not '_PASS' in subst) or (vt == '_PASS' and '_PASS' in subst):
                    stat['SUBSTITUTION'][subst+'_PER'] = round(stat['SUBSTITUTION'][subst] / sumsubst * 100,2)

            if 'DNDS' in stat.keys():
                stat['DNDS']['DNDS'+vt] = round(stat['DNDS']['DN'+vt] / stat['DNDS']['DS'+vt],5)
            if 'RARE' in stat.keys():
                stat['RARE']['RARE1P_PER'+vt] = round(stat['RARE']['RARE1P'+vt] / stat['RARE']['NPOPAF'+vt]*100,2)
                stat['RARE']['RARE01P_PER'+vt] = round(stat['RARE']['RARE01P'+vt] / stat['RARE']['NPOPAF'+vt]*100,2)
                stat['RARE']['SINGLETON_PER'+vt] = round(stat['RARE']['SINGLETON'+vt] / stat['RARE']['NPOPAF'+vt]*100,2)
            if 'REPEAT' in stat.keys():
                stat['REPEAT']['REPEAT_PER'+vt] = round(stat['REPEAT']['REPEAT'+vt] / (stat['REPEAT']['REPEAT'+vt] + stat['REPEAT']['NONREPEAT'+vt])*100,2)

        # print (stat)
        self.qcstat = stat

    def save_stat(self):
        qcutil.savejson(self.out_json, self.qcstat)
        print ('Saved '+self.out_json)

    def run(self):
        self.load_vcf()
        self.save_stat()