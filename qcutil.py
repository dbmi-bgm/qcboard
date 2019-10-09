#!/usr/bin/env python

import os
import os.path
import gzip
from jinja2 import Template
from jinja2 import evalcontextfilter
from jinja2 import Environment, FileSystemLoader, select_autoescape
import conf
import json

# @evalcontextfilter
def percentage(n1, n2):
    if int(n2.replace(',','')) == 0:
        rst = '0'
    else:
        rst = comma(round(int(n1.replace(',',''))/int(n2.replace(',',''))*100,1))
    return rst
    

def add_comma_with_dict(d1):
    for k1 in list(d1.keys()):
        if type(d1[k1]) == type({}):
            for k2 in list(d1[k1].keys()):
                if type(d1[k1][k2]) == type(1) or type(d1[k1][k2]) == type(1.1):
                    d1[k1][k2] = comma(d1[k1][k2])
                else:
                    d1[k1][k2] = d1[k1][k2]
            print (k1, k2, d1[k1][k2])
        else:
            if type(d1[k1]) == type(1) or type(d1[k1]) == type(1.1):
                d1[k1] = comma(d1[k1])
            else:
                d1[k1] = d1[k1]
    return d1

def render(outfile,template_file,data):
    env = Environment(loader=FileSystemLoader(conf.DIRPATH))
    env.filters['percentage'] = percentage
    template = env.get_template(template_file)
    # template = Template(fileOpen(os.path.join(conf.DIRPATH,template_file)))
    fileSave(outfile, template.render(data), 'w')

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

def is_exist(fpath):
    return os.path.exists(fpath)

def gzopen(fname):
    if fname.endswith(".gz") or fname.endswith(".zip"):
        f1 = gzip.GzipFile(fname, "r")
    else:
        f1 = open(fname)
    return f1

def savejson(out, d):
    with open(out, 'w') as outfile:
        json.dump(d, outfile)

def loadjson(json_file):
    data = {}
    with open(json_file) as jf:
        data = json.load(jf)
    return data

def check_key(d, k1, init_value):
    try:
        tmp = d[k1]
    except KeyError:
        tmp_value = {}
        for k2 in init_value.keys():
            tmp_value[k2] = init_value[k2]
            tmp_value[k2+'_PASS'] = init_value[k2]
            tmp_value[k2+'_NONPASS'] = init_value[k2]
        d[k1] = tmp_value
    return d

def get_genotype(gt,ref, alt):
    geno = ""
    for g1 in gt.replace('|','/').split('/'):
        if g1 == '0':
            geno += ref
        elif g1 == '1':
            geno += alt
    return geno


def get_homhet(gt):
    homhet = ''
    tgt = gt.replace('|','/')
    aa = tgt.split('/')
    if tgt == '0/1':
        homhet = 'HET'
    elif tgt == '1/1':
        homhet = 'HOMALT'
    elif tgt == '0/0':
        homhet = 'HOMREF'
    elif tgt == './.':
        homhet = 'MISS'
    else:
        if aa[0] == '0':
            homhet = 'HET'
        else:
            homhet = 'HETALT'
    return homhet

##TODO:
def get_adjusted_refalt(ref,alt,gt):
    adj_ref = ref
    adj_alt = alt
    adj_gt = gt
    tgt = gt.replace('|','/')
    if gt == "1/2":   ##### For the single sample
        altaa = alt.split(',')
        adj_ref = altaa[0]
        adj_alt = altaa[1]
        adj_gt = "0/1"
    return (adj_ref, adj_alt, adj_gt)

def get_vartype(ref, alt):
    vartype = ''
    if ref == '*':
        vartype = 'INS'
    elif alt == '*':
        vartype = 'DEL'
    elif len(ref) == 1 and len(alt) == 1:
        vartype = 'SNV'
    elif len(ref)==1 and ',' in alt:
        altaa = alt.split(',')
        flag = True
        for j in range(len(altaa)):
            if len(altaa[j]) > 1:
                flag = False
        if flag:
            vartype = 'SNV'
        else:
            vartype = 'INS'
    elif len(ref) < len(alt) and not ',' in alt:
        vartype = 'INS'
    elif len(ref) > len(alt) and not ',' in alt:
        vartype = 'DEL'
    elif len(ref) == len(alt):
        vartype = 'MNV'
    else:
        vartype = 'MIXED'
    return vartype

TIMAP={}
TIMAP['A_G']='TI'
TIMAP['G_A']='TI'
TIMAP['C_T']='TI'
TIMAP['T_C']='TI'

def get_titv(ref,alt):
    global TIMAP
    try:
        titv = TIMAP[ref+'_'+alt]
    except KeyError:
        titv = 'TV'
    return titv

COMPLSUBSTITUTION = {}
COMPLSUBSTITUTION['C_A'] = 'C_A'
COMPLSUBSTITUTION['C_G'] = 'C_G'
COMPLSUBSTITUTION['C_T'] = 'C_T'
COMPLSUBSTITUTION['T_A'] = 'T_A'
COMPLSUBSTITUTION['T_C'] = 'T_C'
COMPLSUBSTITUTION['T_G'] = 'T_G'
COMPLSUBSTITUTION['G_T'] = 'C_A'
COMPLSUBSTITUTION['G_C'] = 'C_G'
COMPLSUBSTITUTION['G_A'] = 'C_T'
COMPLSUBSTITUTION['A_T'] = 'T_A'
COMPLSUBSTITUTION['A_G'] = 'T_C'
COMPLSUBSTITUTION['A_C'] = 'T_G'
def get_substitution(ref, alt):
    return (COMPLSUBSTITUTION[ref+'_'+alt])

