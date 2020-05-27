#!/usr/bin/env python
# -*- coding: utf-8 -*-
#### test_vcfqc.py
#### made by Daniel Minseok Kwon
#### 2020-05-26 17:28:46
#########################
import sys
import os
sys.path.append('..')
import qcutil
import qcboard

BASEPATH = os.path.abspath('..')
QCBOARD = os.path.join(BASEPATH, 'qcboard.py')

def test_family_vcfqc():
    global QCBOARD

    for n1 in [10, 100, 1000, 10000]:
        for r1 in range(1,6):
            vcf = "data/test_trio_"+str(n1)+"_"+str(r1)+".vcf.gz"
            out = "out/test_trio_"+str(n1)+"_"+str(r1)+".vcf.gz" + '.vcfqc'
            # out = vcf + '.vcfqc'
            print('-vcf:', vcf)
            
            opt = {} 
            opt['subcommand'] =  'vcfqc'
            opt['vcf']  = [vcf]
            opt['out'] = out
            opt['temp'] = 'qcboard_vcfqc.html'
            opt['relation'] = 'data/test_trio.relation.json'
            opt['silence'] = False
            opt['debug'] = False
            
            print(qcutil.opt2cmd(opt, QCBOARD))
            qcboard.run_with_option(opt)

            this_out = qcutil.loadjson(out + '.famstat.json')
            prev_out = qcutil.loadjson(vcf + '.vcfqc' + '.famstat.json')
            assert this_out == prev_out
        

if __name__ == "__main__":
    test_family_vcfqc()
