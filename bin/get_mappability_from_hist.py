#!/usr/bin/env python
# -*- coding: utf-8 -*-
#### get_mappability_hist.py
#### made by Daniel Minseok Kwon
#### 2019-10-08 06:11:48
#########################
import sys
import os
SVRNAME = os.uname()[1]
if "MBI" in SVRNAME.upper():
    sys_path="/Users/pcaso/bin/python_lib"
elif SVRNAME == "T7":
    sys_path="/ms1/bin/python_lib"
else:
    sys_path="/home/mk446/bin/python_lib"
sys.path.append(sys_path)

import file_util
import proc_util


def get_mappability_hist(histfile):
    flag = True
    cont = []
    for line in open(histfile):
        arr = line.split('\t')
        arr[-1] = arr[-1].strip()
        if arr[1] == "0" or (arr[0]=="MT" and flag):
            mappable_len = int(arr[3]) - int(arr[2])
            chrom = arr[0]
            # print (arr, mappable_len)
            if arr[0] == "MT" or arr[0] == "chrM":
                flag = False
                break

            cont.append("'"+chrom+"':" + str(mappable_len))
    
    print (','.join(cont))

if __name__ == "__main__":
    histfile = sys.argv[1]
    get_mappability_hist(histfile)
