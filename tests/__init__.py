#!/usr/bin/env python
# -*- coding: utf-8 -*-
#### __init__.py
#### made by Daniel Minseok Kwon
#### 2019-10-03 11:31:06
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


def __init__():
    


if __name__ == "__main__":
    __init__()
