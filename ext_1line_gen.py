# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 12:46:30 2022

@author: sasha
"""

import numpy as np
import pandas as pd
from scipy import stats
from matplotlib import pyplot as plt
import glob2 as glob
import seaborn
import os
import sys
import pathlib
from pathlib import Path

#%%
i1 = 0

file_list = glob.glob("/results/sensitivity/*.csv")

add_str = 'oneline_'


for x in file_list:
    i1 = i1+1
    print(i1)
    
    print(x)
    filef = pathlib.Path(add_str + x)
    a = np.zeros(26)
    
    if filef.exists()==0:
            
        # with open(x, 'rb') as f:
        #     try:  # catch OSError in case of a one line file 
        #         f.seek(-2, os.SEEK_END)
        #         while f.read(1) != b'\n':
        #             f.seek(-2, os.SEEK_CUR)
        #     except OSError:
        #         f.seek(0)
        #     last_line = f.readline().decode()
        
        with open(x) as f:
            for line in f:
                pass
            last_line = line
            
        ll=last_line.split(",")
            
        for i2 in range(0,np.size(a)):
            a[i2] = float(ll[i2])
                
                
            np.savetxt(add_str + x,a,delimiter=',')
            
    elif filef.exists() == 1:
        print('done')
        