from multiprocessing import Pool
#from schwimmbad import MPIPool
import os
import numpy as np
import sys

def run(i):
    dic = '/home/yunzheng/PAC/code/hsc_SEDfitting/0p25/'+str(i)+'/'
    os.chdir(dic)
    os.system('pcigale run > run.log')

i = np.arange(72)+72+72
with Pool(72) as p:
    p.map(run,i)
    
print("25_3 is OK~")
