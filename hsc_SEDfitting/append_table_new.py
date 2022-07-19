from astropy.table import Table, Column, vstack
import os
import numpy as np
import multiprocessing as mp
from multiprocessing import Pool
import time

print("start")
print("Start ! :)")
start = time.time()

# os.chdir("/home/yunzheng/PAC/code/hsc_SEDfitting/0p25/")
def select(i):
    photo = Table.read("/home/yunzheng/PAC/code/hsc_SEDfitting/0p25/%s/out/results.fits"%i)
    print("finish:",i)
    return photo


core = 72 * 3
print("The number of core:",core)

pool = mp.Pool(processes=160)
res = pool.map(select, range(core))
photo_all = vstack(res)
photo_all.write("/home/yunzheng/PAC/code/hsc_SEDfitting/results_0p25_new.fits")
pool.close()
pool.join()

end = time.time()
print("The total time:")
print((end - start) / 60)
# -*- coding:utf-8 -*-
