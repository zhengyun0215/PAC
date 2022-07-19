from astropy.io import fits
from astropy.table import Table, vstack, hstack
import numpy as np
import os
from mpi4py import MPI
# from mpi4py import MPI

print("ok!")

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

core = 72*3
z = ['0p35']
for i in range(core):
    print(i)
    t = Table.read('./%s/%d/out/results.fits' % (z[rank],i))
    if i==0:
        t_need = t
    else:
        t_need = vstack([t_need,t])
t_need.write('results_%s.fits' % z[rank],format='fits')
print("Finish!")