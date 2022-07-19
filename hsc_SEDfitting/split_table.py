from astropy.io import fits
import numpy as np
import os
from multiprocessing import Pool

def split(z):
    j = str(z)
    core = 72*3
    os.chdir('./'+j.replace('.','p'))
    data = fits.getdata('./sed_red.fits', 1)
    total = 221208319
    each = total//(core-1)
    for i in range(core):
        os.system('mkdir ' + str(i))
        os.system('cp pci* ./' + str(i) +'/')
        if i != core - 1:
            cut = data[i*each:(i+1)*each]
        else:
            cut = data[i*each:]
        hdu = fits.BinTableHDU(data = cut)
        hdu.writeto('./'+str(i)+'/sed_red.fits')


# z = [0.35,0.45,0.55,0.65]
# z = [0.25,0.35,0.45]
z = [0.25]
with Pool(len(z)) as p:
    p.map(split,z)
    
print("OK~")
