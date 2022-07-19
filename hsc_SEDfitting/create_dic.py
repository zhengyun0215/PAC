import os
import numpy as np
from astropy.io import fits
from multiprocessing import Pool


def table(i):
    j = str(i)
    a = fits.open('sed_red.fits')
    data = a[1].data
    new = data.columns
    new.del_col('redshift')
    redshift = np.zeros(len(data)) + i
    old = fits.ColDefs([fits.Column(name = 'redshift',format = 'D',array = redshift)])
    hdu = fits.BinTableHDU.from_columns(new+old)
   # a.close()
    hdu.writeto('./'+j.replace('.','p')+'/sed_red.fits')

# z = [0.55,0.65]
# z = [0.35,0.45]
z = [0.25]
# z = [0.025,0.075,0.125,0.175,0.225,0.275]
with Pool(len(z)) as p:
    p.map(table,z)
    
print("OK!")

