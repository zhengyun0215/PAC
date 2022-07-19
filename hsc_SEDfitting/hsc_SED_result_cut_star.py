# from astropy import table
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.table import Table,Column,hstack
from astropy.io import fits


#########load data


# star = fits.getdata('/home/kunxu/red_disk/PCF/star_flag_sort.fits', 1)
sat2 = Table.read("./results_0p25.fits")
print("Load file successfully!!!!")



# sat2 = sat2[star['z_extendedness_value'] == 1]
print("The length of target:",len(sat2))


sat2 = sat2['bayes.param.restframe_Lnu(sdss.gp)','bayes.param.restframe_Lnu(sdss.gp)_err','bayes.param.restframe_Lnu(sdss.rp)','bayes.param.restframe_Lnu(sdss.rp)_err','bayes.param.restframe_Lnu(sdss.up)','bayes.param.restframe_Lnu(sdss.up)_err','bayes.stellar.m_star','bayes.stellar.m_star_err','best.param.restframe_sdss.gp-sdss.rp','best.param.restframe_sdss.up-sdss.rp','best.param.restframe_Lnu(sdss.gp)','best.param.restframe_Lnu(sdss.rp)','best.param.restframe_Lnu(sdss.up)','best.stellar.m_star']

Table(sat2).write("/home/yunzheng/PAC/code/hsc_SEDfitting/results_0p25_new.fits")