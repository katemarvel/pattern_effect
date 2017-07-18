import glob
import sys
import cdms2 as cdms
import numpy as np
import MV2 as MV
import difflib
import pickle
#import ENSO_years_piControl as en

global crunchy
import socket
if socket.gethostname().find("crunchy")>=0:
    crunchy = True
else:
    crunchy = False


import cdtime,cdutil,genutil
from eofs.cdms import Eof
from eofs.multivariate.cdms import MultivariateEof
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as patches
### Set classic Netcdf (ver 3)
cdms.setNetcdfShuffleFlag(0)
cdms.setNetcdfDeflateFlag(0)
cdms.setNetcdfDeflateLevelFlag(0)


import scipy.stats as stats
from scipy.interpolate import interp1d

if crunchy:
    sys.path.append("/work/marvel1/python-utils")
else:
    sys.path.append("~/Google Drive/python-utils")
from Plotting import *
import CMIP5_tools as cmip5


def last_ten_minus_first_ten(X):
    axis = X.getAxisIds().index('time')
    for_axes = MV.average(X,axis=axis)
    firstten = [slice(None)] * len(X.shape)
    firstten[axis]=slice(0,10)
    average_firstten = MV.average(X.asma()[firstten],axis=axis)

    lastten = [slice(None)] * len(X.shape)
    nt = X.shape[axis]
    
    lastten[axis]=slice(nt-10,nt)
    average_lastten = MV.average(X.asma()[lastten],axis=axis)

    diff_ten = average_lastten-average_firstten
    diff_ten.setAxisList(for_axes.getAxisList())
    return diff_ten


def low_cloud_diff(fname):
    f=cdms.open(fname)
    clisccp = f("clisccp")
    axes = clisccp.getAxisIds()
    tau_ax = axes.index("tau")
    clisccp_all_od = MV.sum(clisccp,axis=tau_ax)
    plev_ax = clisccp_all_od.getAxisIds().index("plev")
    low = MV.sum(clisccp_all_od(level=(1000*100,440.*100)),axis=plev_ax)
    cdutil.setTimeBoundsMonthly(low)
    return last_ten_minus_first_ten(cdutil.YEAR(low))
    
