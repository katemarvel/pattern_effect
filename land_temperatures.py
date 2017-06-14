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


f = cdms.open("DATA/TIMESERIES/cmip5.amip.ts.nc")

land_ts = MV.zeros(f["ts"].shape)+1.e20
models = eval(f["ts"].getAxis(0).models)
i=0
for model in models:
    print model
    land=cmip5.landfrac(model)
    fland = cdms.open(land)
    lf = fland("sftlf")
    ft = cdms.open(model)
    ts = ft("ts",time=('1979-1-1','2005-12-31'))
    bigland = np.repeat(lf.asma()[np.newaxis,:,:],ts.shape[0],axis=0)
    land_only_mask =cdutil.averager(MV.masked_where(bigland<100.,ts),axis='xy')
    try:
        land_ts[i] = land_only_mask
    except:
        print "problem with "+model
    i+=1
    fland.close()
    ft.close()
land_ts.id = "ts_land"
land_ts.name = "ts_land"
land_ts.setAxisList(f["ts"].getAxisList())
land_ts.long_name ="Land skin temperature"
fw = cdms.open("DATA/TIMESERIES/cmip5.amip.ts_land.nc","w")
fw.write(land_ts)
fw.close()
f.close()
