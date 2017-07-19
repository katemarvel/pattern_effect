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


def low_cloud_diff(clisccp):
    #f=cdms.open(fname)
    #clisccp = f("clisccp")
    axes = clisccp.getAxisIds()
    tau_ax = axes.index("tau")
    clisccp_all_od = MV.sum(clisccp,axis=tau_ax)
    plev_ax = clisccp_all_od.getAxisIds().index("plev")
    low = MV.sum(clisccp_all_od(level=(1000*100,440.*100)),axis=plev_ax)(time=('1979-1-1','2005-12-31'))
    cdutil.setTimeBoundsMonthly(low)
    the_diff= last_ten_minus_first_ten(cdutil.YEAR(low))

    fobs = cdms.open("/work/marvel1/CLOUD_SEASONS/cloud-seasons/CLOUD_OBS/clt_ISCCP_corrected_198301-200912.nc")
    the_grid = fobs["clt"].getGrid()
    the_diff_regrid = the_diff.regrid(the_grid,regridTool='regrid2')
    fobs.close()
    return the_diff_regrid
    
def ensemble_AMIP_LCC():
    direc = "/work/cmip5/amip/atm/mo/clisccp/"
    variable="clisccp"
    AMIP_LCC = cmip5.get_ensemble(direc,variable,func=low_cloud_diff)
    AMIP_LCC.id = "lcc"
    f = cdms.open("AMIP_LCC.nc","w")
    f.write(AMIP_LCC)
    f.close()
    return AMIP_LCC
def ensemble_HISTORICAL_LCC():
    direc = "/work/cmip5/historical/atm/mo/clisccp/"
    variable="clisccp"
    HISTORICAL_LCC = cmip5.get_ensemble(direc,variable,func=low_cloud_diff)
    HISTORICAL_LCC.id = "lcc"
    f = cdms.open("HISTORICAL_LCC.nc","w")
    f.write(HISTORICAL_LCC)
    f.close()
    return HISTORICAL_LCC

def rsut_diff(X):
    
    Xt=X(time=('1979-1-1','2005-12-31'))
    cdutil.setTimeBoundsMonthly(X)
    the_diff= last_ten_minus_first_ten(cdutil.YEAR(X))

    fobs = cdms.open("/work/marvel1/CLOUD_SEASONS/cloud-seasons/CLOUD_OBS/clt_ISCCP_corrected_198301-200912.nc")
    the_grid = fobs["clt"].getGrid()
    the_diff_regrid = the_diff.regrid(the_grid,regridTool='regrid2')
    fobs.close()
    return the_diff_regrid
    
def ensemble_AMIP_RSUT():
    direc = "/work/cmip5/amip/atm/mo/rsut/"
    variable="rsut"
    AMIP_RSUT = cmip5.get_ensemble(direc,variable,func=rsut_diff)
    AMIP_RSUT.id = "rsut"
    f = cdms.open("AMIP_RSUT.nc","w")
    f.write(AMIP_RSUT)
    f.close()
    return AMIP_RSUT

def ensemble_AMIP_RSUTCS():
    direc = "/work/cmip5/amip/atm/mo/rsutcs/"
    variable="rsutcs"
    AMIP_RSUTCS = cmip5.get_ensemble(direc,variable,func=rsut_diff)
    AMIP_RSUTCS.id = "rsutcs"
    f = cdms.open("AMIP_RSUTCS.nc","w")
    f.write(AMIP_RSUTCS)
    f.close()
    return AMIP_RSUTCS

def ensemble_HISTORICAL_RSUT():
    direc = "/work/cmip5/historical/atm/mo/rsut/"
    variable="rsut"
    HISTORICAL_RSUT = cmip5.get_ensemble(direc,variable,func=rsut_diff)
    HISTORICAL_RSUT.id = "rsut"
    f = cdms.open("HISTORICAL_RSUT.nc","w")
    f.write(HISTORICAL_RSUT)
    f.close()
    return HISTORICAL_RSUT

def ensemble_HISTORICAL_RSUTCS():
    direc = "/work/cmip5/historical/atm/mo/rsutcs/"
    variable="rsutcs"
    HISTORICAL_RSUTCS = cmip5.get_ensemble(direc,variable,func=rsut_diff)
    HISTORICAL_RSUTCS.id = "rsutcs"
    f = cdms.open("HISTORICAL_RSUTCS.nc","w")
    f.write(HISTORICAL_RSUTCS)
    f.close()
    return HISTORICAL_RSUTCS

#HISTORICAL_RSUT = ensemble_HISTORICAL_RSUT()
#HISTORICAL_RSUTCS = ensemble_HISTORICAL_RSUTCS()
def historical_SWCRE(HISTORICAL_RSUT,HISTORICAL_RSUTCS):
      allsmod_hist=np.array([x.split("/")[-1].split("rsut")[0] for x in eval(HISTORICAL_RSUT.getAxis(0).models)])
      csmod_hist=np.array([x.split("/")[-1].split("rsutcs")[0] for x in eval(HISTORICAL_RSUTCS.getAxis(0).models)])
      common_models = np.intersect1d(allsmod_hist,csmod_hist)
      nmod=len(common_models)
      HISTORICAL_SWCRE = MV.zeros((nmod,)+HISTORICAL_RSUTCS.shape[1:])
      counter=0
      for model in common_models:
            i = allsmod_hist.tolist().index(model)
            j = csmod_hist.tolist().index(model)
            HISTORICAL_SWCRE[counter]=HISTORICAL_RSUT[i] - HISTORICAL_RSUTCS[j]
            counter+=1
      HISTORICAL_SWCRE.setAxis(1,HISTORICAL_RSUT.getLatitude())
      HISTORICAL_SWCRE.setAxis(2,HISTORICAL_RSUT.getLongitude())
      modax = cmip5.make_model_axis(np.array(common_models).tolist())
      HISTORICAL_SWCRE.setAxis(0,modax)
      HISTORICAL_SWCRE.id="SWCRE"
      f = cdms.open("HISTORICAL_SWCRE.nc","w")
      f.write(HISTORICAL_SWCRE)
      f.close()
      return HISTORICAL_SWCRE
      
      
      

