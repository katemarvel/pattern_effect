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
import cmip5_ECS as c


# historical = c.ECS("historical")
# fobs = cdms.open("sst.mnmean.v4.nc")
# sst_obs = fobs("sst",time=("1979-1-1","2005-12-31"))
# obs_grid = sst_obs.getGrid()
# sst_patt = MV.average(sst_obs(time=('1996-1-1','2005-12-31')),axis=0)-MV.average(sst_obs(time=('1979-1-1','1988-12-31')),axis=0)
# hmodels = eval(historical.dT.getAxis(0).models)
# covariances = np.zeros(len(hmodels))+1.e20
# path = "/work/cmip5/historical/atm/mo/tas/"
# i=0
# for model in hmodels:
#     print model
#     try:
#         f = cdms.open(cmip5.get_latest_version(glob.glob(path+model+"*")))
#         sst_mod = f("tas",time=("1979-1-1","2005-12-31"))
#         sst_mod_patt = MV.average(sst_mod(time=('1996-1-1','2005-12-31')),axis=0)-MV.average(sst_mod(time=('1979-1-1','1988-12-31')),axis=0)
#         sst_mod_patt_regrid = sst_mod_patt.regrid(obs_grid,regridTool='regrid2')
#         sst_mod_patt_regrid_mask = MV.masked_where(sst_patt.mask,sst_mod_patt_regrid)
#         covariances[i] = genutil.statistics.covariance(sst_patt,sst_mod_patt_regrid_mask,axis='xy')
#         f.close()
#     except:
#         continue
#     i+=1

    


def get_warming_pattern(model):
    experiment = amodels[0].split(".")[2]
    path = "/work/cmip5/"+experiment+"/atm/mo/tas/"
    fobs = cdms.open("sst.mnmean.v4.nc")
    sst_obs = fobs("sst",time=("1979-1-1","2005-12-31"))
    obs_grid = sst_obs.getGrid()
    fobs.close()
    f = cdms.open(cmip5.get_latest_version(glob.glob(path+model+"*")))
    sst_mod = f("tas",time=("1979-1-1","2005-12-31"))
    sst_mod_patt = MV.average(sst_mod(time=('1996-1-1','2005-12-31')),axis=0)-MV.average(sst_mod(time=('1979-1-1','1988-12-31')),axis=0)
    sst_mod_patt_regrid = sst_mod_patt.regrid(obs_grid,regridTool='regrid2')
    
    return sst_mod_patt_regrid

def get_regridded_pattern(model,variable):
    fobs = cdms.open("sst.mnmean.v4.nc")
    sst_obs = fobs("sst",time=("1979-1-1","2005-12-31"))
    obs_grid = sst_obs.getGrid()
    fobs.close()
    
    experiment = model.split(".")[2]
    path = "/work/cmip5/"+experiment+"/atm/mo/"+variable+"/"
   
    f = cdms.open(cmip5.get_latest_version(glob.glob(path+model+"*")))
    sst_mod = f(variable,time=("1979-1-1","2005-12-31"))
    sst_mod_patt = MV.average(sst_mod(time=('1996-1-1','2005-12-31')),axis=0)-MV.average(sst_mod(time=('1979-1-1','1988-12-31')),axis=0)
    sst_mod_patt_regrid = sst_mod_patt.regrid(obs_grid,regridTool='regrid2')
    return sst_mod_patt_regrid


def abrupt_patterns(X):
   
    fobs = cdms.open("sst.mnmean.v4.nc")
    sst_obs = fobs("sst",time=("1979-1-1","2005-12-31"))
    obs_grid = sst_obs.getGrid()
    fobs.close()

    last_ten_years = MV.average(X[130*12:140*12],axis=0)-MV.average(X[0*12:10*12],axis=0)
    lty_r = last_ten_years.regrid(obs_grid,regridTool='regrid2')
    return lty_r

    
def get_regridded_zonal(model,variable):
    fobs = cdms.open("sst.mnmean.v4.nc")
    sst_obs = fobs("sst",time=("1979-1-1","2005-12-31"))
    obs_grid = sst_obs.getGrid()
    fobs.close()
    
    experiment = model.split(".")[2]
    path = "/work/cmip5/"+experiment+"/atm/mo/"+variable+"/"
   
    f = cdms.open(cmip5.get_latest_version(glob.glob(path+model+"*")))
    sst_mod = f(variable,time=("1979-1-1","2005-12-31"))
    
    sst_mod_regrid = sst_mod_patt.regrid(obs_grid,regridTool='regrid2')
    return sst_mod_regrid

if __name__ == "__main__":
    fobs = cdms.open("sst.mnmean.v4.nc")
    sst_obs = fobs("sst",time=("1979-1-1","2005-12-31"))
    obs_grid = sst_obs.getGrid()
    fobs.close()


    #TOA = {"rsdt":"TOA Incident Shortwave Radiation",\
            "rsut": "TOA Outgoing Shortwave Radiation",\
            "rlut": "TOA Outgoing Longwave Radiation"}

    amip = c.ECS("amip")
    amodels = eval(amip.dT.getAxis(0).models)
    La = len(amodels)

    historical = c.ECS("historical")
    hmodels = eval(historical.dT.getAxis(0).models)
    Lh = len(hmodels)
    #variables = TOA.keys()+["tas"]
    variables = ["ts"]
    for variable in variables:
        A = MV.zeros((La,)+obs_grid.shape)+1.e20
        for i in range(La):
            try:
                A[i] = get_regridded_pattern(amodels[i],variable)
            except:
                print "BAD: "+amodels[i]

        modax = cmip5.make_model_axis(amodels)
        A = MV.masked_where(A>1.e10,A)
        A.id=variable
        A.setAxis(0,modax)
        A.setAxis(1,obs_grid.getLatitude())
        A.setAxis(2,obs_grid.getLongitude())
        fw = cdms.open("PATTERN_MAPS/amip."+variable+".nc","w")
        fw.write(A)
        fw.close()



        H = MV.zeros((Lh,)+obs_grid.shape)+1.e20
        for i in range(Lh):
            try:
                H[i] = get_regridded_pattern(hmodels[i],variable)
            except:
                print "BAD: "+hmodels[i]

        H = MV.masked_where(H>1.e10,H)
        H.id=variable
        H.name = variable
        modax = cmip5.make_model_axis(hmodels)
        H.setAxis(0,modax)
        H.setAxis(1,obs_grid.getLatitude())
        H.setAxis(2,obs_grid.getLongitude())
        fw = cdms.open("PATTERN_MAPS/historical."+variable+".nc","w")
        fw.write(H)
        fw.close()
