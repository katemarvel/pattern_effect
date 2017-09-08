import glob
import sys
import cdms2 as cdms
import numpy as np
import MV2 as MV
import difflib
import pickle
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap
import string
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
import Plotting
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

def last_thirty_minus_first_thirty(X):
    axis = X.getAxisIds().index('time')
    for_axes = MV.average(X,axis=axis)
    first30 = [slice(None)] * len(X.shape)
    first30[axis]=slice(0,30)
    average_first30 = MV.average(X.asma()[first30],axis=axis)

    last30 = [slice(None)] * len(X.shape)
    nt = X.shape[axis]
    
    last30[axis]=slice(nt-30,nt)
    average_last30 = MV.average(X.asma()[last30],axis=axis)

    diff_30 = average_last30-average_first30
    diff_30.setAxisList(for_axes.getAxisList())
    return diff_30
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
    if crunchy:
        fobs = cdms.open("/work/marvel1/CLOUD_SEASONS/cloud-seasons/CLOUD_OBS/clt_ISCCP_corrected_198301-200912.nc")
    else:
        fobs = cdms.open("/Users/kmarvel/Google Drive/CLOUD_SEASONS/cloud-seasons/CLOUD_OBS/clt_ISCCP_corrected_198301-200912.nc")
    the_grid = fobs["clt"].getGrid()
    the_diff_regrid = the_diff.regrid(the_grid,regridTool='regrid2')
    fobs.close()
    return the_diff_regrid


def low_cloud_diff(clisccp,thresh=680):
    #f=cdms.open(fname)
    #clisccp = f("clisccp")
    print "THRESHOLD = "+str(thresh)
    axes = clisccp.getAxisIds()
    tau_ax = axes.index("tau")
    clisccp_all_od = MV.sum(clisccp,axis=tau_ax)
    plev_ax = clisccp_all_od.getAxisIds().index("plev")
    low = MV.sum(clisccp_all_od(level=(1000*100,thresh*100.)),axis=plev_ax)(time=('1979-1-1','2005-12-31'))
    cdutil.setTimeBoundsMonthly(low)
    the_diff= last_ten_minus_first_ten(cdutil.YEAR(low))
    if crunchy:
        fobs = cdms.open("/work/marvel1/CLOUD_SEASONS/cloud-seasons/CLOUD_OBS/clt_ISCCP_corrected_198301-200912.nc")
    else:
        fobs = cdms.open("/Users/kmarvel/Google Drive/CLOUD_SEASONS/cloud-seasons/CLOUD_OBS/clt_ISCCP_corrected_198301-200912.nc")
    the_grid = fobs["clt"].getGrid()
    the_diff_regrid = the_diff.regrid(the_grid,regridTool='regrid2')
    fobs.close()
    return the_diff_regrid

def low_cloud_abrupt_diff(clisccp,thresh=680):
    #f=cdms.open(fname)
    #clisccp = f("clisccp")
    print "THRESHOLD = "+str(thresh)
    axes = clisccp.getAxisIds()
    tau_ax = axes.index("tau")
    clisccp_all_od = MV.sum(clisccp,axis=tau_ax)
    plev_ax = clisccp_all_od.getAxisIds().index("plev")
    low = MV.sum(clisccp_all_od(level=(1000*100,thresh*100)),axis=plev_ax)#(time=('1979-1-1','2005-12-31'))
    cdutil.setTimeBoundsMonthly(low)
    the_diff= MV.average(low[130*12:140*12],axis=0)-MV.average(low[0*12:10*12],axis=0)
    if crunchy:
        fobs = cdms.open("/work/marvel1/CLOUD_SEASONS/cloud-seasons/CLOUD_OBS/clt_ISCCP_corrected_198301-200912.nc")
    else:
        fobs = cdms.open("/Users/kmarvel/Google Drive/CLOUD_SEASONS/cloud-seasons/CLOUD_OBS/clt_ISCCP_corrected_198301-200912.nc")
    the_grid = fobs["clt"].getGrid()
    the_diff_regrid = the_diff.regrid(the_grid,regridTool='regrid2')
    fobs.close()
    return the_diff_regrid

def ensemble_ABRUPT_LCC(thresh):
    direc = "/work/cmip5/abrupt4xCO2/atm/mo/clisccp/"
    variable="clisccp"
    ABRUPT_LCC = cmip5.get_ensemble(direc,variable,func=low_cloud_abrupt_diff,search_string = "*r1*", thresh=thresh)
    ABRUPT_LCC.id = "lcc"
    f = cdms.open("ABRUPT_LCC.nc","w")
    f.write(ABRUPT_LCC)
    f.close()
    return ABRUPT_LCC

def ensemble_AMIP_LCC(thresh):
    direc = "/work/cmip5/amip/atm/mo/clisccp/"
    variable="clisccp"
    AMIP_LCC = cmip5.get_ensemble(direc,variable,func=low_cloud_diff,thresh=thresh)
    AMIP_LCC.id = "lcc"
    f = cdms.open("AMIP_LCC_"+str(thresh)+".nc","w")
    f.write(AMIP_LCC)
    f.close()
    return AMIP_LCC
def ensemble_HISTORICAL_LCC(thresh):
    direc = "/work/cmip5/historical/atm/mo/clisccp/"
    variable="clisccp"
    HISTORICAL_LCC = cmip5.get_ensemble(direc,variable,func=low_cloud_diff, thresh=thresh)
    HISTORICAL_LCC.id = "lcc"
    f = cdms.open("HISTORICAL_LCC_"+str(thresh)+".nc","w")
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

def rsut_diff_abrupt(X):
    
    Xt=X[0:12*140]
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

def ensemble_ABRUPT_RSUT():
    direc = "/work/cmip5/abrupt4xCO2/atm/mo/rsut/"
    variable="rsut"
    ABRUPT_RSUT = cmip5.get_ensemble(direc,variable,func=rsut_diff_abrupt,search_string="*r1i1*")
    ABRUPT_RSUT.id = "rsut"
    f = cdms.open("ABRUPT_RSUT.nc","w")
    f.write(ABRUPT_RSUT)
    f.close()
    return ABRUPT_RSUT

def ensemble_ABRUPT_RSUTCS():
    direc = "/work/cmip5/abrupt4xCO2/atm/mo/rsutcs/"
    variable="rsutcs"
    ABRUPT_RSUTCS = cmip5.get_ensemble(direc,variable,func=rsut_diff_abrupt,search_string="*r1i1*")
    ABRUPT_RSUTCS.id = "rsutcs"
    f = cdms.open("ABRUPT_RSUTCS.nc","w")
    f.write(ABRUPT_RSUTCS)
    f.close()
    return ABRUPT_RSUTCS

#HISTORICAL_RSUT = ensemble_HISTORICAL_RSUT()
#HISTORICAL_RSUTCS = ensemble_HISTORICAL_RSUTCS()

def abrupt_SWCRE(ABRUPT_RSUT,ABRUPT_RSUTCS):
      allsmod_hist=np.array([x.split("/")[-1].split("rsut")[0] for x in eval(ABRUPT_RSUT.getAxis(0).models)])
      csmod_hist=np.array([x.split("/")[-1].split("rsutcs")[0] for x in eval(ABRUPT_RSUTCS.getAxis(0).models)])
      common_models = np.intersect1d(allsmod_hist,csmod_hist)
      nmod=len(common_models)
      ABRUPT_SWCRE = MV.zeros((nmod,)+ABRUPT_RSUTCS.shape[1:])
      counter=0
      for model in common_models:
            i = allsmod_hist.tolist().index(model)
            j = csmod_hist.tolist().index(model)
            ABRUPT_SWCRE[counter]=ABRUPT_RSUT[i] - ABRUPT_RSUTCS[j]
            counter+=1
      ABRUPT_SWCRE.setAxis(1,ABRUPT_RSUT.getLatitude())
      ABRUPT_SWCRE.setAxis(2,ABRUPT_RSUT.getLongitude())
      modax = cmip5.make_model_axis(np.array(common_models).tolist())
      ABRUPT_SWCRE.setAxis(0,modax)
      ABRUPT_SWCRE.id="SWCRE"
      f = cdms.open("ABRUPT_SWCRE.nc","w")
      f.write(ABRUPT_SWCRE)
      f.close()
      return ABRUPT_SWCRE

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
      
      
#8/6: trying to figure out why the code in cmip5_ECS and averaging over the map yield different results
import cmip5_ECS as c
def compare_timeseries_and_map():
    sw_time = c.SWCRE("amip")
    cdutil.setTimeBoundsMonthly(sw_time)
    from_timeseries = last_ten_minus_first_ten(cdutil.YEAR(sw_time))
    f = cdms.open("MAPS/AMIP_SWCRE.nc")
    sw_map = f("SWCRE")
    from_map = cdutil.averager(sw_map,axis='xy')
    timeseries_models = [x.split(".")[1]+"."+x.split(".")[3] for x in eval(sw_time.getAxis(0).models)]
    map_models = [x.split(".")[1]+"."+x.split(".")[3] for x in eval(sw_map.getAxis(0).models)]
    models_in_common = np.intersect1d(map_models,timeseries_models)
    L = len(models_in_common)
    MAP=MV.zeros(L)
    TS = MV.zeros(L)
    counter=0
    for model in models_in_common:
        time_i = timeseries_models.index(model)
        TS[counter]=from_timeseries[time_i]
        map_i = map_models.index(model)
        MAP[counter]=from_map[map_i]
        counter+=1
    modax = make_model_axis(["cmip5."+x for x in models_in_common.tolist()])
    MAP.setAxis(0,modax)
    MAP.id="swcre_map"
    TS.setAxis(0,modax)
    TS.id = "swcre_timeseries"
    return TS,MAP


def add_giss_to_amip_lcc():
  
    f = cdms.open("MAPS/AMIP_LCC.nc")
    amip_lcc = f("lcc")
    f.close()
    nmodels = amip_lcc.shape[0]
    amip_new = MV.zeros((nmodels+2,)+amip_lcc.shape[1:])
    amip_new[2:]=amip_lcc
    models = eval(amip_lcc.getAxis(0).models)

    giss3model = '/Users/kmarvel/Google Drive/PATTERN_EFFECT/cmip5.GISS-E2-R.amip.r1i1p3.mo.atm.cfMon.clisccp.ver-1.latestX.nc'
    giss1model = '/Users/kmarvel/Google Drive/PATTERN_EFFECT/cmip5.GISS-E2-R.amip.r1i1p1.mo.atm.cfMon.clisccp.ver-1.latestX.nc'
    newmodels = [giss1model,giss3model]+models
    
    fg3 = cdms.open(giss3model)
    clisccp3 = fg3("clisccp")
    gissp3lcc=low_cloud_diff(clisccp3)
    fg3.close()

    fg1 = cdms.open(giss1model)
    clisccp1 = fg1("clisccp")
    gissp1lcc=low_cloud_diff(clisccp1)
    fg1.close()

    amip_new[0]=gissp1lcc
    amip_new[1]=gissp3lcc

    modax = cmip5.make_model_axis(newmodels)

    amip_new.setAxis(0,modax)
    for i in range(len(amip_lcc.shape))[1:]:
        amip_new.setAxis(i,amip_lcc.getAxis(i))
    amip_new.id='lcc'
    return amip_new

    

      
    

#8/7: compare low cloud cover and apparent ECS
#in Tropics, SO, etc.
def AMIP_LCC():
    f = cdms.open("MAPS/AMIP_LCC.nc")
    amip_lcc = f("lcc") 
    amip_estimator = c.ECS("amip")
    amip_estimator.estimate_ECS(typ="annual regression") #Still need to explore why this is systematically lower than decadal
    amip_lcc_ensav = cmip5.ensemble2multimodel(amip_lcc)
    amip_ecs_ensav = cmip5.ensemble2multimodel(amip_estimator.ECS)
    amip_ecs_models = eval(amip_ecs_ensav.getAxis(0).models)
    amip_lcc_models = eval(amip_lcc_ensav.getAxis(0).models)
    models_in_common = np.intersect1d(amip_ecs_models,amip_lcc_models)
    TROPICAL_LCC = cdutil.averager(amip_lcc_ensav(latitude=(-10,10)),axis='xy')
    SO_LCC = cdutil.averager(amip_lcc_ensav(latitude=(-90,-50)),axis='xy')
    TOTAL_LCC = cdutil.averager(amip_lcc_ensav,axis='xy')
    TROPICS={}
    SO={}
    TOTAL = {}
    for model in models_in_common:
        ecs_i = amip_ecs_models.index(model)
        lcc_i = amip_lcc_models.index(model)
        TROPICS[model]=[amip_ecs_ensav[ecs_i],TROPICAL_LCC[lcc_i]]

        SO[model]=[amip_ecs_ensav[ecs_i],SO_LCC[lcc_i]]
        TOTAL[model]=[amip_ecs_ensav[ecs_i],TOTAL_LCC[lcc_i]]

def HISTORICAL_LCC():
    f = cdms.open("MAPS/HISTORICAL_LCC.nc")
    historical_lcc = f("lcc") 
    historical_estimator = c.ECS("historical")
    historical_estimator.estimate_ECS(typ="annual regression") #Still need to explore why this is systematically lower than decadal
    historical_lcc_ensav = cmip5.ensemble2multimodel(historical_lcc)
    historical_ecs_ensav = cmip5.ensemble2multimodel(historical_estimator.ECS)
    historical_ecs_models = eval(historical_ecs_ensav.getAxis(0).models)
    historical_lcc_models = eval(historical_lcc_ensav.getAxis(0).models)
    models_in_common = np.intersect1d(historical_ecs_models,historical_lcc_models)
    historical_TROPICAL_LCC = cdutil.averager(historical_lcc_ensav(latitude=(-10,10)),axis='xy')
    historical_SO_LCC = cdutil.averager(historical_lcc_ensav(latitude=(-90,-50)),axis='xy')
    historical_TROPICS={}
    historical_SO={}
    for model in models_in_common:
        ecs_i = historical_ecs_models.index(model)
        lcc_i = historical_lcc_models.index(model)
        historical_TROPICS[model]=[historical_ecs_ensav[ecs_i],historical_TROPICAL_LCC[lcc_i]]

        historical_SO[model]=[historical_ecs_ensav[ecs_i],historical_SO_LCC[lcc_i]]

        
def models_in_common(X,Y):
    """ X and Y should be ensemble averages """
    Xmodels = eval(X.getAxis(0).models)
    Ymodels = eval(Y.getAxis(0).models)
    models_in_common = np.intersect1d(Xmodels,Ymodels)
    L = len(models_in_common)
    Xnew = MV.zeros((L,)+X.shape[1:])
    Ynew = MV.zeros((L,)+Y.shape[1:])
    counter = 0
    for model in models_in_common:
        xi = Xmodels.index(model)
        Xnew[counter]=X[xi]

        yi = Ymodels.index(model)
        Ynew[counter]=Y[yi]
        counter+=1
    newmodax = cmip5.make_model_axis(models_in_common.tolist())
    Xnew.setAxis(0,newmodax)
    Xnew.id = X.id
    for i in range(len(X.shape))[1:]:
        Xnew.setAxis(i,X.getAxis(i))
    Ynew.setAxis(0,newmodax)
    Ynew.id = Y.id
    for i in range(len(Y.shape))[1:]:
        Ynew.setAxis(i,Y.getAxis(i))
    return Xnew,Ynew
    

def tropical_marine_lcc_klein(X):
    #Stratus regions

    Peru = cdutil.region.domain(latitude=(-20,-10),longitude=(-90,-80))
    Namibia = cdutil.region.domain(latitude=(-20,-10),longitude=(0,10))
    California = cdutil.region.domain(latitude=(20,30),longitude=(-130,-120))
    Australia = cdutil.region.domain(latitude=(-35,-15),longitude=(95,105))
    Canaries = cdutil.region.domain(latitude=(15,25),longitude=(-35,-25))
    return cdutil.averager(X(Peru),axis='xy')+cdutil.averager(X(Namibia),axis='xy')+cdutil.averager(X(California),axis='xy')+cdutil.averager(X(Australia),axis='xy')+cdutil.averager(X(Canaries),axis='xy')


def tropical_marine_lcc_qu(X):
    #Stratus regions as defined by Qu et al 2015

    Peru = cdutil.region.domain(latitude=(-30,-10),longitude=(-110,-70))
    Namibia = cdutil.region.domain(latitude=(-30,-10),longitude=(-25,15))
    California = cdutil.region.domain(latitude=(15,35),longitude=(-155,-115))
    Australia = cdutil.region.domain(latitude=(-35,-15),longitude=(75,115))
    Canaries = cdutil.region.domain(latitude=(10,30),longitude=(-50,-10))
    
    return cdutil.averager(X(Peru),axis='xy')+cdutil.averager(X(Namibia),axis='xy')+cdutil.averager(X(California),axis='xy')+cdutil.averager(X(Australia),axis='xy')+cdutil.averager(X(Canaries),axis='xy')


def plot_tropical_lcc(X,projection='moll'):
    if X.id == "SWCRE":
        v=8
    else:
        v=2
    Peru = cdutil.region.domain(latitude=(-30,-10),longitude=(-110,-70))
    Namibia = cdutil.region.domain(latitude=(-30,-10),longitude=(-25,15))
    California = cdutil.region.domain(latitude=(15,35),longitude=(-155,-115))
    Australia = cdutil.region.domain(latitude=(-35,-15),longitude=(75,115))
    Canaries = cdutil.region.domain(latitude=(10,30),longitude=(-50,-10))
    lon = X.getLongitude().getBounds()[:,0]
    lat = X.getLatitude().getBounds()[:,0]
    
    m = Basemap(lon_0=0,projection=projection)
    #x,y=m(*np.meshgrid(lon,lat))
   # m.pcolormesh(x,y,X,vmin=2,vmax=2,alpha=.3)
    for region in [Peru,Namibia,California,Australia,Canaries]:
        Xr = X(region)
        lon = Xr.getLongitude().getBounds()[:,0]
        lat = Xr.getLatitude().getBounds()[:,0]
        x,y=m(*np.meshgrid(lon,lat))
    #if vmin is None:
        m.pcolormesh(x,y,Xr,vmin=-v,vmax=v)
    m.drawcoastlines()


def match_ensemble_members(X,Y):
    xmodels = cmip5.models(X)
    ymodels = cmip5.models(Y)
    xrips = ["cmip5."+thing.split(".")[1]+"."+thing.split(".")[2]+"."+thing.split(".")[3]+"." for thing in xmodels]
    yrips = ["cmip5."+thing.split(".")[1]+"."+thing.split(".")[2]+"."+thing.split(".")[3]+"." for thing in ymodels]
    in_common = np.intersect1d(xrips,yrips)
    L = len(in_common)
    Xt = MV.zeros((L,)+X.shape[1:])
    Yt = MV.zeros((L,)+Y.shape[1:])
    for i in range(L):
        modrip=in_common[i]
        ix = xrips.index(modrip)
        Xt[i] = X[ix]
        iy = yrips.index(modrip)
        Yt[i]=Y[iy]
    modax = cmip5.make_model_axis(in_common.tolist())
    Xt.id = X.id
    Xt.setAxis(0,modax)
    for axi in range(len(X.shape))[1:]:
        Xt.setAxis(axi,X.getAxis(axi))
    Yt.id = Y.id
    Yt.setAxis(0,modax)
    for axi in range(len(Y.shape))[1:]:
        Yt.setAxis(axi,Y.getAxis(axi))
    return Xt,Yt
        
def correlate_swcre_and_lcc(experiment):
    fl = cdms.open("MAPS/"+string.upper(experiment)+"_LCC.nc")
    lcc = fl("lcc")
    fl.close()
    fs = cdms.open("MAPS/"+string.upper(experiment)+"_SWCRE.nc")
    swcre = fs("SWCRE")
    fs.close()

    lcc_trop = tropical_marine_lcc_qu(lcc)
    swcre_trop = tropical_marine_lcc_qu(swcre)

    x,y= match_ensemble_members(lcc_trop,swcre_trop)
    Plotting.scatterplot_cmip(x,-y)
    plt.legend(ncol=3,fontsize=10,numpoints=1,loc=0)
    print genutil.statistics.correlation(x,-y)
    return x,-y
    
    
def compare_ECS_and_tropcloud(experiment,thing="LCC",ensemble_average=True):
    estimator = c.ECS(experiment)
    f = cdms.open("MAPS/"+string.upper(experiment)+"_"+thing+".nc")
    if thing is "LCC":
        lcc = f("lcc")
    else:
        lcc = f("SWCRE")
    stratocumulus = tropical_marine_lcc_qu(lcc)
    estimator.estimate_ECS(typ="annual regression")
    ecs = estimator.ECS
    if ensemble_average:
        ecs_ensav = cmip5.ensemble2multimodel(ecs)
        stratocum_ensav = cmip5.ensemble2multimodel(stratocumulus)
        X,Y = models_in_common(ecs_ensav,stratocum_ensav)
    else:
        X,Y = match_ensemble_members(ecs,stratocumulus)
    return X,Y
    
def scatterplot_stuff(X,Y,cmap=cm.viridis,unaveraged=False):
    markers = np.tile(["o","s","D","*","p","v"],100)
    if unaveraged:
        ed = cmip5.ensemble_dictionary(X)
        models = sorted(ed.keys())
        for i in range(len(models)):
            c=cmap(float(i)/float(len(models)))
            plt.plot(X.asma()[ed[models[i]]],Y.asma()[ed[models[i]]], markers[i],markersize=10,color=c,label=models[i])
    else:
        models = cmip5.models(X)
        for i in range(len(models)):
            c=cmap(float(i)/float(len(models)))
            plt.plot([X[i]],[Y[i]], markers[i],markersize=10,color=c,label=models[i])

class Sensitivity():
    def __init__(self,ensemble_average = True):
        estimator = c.ECS("amip")
        estimator.estimate_ECS(typ="annual regression")
        if ensemble_average:
            self.amip = cmip5.ensemble2multimodel(estimator.ECS)
        else:
            self.amip = estimator.ECS

        hestimator = c.ECS("historical")
        hestimator.estimate_ECS(typ="annual regression")
        if ensemble_average:
            self.historical = cmip5.ensemble2multimodel(hestimator.ECS)
        else:
            self.historical = hestimator.ECS

        if ensemble_average:
            self.amip,self.historical = models_in_common(self.amip,self.historical)
        if not ensemble_average:
            models = cmip5.models(cmip5.ensemble2multimodel(hestimator.ECS))
        else:
            models = cmip5.models(self.amip)

        equil = np.array([cmip5.clim_sens(model) for model in models])

        equil = MV.masked_where(np.isnan(equil),equil)
        equil.setAxis(0,self.historical.getAxis(0))
        self.equil = equil


def mask_sea_ice(X):
    fi = cdms.open("PATTERN_MAPS/sea_ice_for_amip.nc")
    SI = fi("sea_ice_percent")
    fi.close()
    SIr = SI.regrid(X.getGrid(),regridTool='regrid2')
    mask = SIr != 0
    if len(X.shape)==2:
        return MV.masked_where(mask,X)
    elif len(X.shape)==3:
        L = X.shape[0]
        bigmask = np.repeat(mask.asma()[np.newaxis,:,:],L,axis=0)
        
        return MV.masked_where(bigmask,X)
    else:
        print "input array must have dimensions (time, lat, lon)"
        raise TypeError
def mask_land(X):
    fl = cdms.open("PATTERN_MAPS/obs_climatology.nc")
    obs = fl("sst")
    mask = obs.mask
    fl.close()
    if len(X.shape)==2:
        return MV.masked_where(mask,X)
    elif len(X.shape)==3:
        L = X.shape[0]
        bigmask = np.repeat(mask[np.newaxis,:,:],L,axis=0)
        
        return MV.masked_where(bigmask,X)
    else:
        print "input array must have dimensions (time, lat, lon)"
        raise TypeError   


def sst_patterns(experiment,mask_the_land=True,mask_the_ice=True,remove_global_mean=True):
    variable = "ts"#REPLACE WITH TS!!!!!!!!!
    f = cdms.open("PATTERN_MAPS/"+experiment+"."+variable+".nc") #REPLACE WITH TS!!!!!!!!!
    tspatt = f(variable)
    f.close()
    if mask_the_land:
        tspatt = mask_land(tspatt)
    if mask_the_ice:
        tspatt = mask_sea_ice(tspatt)
    if remove_global_mean:
        globav=cdutil.averager(tspatt,axis='xy')
        tspatt = cmip5.cdms_clone(tspatt.asma()-globav.asma()[:,np.newaxis,np.newaxis],tspatt)
    return tspatt
def compare_ts_patterns(experiment="historical",func="covariance",mask_the_land=True,mask_the_ice=True,remove_global_mean=True,latbounds=(-90,90),lonbounds=(0,360)):
   

    variable = "ts"
    f = cdms.open("PATTERN_MAPS/"+experiment+"."+variable+".nc") 
    ts = f(variable)
    if mask_the_land:
        ts = mask_land(ts)
    if mask_the_ice:
        ts = mask_sea_ice(ts)
    f.close()
    #remove global mean
    if remove_global_mean:
        globav=cdutil.averager(ts,axis='xy')
        tspatt = cmip5.cdms_clone(ts.asma()-globav.asma()[:,np.newaxis,np.newaxis],ts)
    else:
        tspatt=ts

    fa  = cdms.open("PATTERN_MAPS/abrupt4xCO2.ts.nc")
    abrupt = fa("ts")
    if mask_the_land:
        abrupt = mask_land(abrupt)
    if mask_the_ice:
        abrupt=mask_sea_ice(abrupt)
        
    if remove_global_mean:
        globav=cdutil.averager(abrupt,axis='xy')
        abrupt = cmip5.cdms_clone(abrupt.asma()-globav.asma()[:,np.newaxis,np.newaxis],abrupt)
    fa.close()

   
    tspatt = tspatt(latitude=latbounds,longitude=lonbounds)
    abrupt = abrupt(latitude=latbounds,longitude=lonbounds)
    experiment_models = cmip5.models(tspatt)
    
    abrupt_models=cmip5.models(abrupt)
    new_abrupt = []
    for model in abrupt_models:
        trunc = model.split(".")[1]
        if trunc.find("GISS")<0:
            new_abrupt+=[trunc]
        else:
            rip = model.split(".")[3]
            phys_version = rip.split("p")[1]
            new_abrupt += [trunc+"*p"+phys_version]
  
    nmod = len(experiment_models)

    C = MV.zeros(nmod)+1.e20
    bad = []
    for i in range(nmod):
        model = experiment_models[i]
        trunc = model.split(".")[1]
       
        if trunc.find("GISS")<0:
            trunc=trunc
        else:
            rip = model.split(".")[3]
            phys_version = rip.split("p")[1]
            trunc = trunc+"*p"+phys_version
        try:
            corresponding_i = new_abrupt.index(trunc)
           
            C[i] = getattr(genutil.statistics,func)(tspatt[i],abrupt[corresponding_i],axis='xy')
           
        except:
           
            bad += [trunc]
    C = MV.masked_where(C>1.e10,C)
    C.setAxis(0,tspatt.getAxis(0))
    return C,bad

def correlation_histograms(func = "correlation",cmap=cm.magma,mask_the_land=True,mask_the_ice=True):
    Ch,badh = compare_ts_patterns(experiment="historical",func= func,mask_the_land=mask_the_land,mask_the_ice=mask_the_ice)
    Ca,bada = compare_ts_patterns(experiment="amip",func= func,mask_the_land=mask_the_land,mask_the_ice=mask_the_ice)
    
    plt.subplot(211)
    plt.hist(Ch.compressed(),color=cmap(.3),ec=cmap(.3),alpha=.7,label="Historical")
    plt.hist(Ca.compressed(),color=cmap(.6),ec=cmap(.6),alpha=.7,label="AMIP")
    plt.title(func+" with model's own abrupt4xCO2 pattern")
    plt.legend()
    plt.subplot(212)
    plt.title(func+" with multimodel mean abrupt4xCO2 pattern")
    abruptsst = sst_patterns("abrupt4xCO2",mask_the_land=mask_the_land,mask_the_ice=mask_the_ice)
    abrupt_avg = MV.average(abruptsst,axis=0)

    historicalsst = sst_patterns("historical",remove_global_mean=True,mask_the_land=mask_the_land,mask_the_ice=mask_the_ice)
    historical_avg = MV.average(historicalsst,axis=0)
    
    amipsst = sst_patterns("amip",remove_global_mean=True,mask_the_land=mask_the_land,mask_the_ice=mask_the_ice)
    amip_avg = MV.average(amipsst,axis=0)
    
    historical_corr_with_mean=np.array([float(getattr(genutil.statistics,func)(historicalsst[i],abrupt_avg,axis='xy')) for i in range(historicalsst.shape[0])])
    historical_corr_with_mean=MV.array(np.ma.masked_where(np.isnan(historical_corr_with_mean),historical_corr_with_mean))
    historical_corr_with_mean.setAxis(0,Ch.getAxis(0))
     
    amip_corr_with_mean=np.array([float(getattr(genutil.statistics,func)(amipsst[i],abrupt_avg,axis='xy')) for i in range(amipsst.shape[0])])
    amip_corr_with_mean=MV.array(np.ma.masked_where(np.isnan(amip_corr_with_mean),amip_corr_with_mean))
    amip_corr_with_mean.setAxis(0,Ca.getAxis(0))

    plt.hist(historical_corr_with_mean.compressed(),color=cmap(.3),ec=cmap(.3),alpha=.7)
    plt.hist(amip_corr_with_mean.compressed(),color=cmap(.6),ec=cmap(.6),alpha=.7)

    return Ch,historical_corr_with_mean,Ca,amip_corr_with_mean
import Clean_Code as cc
def equil_minus_hist_ecs():
    historical = cc.estimate_ECS("historical")
    equil=cc.estimate_ECS("equil")
    models = cmip5.models(historical)
    eqmodels = cmip5.models(equil)
    L = len(models)
    df = MV.zeros(L)+1.e20
    for i in range(L):
        trunc = models[i].split(".")[1]
        if trunc in eqmodels:
            df[i] = equil[eqmodels.index(trunc)] - historical[i]
    df = MV.masked_where(df>1.e10,df)
    df.setAxis(0,historical.getAxis(0))
    return df
    
    
    
def SO_historical_vs_SO_abrupt(latbounds=(-90,-50),mask_all=False,remove_global_mean=True):

    experiment="historical"
    variable = "ts"
    f = cdms.open("PATTERN_MAPS/"+experiment+"."+variable+".nc") 
    ts = f(variable)
    if mask_all:
        ts = mask_land(mask_sea_ice(ts))
    f.close()
    #remove global mean
    if remove_global_mean:
        globav=cdutil.averager(ts,axis='xy')
        historical = cmip5.cdms_clone(ts.asma()-globav.asma()[:,np.newaxis,np.newaxis],ts)
    else:
        historical=ts

    fa  = cdms.open("PATTERN_MAPS/abrupt4xCO2.ts.nc")
    equil = fa("ts")
    if mask_all:
        equil = mask_land(mask_sea_ice(equil))
    if remove_global_mean:
        globav=cdutil.averager(equil,axis='xy')
        equil = cmip5.cdms_clone(equil.asma()-globav.asma()[:,np.newaxis,np.newaxis],equil)
    fa.close()

    equil=cdutil.averager(equil(latitude=latbounds),axis='xy')
    historical = cdutil.averager(historical(latitude=latbounds),axis='xy')
    

    models = cmip5.models(historical)
    eqmodels = [x.split(".")[1] for x in cmip5.models(equil)]
    L = len(models)
    df = MV.zeros(L)+1.e20
    for i in range(L):
        trunc = models[i].split(".")[1]
       
        if trunc in eqmodels:
            df[i] = equil[eqmodels.index(trunc)] - historical[i]
    df = MV.masked_where(df>1.e10,df)
    df.setAxis(0,historical.getAxis(0))
    return df

def play_with_SO(latbounds):
    sodiff = SO_historical_vs_SO_abrupt(latbounds=latbounds)
    tdiff = equil_minus_hist_ecs()
    plt.subplot(211)
    Plotting.scatterplot_cmip(tdiff,sodiff)
    plt.title(str(genutil.statistics.correlation(tdiff.asma(),sodiff)))

    tdiff_ensav=cmip5.ensemble2multimodel(tdiff)
    sodiff_ensav=cmip5.ensemble2multimodel(sodiff)
    plt.subplot(212)
    Plotting.scatterplot_cmip(tdiff_ensav,sodiff_ensav)
    plt.title(str(genutil.statistics.correlation(tdiff_ensav.asma(),sodiff_ensav)))

    
    
def check_amip_lcc():
    f = cdms.open("REGRIDDED_CLOUD/cmip5.ensemble.amip.r1i1pALL.ann.atm.Amon.cloudbins.ver-1.latestX.nc")
    lcc = 100*(f("low_cloud")+f("mid_cloud"))
    lcc_diff = last_ten_minus_first_ten(lcc)
    
 
    

        
    
