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

def historical_global_mean(x):
    start = '1996-1-1'
    stop = '2006-1-1'
    data = x(time=(start,stop))
    
    return MV.average(cdutil.averager(data,axis='xy'),axis=0)


def generate_global_means(variables):
    hdirec = "/work/cmip5/historical/atm/mo/"
    adirec = "/work/cmip5/amip/atm/mo/"
    
    for variable in variables:
        hpath = hdirec + variable+"/"
        apath = adirec+ variable+"/"
        hwrite = cdms.open("DATA/cmip5.historical."+variable+".nc","w")
        awrite = cdms.open("DATA/cmip5.amip."+variable+".nc","w")
        va = cmip5.get_ensemble(apath,variable,func=historical_global_mean)
        va.id = variable
        va.name = variable
        awrite.write(va)
        
        vh = cmip5.get_ensemble(hpath,variable,func=historical_global_mean)
        vh.id = variable
        vh.name = variable
        hwrite.write(vh)
        hwrite.close()
        awrite.close()

def generate_amip_global_means(variables):
    hdirec = "/work/cmip5/amip4K/atm/mo/"
    adirec = "/work/cmip5/amipFuture/atm/mo/"
    
    for variable in variables:
        hpath = hdirec + variable+"/"
        apath = adirec+ variable+"/"
        hwrite = cdms.open("DATA/cmip5.amip4K."+variable+".nc","w")
        awrite = cdms.open("DATA/cmip5.amipFuture."+variable+".nc","w")
        va = cmip5.get_ensemble(apath,variable,func=historical_global_mean)
        va.id = variable
        va.name = variable
        awrite.write(va)
        
        vh = cmip5.get_ensemble(hpath,variable,func=historical_global_mean)
        vh.id = variable
        vh.name = variable
        hwrite.write(vh)
        hwrite.close()
        awrite.close()

def historical_global_mean_time(x):
    start = '1979-1-1'
    stop = '2006-1-1'
    data = x(time=(start,stop))
    
    return cdutil.averager(data,axis='xy')


def generate_global_mean_timeseries(experiment,variable):
    hdirec = "/work/cmip5/"+experiment+"/atm/mo/"
    
    hpath = hdirec + variable+"/"
    
    hwrite = cdms.open("DATA/TIMESERIES/cmip5."+experiment+"."+variable+".nc","w")
        
    vh = cmip5.get_ensemble(hpath,variable,func=historical_global_mean_time)
    vh.id = variable
    vh.name = variable
    hwrite.write(vh)
    hwrite.close()
    


def TOA_imbalance(typ):
    variable = "rsdt"
    fwrite = cdms.open("DATA/cmip5."+typ+"."+variable+".nc")
    RSDT = fwrite(variable)
    rsdt_models = eval(fwrite(variable).getAxis(0).models)
    rsdt_trunc = [x.split("/")[-1].split(variable)[0] for x in rsdt_models]
    fwrite.close()

    variable = "rsut"
    fwrite = cdms.open("DATA/cmip5."+typ+"."+variable+".nc")
    rsut_models = eval(fwrite(variable).getAxis(0).models)
    rsut_trunc = [x.split("/")[-1].split(variable)[0] for x in rsut_models]
    RSUT = fwrite(variable)
    fwrite.close()

    
    variable = "rlut"
    fwrite = cdms.open("DATA/cmip5."+typ+"."+variable+".nc")
    rlut_models = eval(fwrite(variable).getAxis(0).models)
    rlut_trunc = [x.split("/")[-1].split(variable)[0] for x in rlut_models]
    RLUT = fwrite(variable)
    fwrite.close()

    ok_models = np.intersect1d(np.intersect1d(rlut_trunc,rsut_trunc),rsdt_trunc)
    d={}
    for model in ok_models:
        rsdt_mod = RSDT[rsdt_trunc.index(model)]
        rsut_mod = RSUT[rsut_trunc.index(model)]
        rlut_mod = RLUT[rlut_trunc.index(model)]
        d[model] = rsdt_mod - (rlut_mod+rsut_mod)
    return d
from matplotlib.patches import Ellipse        
def scatterplot(cmap=cm.Set1,xax = "historical",yax = "amip"):
    H = TOA_imbalance(xax)
    A = TOA_imbalance(yax)
    x = []
    y = []
    k=[]
    for Hkey in sorted(H.keys()):
        Akey = Hkey.replace(xax,yax)
        if Akey in A.keys():
            x+=[H[Hkey]]
            y+=[A[Akey]]
            k+=[Akey]
    x=np.ma.array(x)
    x = np.ma.masked_where(np.isnan(x),x)
    y=np.ma.array(y)
    y = np.ma.masked_where(np.isnan(y),y)
    allmodels = np.array([thing.split(".")[1] for thing in k])
    models = np.unique(allmodels)
    L = float(len(models))
    i=0
    dummy=[]
    for model in models:
        I = np.where(allmodels == model)[0]
        
        if not (True in y[I].mask):
            if len(I) == 1:
                plt.plot(x[I],y[I],"o",color=cmap(i/L),label=model,markersize=10)
            else:
                xi=x[I]
                yi=y[I]
                xy=(np.ma.average(xi),np.ma.average(yi))
                stuff,=plt.plot(xi,yi,"o",color=cmap(i/L),label=model,markersize=10)
                dummy+=[stuff]
                
                ell = Ellipse(xy=xy,width=(np.max(xi)-np.min(xi)),height = (np.max(yi)-np.min(yi)),angle=0)
                ell.set_facecolor(cmap(i/L))
                ell.set_label(model)
                plt.gca().add_artist(ell)
        i+=1
    plt.xlim(-2,6)
    plt.ylim(-2,6)
    plt.legend(loc=0,ncol=2,numpoints=1)
    [stuff.set_visible(False) for stuff in dummy]
    
            


if __name__ == "__main__":
    surface = {"hfls": "Surface Upward Latent Heat Flux",\
                "hfss": "Surface Upward Sensible Heat Flux",\
                "rlds":"Surface Downwelling Longwave Radiation",\
                "rlus":"Surface Upwelling Longwave Radiation",\
                "rsds":"Surface Downwelling Shortwave Radiation",\
                "rsus": "Surface Upwelling Shortwave Radiation"}
    TOA = {"rsdt":"TOA Incident Shortwave Radiation",\
            "rsut": "TOA Outgoing Shortwave Radiation",\
            "rlut": "TOA Outgoing Longwave Radiation"}

    experiments = ["historical","amip","amip4K","amipFuture"]
    variables = TOA.keys()
    for variable in variables:
            for experiment in experiments:
                    generate_global_mean_timeseries(experiment,variable)
    
