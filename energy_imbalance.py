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
        
def scatterplot(cmap=cm.viridis):
    H = TOA_imbalance("historical")
    A = TOA_imbalance("amip")
    x = []
    y = []
    k=[]
    for Hkey in sorted(H.keys()):
        Akey = Hkey.replace("historical","amip")
        if Akey in A.keys():
            x+=[H[Hkey]]
            y+=[A[Akey]]
            k+=[Akey]
    x=np.array(x)
    y=np.array(y)
    allmodels = np.array([x.split(".")[1] for x in k])
    models = np.unique(allmodels)
    L = float(len(models))
    i=0
    for model in models:
        I = np.where(allmodels == model)[0]
        plt.plot(x[I],y[I],"o",color=cmap(i/L))
        i+=1
    return x,y,k
            


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
    variables = surface.keys()
    generate_global_means(variables)
