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

def Q(typ):
    """
    Get the TOA radiative imbalance from downward SW - (upward SW + upward LW)
    inputs: typ must be one of ["amip","amip4K","amipFuture","historical"]
    
    """
    variable = "rsdt"
    fwrite = cdms.open("DATA/TIMESERIES/cmip5."+typ+"."+variable+".nc")
    #RSDT_total = fwrite(variable)
    RSDT = fwrite(variable)
    rsdt_models = eval(fwrite(variable).getAxis(0).models)
    rsdt_trunc = [x.split("/")[-1].split(variable)[0] for x in rsdt_models]
    fwrite.close()

    variable = "rsut"
    fwrite = cdms.open("DATA/TIMESERIES/cmip5."+typ+"."+variable+".nc")
    rsut_models = eval(fwrite(variable).getAxis(0).models)
    rsut_trunc = [x.split("/")[-1].split(variable)[0] for x in rsut_models]
    RSUT = fwrite(variable)
    fwrite.close()

    
    variable = "rlut"
    fwrite = cdms.open("DATA/TIMESERIES/cmip5."+typ+"."+variable+".nc")
    rlut_models = eval(fwrite(variable).getAxis(0).models)
    rlut_trunc = [x.split("/")[-1].split(variable)[0] for x in rlut_models]
    RLUT = fwrite(variable)
    fwrite.close()

    ok_models = np.intersect1d(np.intersect1d(rlut_trunc,rsut_trunc),rsdt_trunc)

    toa = MV.zeros((len(ok_models),RSDT.shape[1]))
    i=0
    for model in ok_models:
        rsdt_mod = RSDT[rsdt_trunc.index(model)]
        rsut_mod = RSUT[rsut_trunc.index(model)]
        rlut_mod = RLUT[rlut_trunc.index(model)]
        test = rsdt_mod - (rlut_mod+rsut_mod)
        toa[i] = test
        i+=1
    modax = cdms.createAxis(range(len(ok_models)))
    modax.models = str(ok_models.tolist())
    toa.setAxis(0,modax)
    toa.setAxis(1,RLUT.getTime())
    return toa

def T(typ):
    variable = "tas"
    fwrite = cdms.open("DATA/TIMESERIES/cmip5."+typ+"."+variable+".nc")
    tas_models = eval(fwrite(variable).getAxis(0).models)
    tas_trunc = [x.split("/")[-1].split(variable)[0] for x in tas_models]
    temperature = fwrite(variable)
    toa = Q(typ)
    models = np.intersect1d( tas_trunc, eval(toa.getAxis(0).models))
    tas = MV.zeros((len(models),temperature.shape[1]))
    j=0
    for model in models:
        i = tas_trunc.index(model)
        tas[j] = temperature[i]
        j+=1
    modax = cdms.createAxis(range(len(models)))
    modax.models = str(models.tolist())
    tas.setAxis(0,modax)
    tas.setAxis(1,temperature.getTime())
    return tas



def SWCRE(typ):
    variable = "rsut"
    fwrite = cdms.open("DATA/TIMESERIES/cmip5."+typ+"."+variable+".nc")
    rsut_models = eval(fwrite(variable).getAxis(0).models)
    rsut_trunc = [x.split("/")[-1].split(variable)[0] for x in rsut_models]
    RSUT = fwrite(variable)
    fwrite.close()


    variable = "rsutcs"
    fwrite = cdms.open("DATA/TIMESERIES/cmip5."+typ+"."+variable+".nc")
    rsutcs_models = eval(fwrite(variable).getAxis(0).models)
    rsutcs_trunc = [x.split("/")[-1].split(variable)[0] for x in rsutcs_models]
    RSUTCS = fwrite(variable)
    fwrite.close()

    toa = Q(typ)
    models = np.intersect1d(eval(toa.getAxis(0).models),rsutcs_trunc)
    sw = MV.zeros((len(models),RSUTCS.shape[1]))
    j=0
    for model in models:
        alli = rsut_trunc.index(model)
        csi = rsutcs_trunc.index(model)
        sw[j] = RSUTCS[csi] - RSUT[alli]
        j+=1
    modax = cdms.createAxis(range(len(models)))
    modax.models = str(models.tolist())
    sw.setAxis(0,modax)
    sw.setAxis(1,RSUT.getTime())

    return sw
def get_forcing(typ="giss"):
    if typ == "giss":
    # Read forcing from Ron Miller's forcing calculations (cf Miller et al, JAMES 2013)
        filename="historical.lpl"
        ffile = open(filename)
        lines = ffile.readlines()
        years = [float(x.split()[0]) for x in lines[4:]]
        forcing = [float(x.split()[1]) for x in lines[4:]]
        forcing = MV.array(forcing)
        forcing = forcing - np.average(forcing[:10])
        forcing.units = "W m-2"
        tax = cdms.createAxis(years)
        tax.designateTime()
        tax.units = 'years since 0000-07-1'
        forcing.setAxis(0,tax)
        forcing.id = "forcing"
    else:
        far5 = cdms.open("AR5_forcings.nc","r")
        forcing = far5("total")
        forcing.id="forcing"
    return forcing



class ECS():
    def __init__(self,typ):
        self.typ = typ
        self.T = cdutil.YEAR(T(typ)) 
        self.Q = cdutil.YEAR(Q(typ))
        if typ != "piControl":
            self.F = get_forcing(typ="AR5")(time=('1979-1-1','2005-12-31'))
        else:
            self.F = MV.zeros(27)
            Qamip=cdutil.YEAR(Q("amip"))
            nmod,nt = self.Q.shape
            num_each = nt/27
            Qre = self.Q[:,:num_each*27].reshape(nmod*num_each,27)
            Tre = self.T[:,:num_each*27].reshape(nmod*num_each,27)
            self.Q = Qre
            self.T = Tre
            self.Q.setAxis(1,Qamip.getTime())
            self.T.setAxis(1,Qamip.getTime())
            
        self.dT = cmip5.cdms_clone(self.T.asma() - MV.average(self.T[:,:10],axis=1).asma()[:,np.newaxis],self.T)
        self.dF = self.F - MV.average(self.F[:10])
        self.dQ =  cmip5.cdms_clone(self.Q.asma() - MV.average(self.Q[:,:10],axis=1).asma()[:,np.newaxis],self.Q)

    def estimate_ECS(self,decadal=True):
        if decadal:
            y = genutil.filters.runningaverage(self.dF-self.dQ,10,axis=1)
            x = genutil.filters.runningaverage(self.dT,10,axis=1)
        else:
            y = self.dF-self.dQ
            x = self.dT
        F2 = 3.7
        lam = np.ma.zeros(x.shape[0])+1.e20
        for i in range(lam.shape[0]):
            try:
                L = np.polyfit(x[i],y[i],1)[0]
                lam[i] = F2/L
            except:
                continue
        ECS = MV.array(np.ma.masked_where(np.abs(lam)>1.e10,lam))
        ECS.setAxis(0,x.getAxis(0))
        self.ECS = ECS
            
        

                
def correspond_models():
     historical = ECS("historical")
     amip = ECS("amip")
     #hmodels = historical.
