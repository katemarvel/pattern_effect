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


def piC_global_mean_time(x):
    
    data = x[:200*12]
    
    return cdutil.averager(data,axis='xy')

def generate_global_mean_timeseries(experiment,variable):
    hdirec = "/work/cmip5/"+experiment+"/atm/mo/"
    
    hpath = hdirec + variable+"/"
    
    hwrite = cdms.open("DATA/TIMESERIES/cmip5."+experiment+"."+variable+".nc","w")
    if (experiment == "historical" or experiment == "amip"):
        func = historical_global_mean_time
    elif experiment == "piControl":
        func = piC_global_mean_time
    vh = cmip5.get_ensemble(hpath,variable,func=func)
    vh.id = variable
    vh.name = variable
    hwrite.write(vh)
    hwrite.close()
    


def TOA_imbalance_dec(typ):
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

def TOA_imbalance(typ):
    variable = "rsdt"
    fwrite = cdms.open("DATA/TIMESERIES/cmip5."+typ+"."+variable+".nc")
    #RSDT_total = fwrite(variable)
    RSDT = MV.average(fwrite(variable,time=('1996-1-1','2005-12-31')),axis=1)-MV.average(fwrite(variable,time=('1979-1-1','1988-12-31')),axis=1)
    rsdt_models = eval(fwrite(variable).getAxis(0).models)
    rsdt_trunc = [x.split("/")[-1].split(variable)[0] for x in rsdt_models]
    fwrite.close()

    variable = "rsut"
    fwrite = cdms.open("DATA/TIMESERIES/cmip5."+typ+"."+variable+".nc")
    rsut_models = eval(fwrite(variable).getAxis(0).models)
    rsut_trunc = [x.split("/")[-1].split(variable)[0] for x in rsut_models]
    RSUT = MV.average(fwrite(variable,time=('1996-1-1','2005-12-31')),axis=1)-MV.average(fwrite(variable,time=('1979-1-1','1988-12-31')),axis=1)#fwrite(variable)
    fwrite.close()

    
    variable = "rlut"
    fwrite = cdms.open("DATA/TIMESERIES/cmip5."+typ+"."+variable+".nc")
    rlut_models = eval(fwrite(variable).getAxis(0).models)
    rlut_trunc = [x.split("/")[-1].split(variable)[0] for x in rlut_models]
    RLUT = MV.average(fwrite(variable,time=('1996-1-1','2005-12-31')),axis=1)-MV.average(fwrite(variable,time=('1979-1-1','1988-12-31')),axis=1)#fwrite(variable)
    fwrite.close()

    ok_models = np.intersect1d(np.intersect1d(rlut_trunc,rsut_trunc),rsdt_trunc)
    d={}
    for model in ok_models:
        rsdt_mod = RSDT[rsdt_trunc.index(model)]
        rsut_mod = RSUT[rsut_trunc.index(model)]
        rlut_mod = RLUT[rlut_trunc.index(model)]
        d[model] = rsdt_mod - (rlut_mod+rsut_mod)
    return d
def TOA_imbalance_whole_time(typ):
    variable = "rsdt"
    fwrite = cdms.open("DATA/TIMESERIES/cmip5."+typ+"."+variable+".nc")
    #RSDT_total = fwrite(variable)
    RSDT = MV.average(fwrite(variable,time=('1979-1-1','2005-12-31')),axis=1)
    rsdt_models = eval(fwrite(variable).getAxis(0).models)
    rsdt_trunc = [x.split("/")[-1].split(variable)[0] for x in rsdt_models]
    fwrite.close()

    variable = "rsut"
    fwrite = cdms.open("DATA/TIMESERIES/cmip5."+typ+"."+variable+".nc")
    rsut_models = eval(fwrite(variable).getAxis(0).models)
    rsut_trunc = [x.split("/")[-1].split(variable)[0] for x in rsut_models]
    RSUT = MV.average(fwrite(variable,time=('1979-1-1','2005-12-31')),axis=1)
    fwrite.close()

    
    variable = "rlut"
    fwrite = cdms.open("DATA/TIMESERIES/cmip5."+typ+"."+variable+".nc")
    rlut_models = eval(fwrite(variable).getAxis(0).models)
    rlut_trunc = [x.split("/")[-1].split(variable)[0] for x in rlut_models]
    RLUT = MV.average(fwrite(variable,time=('1979-1-1','2005-12-31')),axis=1)
    fwrite.close()

    ok_models = np.intersect1d(np.intersect1d(rlut_trunc,rsut_trunc),rsdt_trunc)
    d={}
    for model in ok_models:
        rsdt_mod = RSDT[rsdt_trunc.index(model)]
        rsut_mod = RSUT[rsut_trunc.index(model)]
        rlut_mod = RLUT[rlut_trunc.index(model)]
        d[model] = rsdt_mod - (rlut_mod+rsut_mod)
    return d


def scatterplot_whole_time():
    xax = "amipFuture"
    yax = "amip4K"
    H = TOA_imbalance_whole_time(xax)
    badH = {}
    He = {}
    A = TOA_imbalance_whole_time(yax)
    badA = {}
    Ae = {}
    
    
    scatterx = []
    scattery = []
    k=[]
    for Hkey in sorted(H.keys()):
        Akey = Hkey.replace(xax,yax)
        if not (Akey in A.keys()):
            badH[Hkey] = H.pop(Hkey)
    for Akey in sorted(A.keys()):
        Hkey = Akey.replace(yax,xax)
        if not (Hkey in H.keys()):
            badA[Akey] = A.pop(Akey)
    Hnew = {}
    Anew = {}
    amip_for_subtraction = TOA_imbalance_whole_time("amip")
    for key in H.keys():
        amip_key = key.replace("Future","")
        if amip_key in amip_for_subtraction.keys():
            Hnew[key] = H[key]-amip_for_subtraction[amip_key]
            Anew[key.replace("Future","4K")] = A[key.replace("Future","4K")] - amip_for_subtraction[amip_key]
    return Hnew,Anew
   
    
   
def scatterplot_ensemble_average(cmap=cm.Set1,xax = "historical",yax = "amip"):
    H = TOA_imbalance(xax)
    badH = {}
    He = {}
    A = TOA_imbalance(yax)
    badA = {}
    Ae = {}
    
    
    scatterx = []
    scattery = []
    k=[]
    for Hkey in sorted(H.keys()):
        Akey = Hkey.replace(xax,yax)
        if not (Akey in A.keys()):
            badH[Hkey] = H.pop(Hkey)
    for Akey in sorted(A.keys()):
        Hkey = Akey.replace(yax,xax)
        if not (Hkey in H.keys()):
            badA[Akey] = A.pop(Akey)
    Hmodels = np.array([x.split(".")[1]+".p"+x.split(".")[3].split("p")[-1]  for x in sorted(H.keys())])
    Amodels = np.array([x.split(".")[1]+".p"+x.split(".")[3].split("p")[-1]  for x in sorted(A.keys())])
    nmod = len(np.unique(Hmodels))
    Hkeys_good = np.array(sorted(H.keys()))
    
    Akeys_good = np.array(sorted(A.keys()))
    for mod in np.unique(Hmodels):
        I = np.where(Hmodels == mod)[0]
        
        Hensemble = np.ma.zeros(len(I))
        counter = 0
        for key in Hkeys_good[I]:
            Hensemble[counter] = H[key]
            counter +=1
        He[mod]=Hensemble
        scatterx += [np.ma.average(Hensemble)]
        I = np.where(Amodels == mod)[0]
        Aensemble = np.ma.zeros(len(I))
        counter = 0
        for key in Akeys_good[I]:
            Aensemble[counter] = A[key]
            counter +=1
        scattery += [float(np.ma.average(Aensemble))]
        Ae[mod]=Aensemble
        
    X = np.ma.masked_where(np.isnan(scatterx),scatterx)
    Y = np.ma.masked_where(np.isnan(scattery),scattery)
    markers = ["o","s","d","8","*","H","<","^","v","h","","o","s","d","8","*","H","<","^","v","h","D"]
    for i in range(len(X)):
        plt.scatter([X[i]],[Y[i]],marker=markers[i],color=cmap(i/float(len(X))),label=np.unique(Hmodels)[i],s=100)
    plt.legend(numpoints=1,fontsize=8,ncol=3,loc=0)
    axmax=max(max(plt.xlim()),max(plt.ylim()))
    axmin=min(min(plt.xlim()),min(plt.ylim()))
    plt.xlim((axmin,axmax))
    plt.ylim((axmin,axmax))
    sx = np.linspace(axmax,axmin)
    plt.plot(sx,sx,"k:")
    plt.xlabel("TOA energy imbalance (1996-2005)-(1979-1988) in W/m2: "+ xax)
    plt.ylabel("TOA energy imbalance (1996-2005)-(1979-1988) in W/m2: "+ yax)
    return X,Y,He,Ae
    #plt.xlim(-.1,.6)
    #plt.ylim(-.1,.6)
                
def scatterplot_r1(cmap=cm.Set1,xax = "historical",yax = "amip"):
    H = TOA_imbalance(xax)
    badH = {}
    A = TOA_imbalance(yax)
    badA = {}
    
    scatterx = []
    scattery = []
    k=[]
    for Hkey in sorted(H.keys()):
        Akey = Hkey.replace(xax,yax)
        if not (Akey in A.keys()):
            badH[Hkey] = H.pop(Hkey)
    for Akey in sorted(A.keys()):
        Hkey = Akey.replace(yax,xax)
        if not (Hkey in H.keys()):
            badA[Akey] = A.pop(Akey)
    Hmodels = np.array([x.split(".")[1]+".p"+x.split(".")[3].split("p")[-1]  for x in sorted(H.keys())])
    Amodels = np.array([x.split(".")[1]+".p"+x.split(".")[3].split("p")[-1]  for x in sorted(A.keys())])
    nmod = len(np.unique(Hmodels))
    Hkeys_good = np.array(sorted(H.keys()))
    
    Akeys_good = np.array(sorted(A.keys()))
    for mod in np.unique(Hmodels):
        I = np.where(Hmodels == mod)[0]
        
       
        
        scatterx += [H[Hkeys_good[I[0]]]]
        I = np.where(Amodels == mod)[0]
        Aensemble = np.ma.zeros(len(I))
        counter = 0
        

        
        
    X = np.ma.masked_where(np.isnan(scatterx),scatterx)
    Y = np.ma.masked_where(np.isnan(scattery),scattery)
    markers = ["o","s","d","8","*","H","<","^","v","h","X","o","s","d","8","*","H","<","^","v","h","X"]
    for i in range(len(X)):
        plt.scatter([X[i]],[Y[i]],marker=markers[i],color=cmap(i/float(len(X))),label=np.unique(Hmodels)[i],s=100)
    plt.legend(numpoints=1,fontsize=8,ncol=3,loc=0)
    axmax=max(max(plt.xlim()),max(plt.ylim()))
    axmin=min(min(plt.xlim()),min(plt.ylim()))
    plt.xlim((axmin,axmax))
    plt.ylim((axmin,axmax))
    sx = np.linspace(axmax,axmin)
    plt.plot(sx,sx,"k:")
    plt.xlabel("TOA energy imbalance (1996-2005)-(1979-1988) in W/m2: "+ xax)
    plt.ylabel("TOA energy imbalance (1996-2005)-(1979-1988) in W/m2: "+ yax)
    return X,Y,H,A    
from matplotlib.patches import Ellipse        
def scatterplot(cmap=cm.Set1,xax = "historical",yax = "amip"):
    H = TOA_imbalance_dec(xax)
    A = TOA_imbalance_dec(yax)
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
                
                ell = Ellipse(xy=xy,width=(np.max(xi)-np.min(xi)),height = (np.max(yi)-np.min(yi)),angle=0,alpha=.5)
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

    experiments = ["historical","amip","amip4K","amipFuture","piControl"]
    variables = TOA.keys()+surface.keys()+["tas"]
    for variable in variables:
            for experiment in experiments:
                    generate_global_mean_timeseries(experiment,variable)
    
def get_forcing(typ="giss"):
    if typ == "giss":
    # Read forcing from Ron Miller's forcing calculations (cf Miller et al, JAMES 2013)
        filename="/Users/kmarvel/Google Drive/ECS/FORCINGS/historical.lpl"
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


def giss_only(typ="amip",year=True):
    variable = "rsdt"
    fwrite = cdms.open("DATA/TIMESERIES/cmip5."+typ+"."+variable+".nc")
    rsdt_models = eval(fwrite(variable).getAxis(0).models)
    giss = np.where([x.find("GISS-E2-R.")>0 for x in rsdt_models])[0]
    rips = [x.split(".")[3] for x in rsdt_models]
    p1 = np.where([x[-2:]=='p1' for x in rips])[0]
    gissp1=np.intersect1d(giss,p1)
    
   
    RSDT_all = fwrite(variable)
    RSDT = MV.array(RSDT_all.asma()[gissp1])
    for i in range(len(RSDT_all.getAxisIds()))[1:]:
        RSDT.setAxis(i,RSDT_all.getAxis(i))
    
    fwrite.close()

    variable = "rsut"
    fwrite = cdms.open("DATA/TIMESERIES/cmip5."+typ+"."+variable+".nc")
    rsut_models = eval(fwrite(variable).getAxis(0).models)
    giss = np.where([x.find("GISS-E2-R.")>0 for x in rsut_models])[0]
    rips = [x.split(".")[3] for x in rsut_models]
    p1 = np.where([x[-2:]=='p1' for x in rips])[0]
    gissp1=np.intersect1d(giss,p1)
    
   
    RSUT_all = fwrite(variable)
    RSUT = MV.array(RSUT_all.asma()[gissp1])
    for i in range(len(RSUT_all.getAxisIds()))[1:]:
        RSUT.setAxis(i,RSUT_all.getAxis(i))
    
    fwrite.close()

    
    variable = "rlut"
    fwrite = cdms.open("DATA/TIMESERIES/cmip5."+typ+"."+variable+".nc")
    rlut_models = eval(fwrite(variable).getAxis(0).models)
    giss = np.where([x.find("GISS-E2-R.")>0 for x in rlut_models])[0]
    rips = [x.split(".")[3] for x in rlut_models]
    p1 = np.where([x[-2:]=='p1' for x in rips])[0]
    gissp1=np.intersect1d(giss,p1)
    
   
    RLUT_all = fwrite(variable)
    RLUT = MV.array(RLUT_all.asma()[gissp1])
    for i in range(len(RLUT_all.getAxisIds()))[1:]:
        RLUT.setAxis(i,RLUT_all.getAxis(i))
    
    fwrite.close()


    variable = "tas"
    fwrite = cdms.open("DATA/TIMESERIES/cmip5."+typ+"."+variable+".nc")
    tas_models = eval(fwrite(variable).getAxis(0).models)
    giss = np.where([x.find("GISS-E2-R.")>0 for x in tas_models])[0]
    rips = [x.split(".")[3] for x in tas_models]
    p1 = np.where([x[-2:]=='p1' for x in rips])[0]
    gissp1=np.intersect1d(giss,p1)
    
   
    TAS_all = fwrite(variable)
    TAS = MV.array(TAS_all.asma()[gissp1])
    for i in range(len(TAS_all.getAxisIds()))[1:]:
        TAS.setAxis(i,TAS_all.getAxis(i))
    
    fwrite.close()

    toa = RSDT-(RSUT+RLUT)
    cdutil.setTimeBoundsMonthly(toa)
    cdutil.setTimeBoundsMonthly(TAS)
    if year:
        return cdutil.YEAR(toa), cdutil.YEAR(TAS)
 
        
    else:
        return toa,TAS


def write_ar5_netcdf():
    f = open("../ECS/AR5Forcings.csv")
    data = f.readlines()
    stuff = data[0].split("\r")
    headings = stuff[0].split(",")
    years = [float(x.split(",")[0]) for x in stuff[1:]]
    tax = cdms.createAxis(years)
    tax.designateTime()
    tax.units='years since 0000-1-1'
    far5 = cdms.open("AR5_forcings.nc","w")
    for forcingtype in headings:
        i = headings.index(forcingtype)
        fdata = MV.array([float(x.split(",")[i]) for x in stuff[1:]])
        fdata.units = "W/m2"
        fdata.id = ".".join(forcingtype.split(" "))
        fdata.name = forcingtype
        fdata.setAxis(0,tax)
        far5.write(fdata)
    totalforcings = MV.array([sum(map(float,x.split(","))[1:]) for x in stuff[1:]])
    totalforcings.id="total"
    totalforcings.name = "total"
    totalforcings.setAxis(0,tax)
    far5.write(totalforcings)
    f.close()
    far5.close()
    
    
def get_ecs_estimates(typ="amip",returndt=False,returndq=False):
    toa,t = giss_only(typ,year=True)
    dt =  MV.average(t[:,-10:],axis=1) - MV.average(t[:,:10],axis=1)
    dq = MV.average(toa[:,-10:],axis=1) - MV.average(toa[:,:10],axis=1)
    forcing = get_forcing(typ="giss")(time=('1979-1-1','2005-12-31'))
    df = MV.average(forcing[-10:]) - MV.average(forcing[:10])
    if returndt:
        return dt
    if returndq:
        return dq
    else:
        return dt/(df-dq)

def get_all_ecs_estimates(typ="amip"):
    f = cdms.open("DATA/TIMESERIES/cmip5."+typ+".tas.nc")
    TAS = f("tas")(time=('1979-1-1','2005-12-31'))
    tas_keys = [x.split("/")[-1].split("tas")[0] for x in eval(TAS.getAxis(0).models)]

    toa = TOA_imbalance(typ)
    forcing = get_forcing(typ="ar5")(time=('1979-1-1','2005-12-31'))
    df = MV.average(forcing[-10:]) - MV.average(forcing[:10])
    sens={}
    for key in np.intersect1d(tas_keys,toa.keys()):
        i = tas_keys.index(key)
        t = TAS[i]
        dt =  MV.average(t[-10:]) - MV.average(t[:10])
        dq = toa[key]
        sens[key] = [dt,dq]
    return sens
        
    
    
