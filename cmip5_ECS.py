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


def models_in_common(dict):

    keys = dict.keys()
    typ = dict[keys[0]][0].split(".")[2]
    toa = Q(typ)

    ok_models = np.intersect1d(eval(toa.getAxis(0).models),dict[keys[0]])
    for key in keys[1:]:
        ok_models = np.intersect1d(ok_models,dict[key])
    return ok_models
    

def OcnQFlux(typ):
    """
    Get the surface radiative imbalance from (downward SW+downward LW) - (upward SW + upward LW+ sensible heat + latent heat)
    inputs: typ must be one of ["amip","amip4K","amipFuture","historical"]
    
    """

    surface = {"hfls": "Surface Upward Latent Heat Flux",\
                "hfss": "Surface Upward Sensible Heat Flux",\
                "rlds":"Surface Downwelling Longwave Radiation",\
                "rlus":"Surface Upwelling Longwave Radiation",\
                "rsds":"Surface Downwelling Shortwave Radiation",\
                "rsus": "Surface Upwelling Shortwave Radiation"}

    upward = {"hfls": "Surface Upward Latent Heat Flux",\
                "hfss": "Surface Upward Sensible Heat Flux",\
                "rlus":"Surface Upwelling Longwave Radiation",\
                "rsus": "Surface Upwelling Shortwave Radiation"}
                
    downward = {"rlds":"Surface Downwelling Longwave Radiation",\
                "rsds":"Surface Downwelling Shortwave Radiation"}
                

                
    TRUNC = {}
    DATA = {}
    for variable in surface.keys():

        fwrite = cdms.open("DATA/TIMESERIES/cmip5."+typ+"."+variable+".nc")
   
        DATA[variable] = fwrite(variable)
        all_models = eval(fwrite(variable).getAxis(0).models)
        TRUNC[variable] = [x.split("/")[-1].split(variable)[0] for x in all_models]
        fwrite.close()
   
    ok_models = models_in_common(TRUNC)

    qflux = MV.zeros((len(ok_models),DATA[variable].shape[1]))
    i=0
    for model in ok_models:
        test = MV.zeros(DATA["hfls"].shape[1])
        for variable in downward.keys():
            test=test + DATA[variable][TRUNC[variable].index(model)]
        for variable in upward.keys():
            test = test - DATA[variable][TRUNC[variable].index(model)]
        qflux[i] = test
        i+=1
    modax = cdms.createAxis(range(len(ok_models)))
    modax.models = str(ok_models.tolist())
    qflux.setAxis(0,modax)
    qflux.setAxis(1,DATA["hfls"].getTime())
    return qflux

    
    

   
    # ok_models = np.intersect1d(np.intersect1d(rlut_trunc,rsut_trunc),rsdt_trunc)

    # toa = MV.zeros((len(ok_models),RSDT.shape[1]))
    # i=0
    # for model in ok_models:
    #     rsdt_mod = RSDT[rsdt_trunc.index(model)]
    #     rsut_mod = RSUT[rsut_trunc.index(model)]
    #     rlut_mod = RLUT[rlut_trunc.index(model)]
    #     test = rsdt_mod - (rlut_mod+rsut_mod)
    #     toa[i] = test
    #     i+=1
    # modax = cdms.createAxis(range(len(ok_models)))
    # modax.models = str(ok_models.tolist())
    # toa.setAxis(0,modax)
    # toa.setAxis(1,RLUT.getTime())
    # return toa
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

def TS(typ):
    variable = "ts"
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
        sw[j] = -(RSUTCS[csi] - RSUT[alli])
        j+=1
    modax = cdms.createAxis(range(len(models)))
    modax.models = str(models.tolist())
    sw.setAxis(0,modax)
    sw.setAxis(1,RSUT.getTime())

    return sw
def get_forcing(typ):
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
        self.SWCRE = cdutil.YEAR(SWCRE(typ))
        self.surf = cdutil.YEAR(OcnQFlux(typ))
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
        self.dSWCRE =  cmip5.cdms_clone(self.SWCRE.asma() - MV.average(self.SWCRE[:,:10],axis=1).asma()[:,np.newaxis],self.SWCRE)
        self.dsurf =  cmip5.cdms_clone(self.surf.asma() - MV.average(self.surf[:,:10],axis=1).asma()[:,np.newaxis],self.surf)
    def estimate_ECS(self,typ = "decades"):
        F2 = 3.7
        if typ.find("regression")>=0:
            if typ.find("decadal")>=0:
                y = genutil.filters.runningaverage(self.dF-self.dQ,10,axis=1)
                x = genutil.filters.runningaverage(self.dT,10,axis=1)
            else:
                y = self.dF-self.dQ
                x = self.dT
        
            lam = np.ma.zeros(x.shape[0])+1.e20
            for i in range(lam.shape[0]):
                try:
                    L = np.polyfit(x[i],y[i],1)[0]
                    lam[i] = F2/L
                except:
                    continue
        else:
            FmQ = self.F-self.Q
            x = MV.average(FmQ[:,-10:],axis=1)- MV.average(FmQ[:,:10],axis=1)
            y = MV.average(self.T[:,-10:],axis=1)-MV.average(self.T[:,:10],axis=1)
            lam = F2*y/x
        ECS = MV.array(np.ma.masked_where(np.abs(lam)>1.e10,lam))
        ECS.setAxis(0,x.getAxis(0))
        self.ECS = ECS
            
        


    
                
def correspond_models(model,typ = "decades"):
    historical = ECS("historical")
    historical.estimate_ecs(typ=typ)
    amip = ECS("amip")
    amip.estimate_ecs(typ=typ)
     
    hmodels = np.array(eval(historical.ECS.getAxis(0).models))
    amodels = np.array(eval(amip.ECS.getAxis(0).models))
    
    hist_i=np.where(np.array([x.split(".")[1]==model for x in hmodels]))[0]
    amip_i=np.where(np.array([x.split(".")[1]==model for x in amodels]))[0]

    return hist_i,amip_i


class Comparison():
    def __init__(self,typ = "decades"):
        self.historical = ECS("historical")
        self.historical.estimate_ECS(typ=typ)
        self.amip = ECS("amip")
        self.amip.estimate_ECS(typ=typ)

        hmodels = np.array(eval(self.historical.ECS.getAxis(0).models))
        amodels = np.array(eval(self.amip.ECS.getAxis(0).models))
        self.unique_hmodels = []
        for x in hmodels:
            if x.find("GISS-E2-R")>=0: #Split GISS p1 and p3 runs 
                self.unique_hmodels += [x.split(".")[1]+"."+x.split(".")[3].split("i1")[-1]]
            else:
                self.unique_hmodels += [x.split(".")[1]]
        self.unique_hmodels = np.unique(np.array(self.unique_hmodels))


        self.unique_amodels = []
        for x in amodels:
            if x.find("GISS-E2-R")>=0:
                self.unique_amodels += [x.split(".")[1]+"."+x.split(".")[3].split("i1")[-1]]
            else:
                self.unique_amodels += [x.split(".")[1]]
        self.unique_amodels = np.unique(np.array(self.unique_amodels))
        
        self.models = np.intersect1d(self.unique_hmodels,self.unique_amodels)
        self.models = np.array(sorted(np.append(['HadGEM2-A',"CanAM4"],self.models))) #Add HadGEM back in
    
    def correspond_models(self,model,typ="ECS"):
        if typ.find("zonal")<0:
            hmodels = np.array(eval(getattr(self.historical,typ).getAxis(0).models))
            amodels = np.array(eval(getattr(self.amip,typ).getAxis(0).models))
        else:
            variable = typ.split("_zonal")[0]
            fa = cdms.open("DATA/TIMESERIES/ZONAL/cmip5.amip."+variable+".nc")
            amodels = eval(fa[variable].getAxis(0).models)
            fh = cdms.open("DATA/TIMESERIES/ZONAL/cmip5.historical."+variable+".nc")
            hmodels = eval(fh[variable].getAxis(0).models)
            fa.close()
            fh.close()
        if model.find("HadGEM")>=0:
            amodel = 'HadGEM2-A'
            hmodel = 'HadGEM2-AO'
        elif model.find("Can")>=0:
            amodel = 'CanAM4'
            hmodel = 'CanESM2'
        else:
            amodel = model
            hmodel = model
        if model.find("GISS")<0:
            hist_i=np.where(np.array([x.split(".")[1]==hmodel for x in hmodels]))[0]
            amip_i=np.where(np.array([x.split(".")[1]==amodel for x in amodels]))[0]
        else:
            hist_i=np.where(np.array([x.split(".")[1]+"."+x.split(".")[3].split("i1")[-1]==hmodel for x in hmodels]))[0]
            amip_i=np.where(np.array([x.split(".")[1]+"."+x.split(".")[3].split("i1")[-1]==amodel for x in amodels]))[0]
            

        return hist_i,amip_i

    def ensemble_compare_zonal(self,model,variable):
        fa = cdms.open("DATA/TIMESERIES/ZONAL/cmip5.amip."+variable+".nc")
        fh = cdms.open("DATA/TIMESERIES/ZONAL/cmip5.historical."+variable+".nc")
        hist_i,amip_i = self.correspond_models(model,variable+"_zonal")
        amip_variable = fa(variable)
        A = MV.array(amip_variable.asma()[amip_i])
        A.setAxis(1,amip_variable.getTime())
        A.setAxis(2,amip_variable.getLatitude())

        historical_variable = fh(variable)
        H = MV.array(historical_variable.asma()[hist_i])
        H.setAxis(1,historical_variable.getTime())
        H.setAxis(2,historical_variable.getLatitude())
        return A,H
    def ensemble_compare_ECS(self,model):
        hist_i,amip_i = self.correspond_models(model)
        A = self.amip.ECS.asma()[amip_i]
        H = self.historical.ECS.asma()[hist_i]
        return A,H

    def ensemble_compare_Q(self,model):
        historical_i,amip_i = self.correspond_models(model,"Q")
        amip_lastdecade = self.amip.Q[:,-10:]
        amip_firstdecade = self.amip.Q[:,:10]
        amip = MV.average(amip_lastdecade,axis=1)-MV.average(amip_firstdecade, axis=1)
        A = amip.asma()[amip_i]

        historical_lastdecade = self.historical.Q[:,-10:]
        historical_firstdecade = self.historical.Q[:,:10]
        historical = MV.average(historical_lastdecade,axis=1)-MV.average(historical_firstdecade, axis=1)
        H = historical.asma()[historical_i]
        
        return A,H

    def ensemble_compare_SWCRE(self,model):
        historical_i,amip_i = self.correspond_models(model,"SWCRE")
        amip_lastdecade = self.amip.SWCRE[:,-10:]
        amip_firstdecade = self.amip.SWCRE[:,:10]
        amip = MV.average(amip_lastdecade,axis=1)-MV.average(amip_firstdecade, axis=1)
        A = amip.asma()[amip_i]

        historical_lastdecade = self.historical.SWCRE[:,-10:]
        historical_firstdecade = self.historical.SWCRE[:,:10]
        historical = MV.average(historical_lastdecade,axis=1)-MV.average(historical_firstdecade, axis=1)
        H = historical.asma()[historical_i]
        
        return A,H

    def ensemble_compare_T(self,model):
        historical_i,amip_i = self.correspond_models(model)
        amip_lastdecade = self.amip.T[:,-10:]
        amip_firstdecade = self.amip.T[:,:10]
        amip = MV.average(amip_lastdecade,axis=1)-MV.average(amip_firstdecade, axis=1)
        A = amip.asma()[amip_i]

        historical_lastdecade = self.historical.T[:,-10:]
        historical_firstdecade = self.historical.T[:,:10]
        historical = MV.average(historical_lastdecade,axis=1)-MV.average(historical_firstdecade, axis=1)
        H = historical.asma()[historical_i]
        
        return A,H

    def ECS_histograms(self,cmap=cm.Dark2):
        """ Distributions of ECS estimated from AMIP, historical, and Forster et al """
        bigA = []
        bigH = []
        for model in self.models:
            A,H = self.ensemble_compare_ECS(model)
            bigA += A.filled().tolist()
            bigH += H.filled().tolist()
        bigA = np.ma.array(bigA)
        bigA = np.ma.masked_where(bigA>1.e10,bigA)
        bigH = np.ma.array(bigH)
        bigH = np.ma.masked_where(bigH>1.e10,bigH)
        x_fit = np.linspace(0,7)
        plt.hist(bigH.compressed(),20,color=cmap(.3),ec=cmap(.3),alpha=.5,normed=True,label="Historical")
        scatter,loc,mean = stats.lognorm.fit(bigH.compressed())
        pdf_fitted = stats.lognorm.pdf(x_fit,scatter,loc,mean)
        plt.plot(x_fit,pdf_fitted,color=cmap(.3),lw=3)
        
        plt.hist(bigA.compressed(),20,color=cmap(.6),ec=cmap(0.6),alpha=.5,normed=True,label="AMIP")
        scatter,loc,mean = stats.lognorm.fit(bigA.compressed())
        pdf_fitted = stats.lognorm.pdf(x_fit,scatter,loc,mean)
        plt.plot(x_fit,pdf_fitted,color=cmap(.6),lw=3)
        
        plt.hist(cmip5.all_clim_sens(),color=cmap(.9),ec=cmap(.9),alpha=.5,normed=True,label="Forster et al.")
        scatter,loc,mean = stats.lognorm.fit(cmip5.all_clim_sens())
        pdf_fitted = stats.lognorm.pdf(x_fit,scatter,loc,mean)
        plt.plot(x_fit,pdf_fitted,color=cmap(.9),lw=3)
        plt.xlabel(r'Estimated ECS ($^{\circ}$C)')
        plt.ylabel("PDF")
        return bigA,bigH

    def historical_vs_amip_figure(self,to_compare,cmap=cm.magma):
        i=0
        func = getattr(self,"ensemble_compare_"+to_compare)
        for model in self.models:
            A,H = func(model)
            x = np.zeros_like(A.compressed())+i+.33
            plt.plot(x,A.compressed(),"o",color=cmap(.6),mec=cmap(.7),markersize=10)
            x = np.zeros_like(H.compressed())+i+.66
            plt.plot(x,H.compressed(),"s",color=cmap(.3),mec=cmap(.4),markersize=10)
            plt.axvline(i+1,c="k",ls=":")
            i+=1
        plt.xticks(np.arange(len(self.models))+.5,self.models,rotation=90)

    def scatterplot_averages(self,to_compare="ECS",cmap=cm.viridis):
        Aecs = []
        Hecs = []
        for model in self.models:
            func = getattr(self,"ensemble_compare_"+to_compare)
            A,H = func(model)
            Aecs+=[np.ma.average(A)]
            Hecs+=[np.ma.average(H)]
        modax = cmip5.make_model_axis(self.models)
        Aecs = MV.array(Aecs)
        Aecs.setAxis(0,modax)

        Hecs = MV.array(Hecs)
        Hecs.setAxis(0,modax)

        markers = np.tile(["o","s","D","*","p","v"],7)
        for i in range(len(self.models)):
            c=cmap(float(i)/float(len(self.models)))
            plt.plot([Aecs[i]],[Hecs[i]], markers[i],markersize=10,color=c,label=self.models[i])
        #plt.xlim(0,5)
        #plt.ylim(0,5)
        plt.xlabel(to_compare+" inferred from AMIP")
        plt.ylabel(to_compare+" inferred from historical")
        plt.legend(loc=0,numpoints=1,ncol=4,fontsize=10)
        x = np.linspace(*plt.xlim())
        plt.plot(x,x,ls=":",c="k")
                    
        

           
            
        

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

sys.path.append("/Users/kmarvel/Google Drive/SEASONAL/")
import peakfinder as pf
                
def play_with_zonal_swcre(C,model,season="YEAR",typ="amip",**kwargs):

    A,H = C.ensemble_compare_zonal(model,"rsut")
    Aclear,Hclear = C.ensemble_compare_zonal(model,"rsutcs")
    Aswcre=A-Aclear
    Hswcre = H-Hclear

    cdutil.setTimeBoundsMonthly(Aswcre)
    cdutil.setTimeBoundsMonthly(Hswcre)

    Aswcre_year = getattr(cdutil,season)(Aswcre)
    Hswcre_year =getattr(cdutil,season)(Hswcre)
    

    
    lat_extrema = MV.zeros(Aswcre_year.shape[:-1]+(5,))
    sw_extrema = MV.zeros(Aswcre_year.shape[:-1]+(5,))

    for i in range(Aswcre_year.shape[0]):
        X = pf.spatially_smooth(Aswcre_year[i])
        x,y = pf.thermo_and_dynamic(X)
        lat_extrema[i] = x
        sw_extrema[i]=y
    print "OK"
    lat_extrema.setAxis(0,Aswcre_year.getAxis(0))
    lat_extrema.setAxis(1,Aswcre_year.getAxis(1))
    lat_extrema.id='lat'
    
    sw_extrema.setAxis(0,Aswcre_year.getAxis(0))
    sw_extrema.setAxis(1,Aswcre_year.getAxis(1))
    sw_extrema.id='swcre'
    lat_plot(MV.average(MV.average(Hswcre_year,axis=0),axis=0),**kwargs)
    clim = cdutil.YEAR.climatology(MV.average(Hswcre_year,axis=0))
    sm = pf.spatially_smooth(clim)
    lat_plot(sm[0],ls="--")
    x,y=pf.thermo_and_dynamic(sm)
    x=x[0]
    y=y[0]
    plt.plot(x.asma(),y.asma(),"ko")
    
    lat_ensav=MV.average(lat_extrema,axis=0)
    lat_trends = cmip5.get_linear_trends(lat_ensav)

    sw_ensav=MV.average(sw_extrema,axis=0)
    sw_trends = cmip5.get_linear_trends(sw_ensav)

    for i in range(5):
        plt.arrow(x[i],y[i],lat_trends[i]*100,0)
        plt.arrow(x[i],y[i],0,sw_trends[i]*100)
        
    
    
    return lat_trends,sw_trends
   
    
def zonal_prettyplot_swcre(C,model,**kwargs):
    A,H = C.ensemble_compare_zonal(model,"rsut")
    Aclear,Hclear = C.ensemble_compare_zonal(model,"rsutcs")
    Aswcre=A-Aclear
    Hswcre = H-Hclear

    cdutil.setTimeBoundsMonthly(Aswcre)
    cdutil.setTimeBoundsMonthly(Hswcre)

    Aswcre_year = cdutil.YEAR(Aswcre)
    Hswcre_year = cdutil.YEAR(Hswcre)
    
    lat_plot(MV.average(MV.average(Hswcre_year,axis=0),axis=0),label=model,**kwargs)

    
def zonal_compare_swcre(C,model):
    A,H = C.ensemble_compare_zonal(model,"rsut")
    Aclear,Hclear = C.ensemble_compare_zonal(model,"rsutcs")
    Aswcre=A-Aclear
    Hswcre = H-Hclear

    cdutil.setTimeBoundsMonthly(Aswcre)
    cdutil.setTimeBoundsMonthly(Hswcre)

    Aswcre_year = cdutil.YEAR(Aswcre)
    Hswcre_year = cdutil.YEAR(Hswcre)

    Adiff=last_ten_minus_first_ten(Aswcre_year)
    Hdiff=last_ten_minus_first_ten(Hswcre_year)

    return Adiff,Hdiff
