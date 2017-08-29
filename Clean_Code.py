import cdms2 as cdms
import genutil,cdutil
import numpy as np
import MV2 as MV
import matplotlib.pyplot as plt
import scipy.stats as stats
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap
import string
#Modules I wrote
import CMIP5_tools as cmip5
import Plotting


### TO DO 8/18:


#### File Management Functions ####
def match_ensemble_members(X,Y):
    """
    Given two arrays with model axes, find ensemble members in common
    Returns truncated X and Y arrays corresponding to common members only
    """
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

def models_in_common(X,Y):
    """
    X and Y should be ensemble averages
    Return arrays with models in common
    """
    Xmodels = eval(X.getAxis(0).models)
    
     #Hadley Centre and CanEsm Models
    if len(np.where(np.array([x.find("Can")==0 for x in Xmodels]))[0])>0:
        cani=int(np.where(np.array([x.find("Can")==0 for x in Xmodels]))[0])
        
        Xmodels[cani] = "Can*"
    if len(np.where(np.array([x.find("HadGEM2-A")==0 for x in Xmodels]))[0])>0:
        cani=int(np.where(np.array([x.find("HadGEM2-A")==0 for x in Xmodels]))[0])
        Xmodels[cani] = "HadGEM2-A*"
    Ymodels = eval(Y.getAxis(0).models)
    if len(np.where(np.array([x.find("Can")==0 for x in Ymodels]))[0])>0:
        cani=int(np.where(np.array([x.find("Can")==0 for x in Ymodels]))[0])
        Ymodels[cani] = "Can*"
    if len(np.where(np.array([x.find("HadGEM2-A")==0 for x in Ymodels]))[0])>0:
        cani=int(np.where(np.array([x.find("HadGEM2-A")==0 for x in Ymodels]))[0])
        Ymodels[cani] = "HadGEM2-A*"
    common = np.intersect1d(Xmodels,Ymodels)
   
    
    L = len(common)
    Xnew = MV.zeros((L,)+X.shape[1:])
    Ynew = MV.zeros((L,)+Y.shape[1:])
    counter = 0
    for model in common:
        xi = Xmodels.index(model)
        Xnew[counter]=X[xi]

        yi = Ymodels.index(model)
        Ynew[counter]=Y[yi]
        counter+=1
    newmodax = cmip5.make_model_axis(common.tolist())
    Xnew.setAxis(0,newmodax)
    Xnew.id = X.id
    for i in range(len(X.shape))[1:]:
        Xnew.setAxis(i,X.getAxis(i))
    Ynew.setAxis(0,newmodax)
    Ynew.id = Y.id
    for i in range(len(Y.shape))[1:]:
        Ynew.setAxis(i,Y.getAxis(i))
    return Xnew,Ynew

### Raw ingredients: calculate TOA energy balance, temperature change, and forcint ###
    
def Q(typ):
    """
    Get the global mean TOA radiative imbalance from downward SW - (upward SW + upward LW)
    inputs: typ must be one of ["amip","amip4K","amipFuture","historical"]
    
    """
    #Downwelling shortwave
    variable = "rsdt"
    fwrite = cdms.open("DATA/TIMESERIES/cmip5."+typ+"."+variable+".nc")
    RSDT = fwrite(variable)
    rsdt_models = eval(fwrite(variable).getAxis(0).models)
    rsdt_trunc = [x.split("/")[-1].split(variable)[0] for x in rsdt_models]
    fwrite.close()

    #Upwelling shortwave
    variable = "rsut"
    fwrite = cdms.open("DATA/TIMESERIES/cmip5."+typ+"."+variable+".nc")
    rsut_models = eval(fwrite(variable).getAxis(0).models)
    rsut_trunc = [x.split("/")[-1].split(variable)[0] for x in rsut_models]
    RSUT = fwrite(variable)
    fwrite.close()

    #Upwelling longwave
    variable = "rlut"
    fwrite = cdms.open("DATA/TIMESERIES/cmip5."+typ+"."+variable+".nc")
    rlut_models = eval(fwrite(variable).getAxis(0).models)
    rlut_trunc = [x.split("/")[-1].split(variable)[0] for x in rlut_models]
    RLUT = fwrite(variable)
    fwrite.close()

    #Data axis 0 specifies CMIP5 models.  Get the models for which we have all three required components.
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
    #Set the axes (models, time)
    modax = cdms.createAxis(range(len(ok_models)))
    modax.models = str(ok_models.tolist())
    toa.setAxis(0,modax)
    toa.setAxis(1,RLUT.getTime())
    return toa

def TS(typ,variable = "ts"):
    """
    Get the global mean surface skin temperature
    inputs: typ must be one of ["amip","amip4K","amipFuture","historical"]
    
    """
    #variable = "ts"
    fwrite = cdms.open("DATA/TIMESERIES/cmip5."+typ+"."+variable+".nc")
    ts_models = eval(fwrite(variable).getAxis(0).models)
    ts_trunc = [x.split("/")[-1].split(variable)[0] for x in ts_models]
    temperature = fwrite(variable)
    toa = Q(typ)
    models = np.intersect1d( ts_trunc, eval(toa.getAxis(0).models))
    ts = MV.zeros((len(models),temperature.shape[1]))
    j=0
    for model in models:
        i = ts_trunc.index(model)
        ts[j] = temperature[i]
        j+=1
    modax = cdms.createAxis(range(len(models)))
    modax.models = str(models.tolist())
    ts.setAxis(0,modax)
    ts.setAxis(1,temperature.getTime())
    return ts

def get_forcing(forcing_from = "AR5"):
    """
    Return the global mean radiative forcing timeseries
    Inputs:forcing_from = [AR5, giss]
    
    """
    if forcing_from == "giss":
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
        # Get AR5 forcing data from IPCC AR5
        far5 = cdms.open("AR5_forcings.nc","r")
        forcing = far5("total")
        forcing.id="forcing"
    return forcing

def estimate_ECS(typ,temperature_variable="ts"):
    """
    Estimate the ECS over the AMIP time period (1979-2005 inclusive) as the regression coefficient of annual means.
    """
    if typ == "equil":
        return cmip5.all_clim_sens()
    
    T = cdutil.YEAR(TS(typ,variable=temperature_variable)) #Annual mean surface temperature
    QTOA = cdutil.YEAR(Q(typ)) #Annual mean TOA radiative imbalance
    F = get_forcing(forcing_from="AR5")(time=('1979-1-1','2005-12-31')) #Assume AR5 forcing
    
    # Differences with respect to the first decade
    dT = cmip5.cdms_clone(T.asma() - MV.average(T[:,:10],axis=1).asma()[:,np.newaxis],T)
    dF = F - MV.average(F[:10])
    dQ =  cmip5.cdms_clone(QTOA.asma() - MV.average(QTOA[:,:10],axis=1).asma()[:,np.newaxis],QTOA)

    #Assume standard value for CO2 doubling
    F2 = 3.7
    
    #Regress over annual means
    y = dF-dQ
    x = dT
        
    lam = np.ma.zeros(x.shape[0])+1.e20
    for i in range(lam.shape[0]):
        try:
            L = np.polyfit(x[i],y[i],1)[0] #Fit a polynomial with degree 1
            lam[i] = F2/L
        except:
            continue #If data is missing, do nothing and mask below
    ECS = MV.array(np.ma.masked_where(np.abs(lam)>1.e10,lam))
    ECS.setAxis(0,x.getAxis(0))
    return ECS


####### FIGURE CODE #########

def Figure1(cmap=cm.magma):
    """
    Figure 1 of proposed paper
    """
    
    historical = estimate_ECS("historical")
    amip = estimate_ECS("amip")
    equil = estimate_ECS("equil")
    #Ensemble averages
    amip_ensav = cmip5.ensemble2multimodel(amip)
    historical_ensav = cmip5.ensemble2multimodel(historical)
    #Make sure we're looking at the same set of models
    amip2,hist2=models_in_common(amip_ensav,historical_ensav)
    amip_comp,eq_comp = models_in_common(amip2,equil)
    historical_comp,eq_comp_h = models_in_common(hist2,equil)
    
    #First panel: histograms of inferred ECS
    plt.subplot(221)
    x_fit = np.linspace(0,7)
    plt.hist(historical.compressed(),20,color=cmap(.3),ec=cmap(.3),alpha=.5,normed=True,label="Historical")
    scatter,loc,mean = stats.lognorm.fit(historical.compressed())
    pdf_fitted = stats.lognorm.pdf(x_fit,scatter,loc,mean)
    plt.plot(x_fit,pdf_fitted,color=cmap(.3),lw=3)
        
    plt.hist(amip.compressed(),20,color=cmap(.6),ec=cmap(0.6),alpha=.5,normed=True,label="AMIP")
    scatter,loc,mean = stats.lognorm.fit(amip.compressed())
    pdf_fitted = stats.lognorm.pdf(x_fit,scatter,loc,mean)
    plt.plot(x_fit,pdf_fitted,color=cmap(.6),lw=3)
        
    plt.hist(equil,color=cmap(.9),ec=cmap(.9),alpha=.5,normed=True,label="Forster et al.")
    scatter,loc,mean = stats.lognorm.fit(cmip5.all_clim_sens())
    pdf_fitted = stats.lognorm.pdf(x_fit,scatter,loc,mean)
    plt.plot(x_fit,pdf_fitted,color=cmap(.9),lw=3)
    plt.xlabel(r'Estimated ECS ($^{\circ}$C)')
    plt.ylabel("PDF")
    plt.legend(loc=0)
    plt.title("(a): Inferred ECS")

    #Second panel: correlate historical and amip
    plt.subplot(222)
    Plotting.scatterplot_cmip(amip_comp,historical_comp)
    plt.xlabel(r"Ensemble mean ECS inferred from AMIP ($^{\circ}$)")
    plt.ylabel(r"Ensemble mean ECS inferred from historical ($^{\circ}$)")
    x = np.linspace(1,4.5)
    plt.plot(x,x,"k--")
    plt.legend(loc=0,numpoints=1,ncol=2,fontsize=10)
    corrcoef=str(np.round(float(genutil.statistics.correlation(amip_comp,historical_comp)),2))
    plt.title("(b) AMIP and historical: R = "+corrcoef)

    #Third panel: correlate historical and equilibrium
    plt.subplot(223)
    Plotting.scatterplot_cmip(historical_comp,eq_comp)
    plt.ylabel(r"Equilibrium ECS ($^{\circ}$)")
    plt.xlabel(r"Ensemble mean ECS inferred from historical ($^{\circ}$)")
    x = np.linspace(1,4.5)
    plt.plot(x,x,"k--")
    plt.legend(loc=0,numpoints=1,ncol=2,fontsize=10)
    corrcoef=str(np.round(float(genutil.statistics.correlation(eq_comp,historical_comp)),2))
    plt.title("(c) Historical and equilibrium: R = "+corrcoef)

    
    #Fourth panel: correlate amip and equilibrium
    plt.subplot(224)
    Plotting.scatterplot_cmip(amip_comp,eq_comp)
    plt.ylabel(r"Equilibrium ECS ($^{\circ}$)")
    plt.xlabel(r"Ensemble mean ECS inferred from AMIP ($^{\circ}$)")
    x = np.linspace(1,4.5)
    plt.plot(x,x,"k--")
    plt.legend(loc=0,numpoints=1,ncol=2,fontsize=10)
    corrcoef=str(np.round(float(genutil.statistics.correlation(eq_comp,amip_comp)),2))
    plt.title("(d) AMIP and equilibrium: R = "+corrcoef)


    
def Figure2(cmap = cm.Dark2):
    """
    Figure 2 of proposed paper.  Shows the inferred ensemble spread for each CMIP5 model in historical vs amip
    """
    #Estimate the ECS using linear regression over annual means
    historical = estimate_ECS("historical")
    amip = estimate_ECS("amip")
    
    #Sort ensemble members into dictionaries labeled by model
    hdict = cmip5.ensemble_dictionary(historical)
    adict = cmip5.ensemble_dictionary(amip)
    i = 0
    models = sorted(np.append(np.intersect1d(hdict.keys(),adict.keys()),(np.array(['HadGEM2-AO','CanESM2']))))

  
    for model in models:
        #Hadley Centre and Canadian models have different names for fully coupled and atmosphere-only
        if model.find("HadGEM")>=0:
            amodel = 'HadGEM2-A'
            hmodel = 'HadGEM2-AO'
        elif model.find("Can")>=0:
            amodel = 'CanAM4'
            hmodel = 'CanESM2'
        else:
            amodel = model
            hmodel = model
        #Get the estimated ECS corresponding to the model
        A = amip.asma()[adict[amodel]]
        H = historical.asma()[hdict[hmodel]]
        #Now plot
        x = np.zeros_like(A.compressed())+i+.33
        plt.plot(x,A.compressed(),"o",color=cmap(.6),mec=cmap(.7),markersize=5)
        x = np.zeros_like(H.compressed())+i+.66
        plt.plot(x,H.compressed(),"s",color=cmap(.3),mec=cmap(.4),markersize=5)
        plt.axvline(i+1,c="k",ls=":")
       
        i+=1
    plt.xticks(np.arange(len(models))+.5,models,rotation=90)
    plt.xlim(0,len(models))



#### CLOUDS AND ECS #####             
def tropical_marine_lcc(X):
    """
    Average the variable X over tropical marine stratocumulus regions
    """
    #Stratus regions as defined by Qu et al 2015
 
    Peru = cdutil.region.domain(latitude=(-30,-10),longitude=(-110,-70))
    Namibia = cdutil.region.domain(latitude=(-30,-10),longitude=(-25,15))
    California = cdutil.region.domain(latitude=(15,35),longitude=(-155,-115))
    Australia = cdutil.region.domain(latitude=(-35,-25),longitude=(75,115))
    Canaries = cdutil.region.domain(latitude=(10,30),longitude=(-50,-10))
    
    return cdutil.averager(X(Peru),axis='xy')+cdutil.averager(X(Namibia),axis='xy')+cdutil.averager(X(California),axis='xy')+cdutil.averager(X(Australia),axis='xy')+cdutil.averager(X(Canaries),axis='xy')




def compare_ECS_and_tropicalcloud(experiment,compare_to="LCC",ensemble_average=True,plot=True):
    """
    On a model-by-model basis, compare low cloud cover or SWCRE to ecs
    compare_to should be LCC or SWCRE
    ensemble_average = True to average over ensemble members
    """
    #Get the estimated ECS 
    ecs = estimate_ECS(experiment)
    #Get either LCC or SWCRE to compare to
    f = cdms.open("MAPS/"+string.upper(experiment)+"_"+compare_to+".nc")
    if compare_to is "LCC":
        lcc = f("lcc")
    else:
        lcc = f("SWCRE")
    stratocumulus = tropical_marine_lcc(lcc)
    
    
    if ensemble_average:
        ecs_ensav = cmip5.ensemble2multimodel(ecs)
        stratocum_ensav = cmip5.ensemble2multimodel(stratocumulus)
        X,Y = models_in_common(ecs_ensav,stratocum_ensav)
    else:
        X,Y = match_ensemble_members(ecs,stratocumulus)
    if plot:
        Plotting.scatterplot_cmip(X,Y)
        plt.title("R = "+str(np.round(float(genutil.statistics.correlation(X,Y)),2)))
        plt.legend(loc=0,numpoints=1,ncol=1,fontsize=6)
        plt.xlabel(r"ECS ($^{\circ}$) estimated from "+experiment)
        plt.ylabel(r'$\Delta$ '+compare_to)
    return X,Y

def Figure3():
    """
    Plot low cloud cover vs amip
    """
    f = cdms.open("MAPS/AMIP_LCC.nc")
    lcc = f("lcc")
    X = MV.average(lcc,axis=0)
    plt.subplot(211)
    plot_tropical_lcc(X,projection="cyl")
    plt.title("(a): Tropical Marine Cloud Regions")
    ax = plt.subplot(212)
    x,y=compare_ECS_and_tropicalcloud("amip",ensemble_average=False,compare_to="LCC")
    ax.set_title("(b) "+ax.get_title())
    #plt.title("(b): AMIP ECS vs tropical LCC")
    


####TEMPORARY EXPERIMENTAL CODE #####
def plot_tropical_lcc(X,projection='moll'):
    """
    Plot stratocumulus regions off the western coasts of continents
    """
    if X.id == "SWCRE":
        v=8
        label=r'$\Delta$ SWCRE (Wm$^2$)'
    else:
        v=1.5
        label=r'$\Delta$ LCC (%)'
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
    plt.colorbar(orientation='horizontal',label=label)

def mask_sea_ice(X):
    """
    Mask areas covered by sea ice at some point (for AMIP)
    """
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
    """
    use observational SST land mask 
    """
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


def sst_patterns(experiment,land=False,sea_ice=False):
    variable = "ts"#REPLACE WITH TS!!!!!!!!!
    f = cdms.open("PATTERN_MAPS/"+experiment+"."+variable+".nc") #REPLACE WITH TS!!!!!!!!!
    tspatt = f(variable)
    f.close()
    if land is False:
        tspatt = mask_land(tspatt)
    if sea_ice is False:
        tspatt = mask_sea_ice(tspatt)
    return tspatt

def plot_single_model(X,Y,model):
    Xdict = cmip5.ensemble_dictionary(X)
    Ydict = cmip5.ensemble_dictionary(Y)
    xi = Xdict[model]
    yi = Ydict[model]
    Xplot = X.asma()[xi]
    Yplot = Y.asma()[yi]
    m=Plotting.model_dictionary()[model]["marker"]
    c=Plotting.model_dictionary()[model]["color"]
    plt.plot(Xplot,Yplot,m,color=c)

def highest_sensitivity_patterns(model):
    historicalsst=sst_patterns("historical")
    dict = cmip5.ensemble_dictionary(historicalsst)
    historical = estimate_ECS("historical")
    xi=dict[model]
    imax=np.argmax(historical.asma()[xi])
    print cmip5.models(historicalsst)[xi[imax]]
    return historicalsst[xi[imax]]

def average_high_sens_patterns():
    historicalsst=sst_patterns("historical",land=True,sea_ice=True)
    dict = cmip5.ensemble_dictionary(historicalsst)
    models = sorted(dict.keys())
    historical = estimate_ECS("historical")
    H = MV.zeros((len(models),)+historicalsst.shape[1:])+1.e20
    for i in range(len(models)):
        model=models[i]
        xi=dict[model]
        if len(xi)>3:
            print i
            imax=np.argmax(historical.asma()[xi])
            H[i]=historicalsst[xi[imax]]
            
    H = MV.masked_where(H>1.e10,H)
    #Ha = MV.average(H,axis=0)
    H.setAxis(0,cmip5.make_model_axis(models))
    H.setAxis(1,historicalsst.getLatitude())
    H.setAxis(2,historicalsst.getLongitude())
    return H

def average_low_sens_patterns():
    historicalsst=sst_patterns("historical",land=True,sea_ice=True)
    dict = cmip5.ensemble_dictionary(historicalsst)
    models = sorted(dict.keys())
    historical = estimate_ECS("historical")
    H = MV.zeros((len(models),)+historicalsst.shape[1:])+1.e20
    for i in range(len(models)):
        model=models[i]
        xi=dict[model]
        if len(xi)>3:
            print i
            imin=np.argmin(historical.asma()[xi])
            H[i]=historicalsst[xi[imin]]
            
    H = MV.masked_where(H>1.e10,H)
    #Ha = MV.average(H,axis=0)
    H.setAxis(0,cmip5.make_model_axis(models))
    H.setAxis(1,historicalsst.getLatitude())
    H.setAxis(2,historicalsst.getLongitude())
    return H
    

def lowest_sensitivity_patterns(model):
    historicalsst=sst_patterns("historical")
    dict = cmip5.ensemble_dictionary(historicalsst)
    historical = estimate_ECS("historical")
    xi=dict[model]
    imin=np.argmin(historical.asma()[xi])
    print cmip5.models(historicalsst)[xi[imin]]
    return historicalsst[xi[imin]]
        
def southern_ocean_lcc(X):
    return cdutil.averager(X(latitude=(-90,-50)),axis='xy')
        
def compare_ECS_and_southern_cloud(experiment,compare_to="LCC",ensemble_average=True):
    """
    On a model-by-model basis, compare low cloud cover or SWCRE to ecs
    compare_to should be LCC or SWCRE
    ensemble_average = True to average over ensemble members
    """
    #Get the estimated ECS 
    ecs = estimate_ECS(experiment)
    #Get either LCC or SWCRE to compare to
    f = cdms.open("MAPS/"+string.upper(experiment)+"_"+compare_to+".nc")
    if compare_to is "LCC":
        lcc = f("lcc")
    else:
        lcc = f("SWCRE")
    stratocumulus = southern_ocean_lcc(lcc)
    
    
    if ensemble_average:
        ecs_ensav = cmip5.ensemble2multimodel(ecs)
        stratocum_ensav = cmip5.ensemble2multimodel(stratocumulus)
        X,Y = models_in_common(ecs_ensav,stratocum_ensav)
    else:
        X,Y = match_ensemble_members(ecs,stratocumulus)
    return X,Y
    
    
####### DEPRECATED ########
## Use Plotting.scatterplot_cmip instead
def scatterplot_stuff(X,Y,cmap=cm.viridis,ensemble_average=True):
    """
    Scatterplot the arrays X and Y.  If ensemble_average, just plot.  Otherwise, group ensemble members by model.
    """
    markers = np.tile(["o","s","D","*","p","v"],100)
    if not ensemble_average:
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

