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
import CMIP5_tools as cmip5

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


def bin_cloud_cover(clisccp):
    if crunchy:
        fobs = cdms.open("/work/marvel1/CLOUD_SEASONS/cloud-seasons/CLOUD_OBS/clt_ISCCP_corrected_198301-200912.nc")
    else:
        fobs = cdms.open("/Users/kmarvel/Google Drive/CLOUD_SEASONS/cloud-seasons/CLOUD_OBS/clt_ISCCP_corrected_198301-200912.nc")
    the_grid = fobs["clt"].getGrid()

    cdutil.setTimeBoundsMonthly(clisccp)
    ann_original_grid = cdutil.YEAR(clisccp)
    ann = ann_original_grid.regrid(the_grid,regridTool='regrid2')
    low = MV.sum(ann(level=(1000*100,680.*100)),axis=1)
    mid = MV.sum(ann(level=(680*100,440.*100)),axis=1)
    high= MV.sum(ann(level=(440*100,0.*100)),axis=1)
    low.id = "low_cloud"
    mid.id="mid_cloud"
    high.id="high_cloud"
    fobs.close()
    return low,mid,high

def abrupt_cloud_cover():
    path = "/kate/cl_regrid_isccp/abrupt4xCO2/"
    files = glob.glob(path+"*r1i**xml")
    afiles=sorted(cmip5.only_most_recent(files))
    if crunchy:
        fobs = cdms.open("/work/marvel1/CLOUD_SEASONS/cloud-seasons/CLOUD_OBS/clt_ISCCP_corrected_198301-200912.nc")
    else:
        fobs = cdms.open("/Users/kmarvel/Google Drive/CLOUD_SEASONS/cloud-seasons/CLOUD_OBS/clt_ISCCP_corrected_198301-200912.nc")
    the_grid = fobs["clt"].getGrid()
    L = len(afiles)
    LOW = MV.zeros((L,140)+the_grid.shape)+1.e20
    MID = MV.zeros((L,140)+the_grid.shape)+1.e20
    HIGH = MV.zeros((L,140)+the_grid.shape)+1.e20
    fobs.close()
    i=0
    
    for i in range(L):
        fname = afiles[i]
        try:
            f = cdms.open(fname)
            clisccp =  f("pgrid_cl")
            low,mid,high = bin_cloud_cover(clisccp)
            LOW[i] = low[:140]
            MID[i] = mid[:140]
            HIGH[i] = high[:140]
            f.close()
        except:
            continue
    axesList = [cmip5.make_model_axis(afiles),low[:140].getTime(),low.getLatitude(),low.getLongitude()]
    
    LOW = MV.masked_where(LOW>1.e10,LOW)
    LOW.id = "low_cloud"
    LOW.setAxisList(axesList)

    MID = MV.masked_where(MID>1.e10,MID)
    MID.id = "mid_cloud"
    MID.setAxisList(axesList)

    HIGH = MV.masked_where(HIGH>1.e10,HIGH)
    HIGH.id = "high_cloud"
    HIGH.setAxisList(axesList)
        
    fw = cdms.open('/work/marvel1/PATTERN_EFFECT/pattern_effect/REGRIDDED_CLOUD/cmip5.ensemble.abrupt4xCO2.r1i1pALL.ann.atm.Amon.cloudbins.ver-1.latestX.nc','w')
    fw.write(LOW)
    fw.write(MID)
    fw.write(HIGH)
    fw.history="Generated using function abrupt_cloud_cover in generate_lcc_from_regridded.py by Kate Marvel 8/27/17"
    fw.close()

def make_xmls(rn):
    path =  "/kate/cl_regrid_isccp/"+rn+"/"
    prefixes = np.unique([x.split("latestX")[0] for x in glob.glob(path+"*")])
    
    for prefix in prefixes:
        
        xml = prefix+"latestX.xml"
        
        cmd = "cdscan -x "+xml+" "+prefix+"*"
        os.system(cmd)
    
def historical_cloud_cover(rn):
    path =  "/kate/cl_regrid_isccp/"+rn+"/"
    files = glob.glob(path+"*xml")
    afiles=sorted(cmip5.only_most_recent(files))
    if crunchy:
        fobs = cdms.open("/work/marvel1/CLOUD_SEASONS/cloud-seasons/CLOUD_OBS/clt_ISCCP_corrected_198301-200912.nc")
    else:
        fobs = cdms.open("/Users/kmarvel/Google Drive/CLOUD_SEASONS/cloud-seasons/CLOUD_OBS/clt_ISCCP_corrected_198301-200912.nc")
    the_grid = fobs["clt"].getGrid()
    L = len(afiles)
    LOW = MV.zeros((L,27)+the_grid.shape)+1.e20
    MID = MV.zeros((L,127)+the_grid.shape)+1.e20
    HIGH = MV.zeros((L,27)+the_grid.shape)+1.e20
    fobs.close()
    
    
    for i in range(L):
        fname = afiles[i]
        try:
            f = cdms.open(fname)
            clisccp =  f("pgrid_cl")(time=('1979-1-1','2005-12-31'))
            low,mid,high = bin_cloud_cover(clisccp)
            LOW[i] = low
            MID[i] = mid
            HIGH[i] = high
            f.close()
        except:
            continue
    axesList = [cmip5.make_model_axis(afiles),low.getTime(),low.getLatitude(),low.getLongitude()]
    
    LOW = MV.masked_where(LOW>1.e10,LOW)
    LOW.id = "low_cloud"
    LOW.setAxisList(axesList)

    MID = MV.masked_where(MID>1.e10,MID)
    MID.id = "mid_cloud"
    MID.setAxisList(axesList)

    HIGH = MV.masked_where(HIGH>1.e10,HIGH)
    HIGH.id = "high_cloud"
    HIGH.setAxisList(axesList)
        
    fw = cdms.open('/work/marvel1/PATTERN_EFFECT/pattern_effect/REGRIDDED_CLOUD/cmip5.ensemble.'+rn+'.r1i1pALL.ann.atm.Amon.cloudbins.ver-1.latestX.nc','w')
    fw.write(LOW)
    fw.write(MID)
    fw.write(HIGH)
    fw.history="Generated using function historical_cloud_cover in generate_lcc_from_regridded.py by Kate Marvel 8/27/17"
    fw.close()
    
