import os
cwd = os.getcwd()
os.chdir('/Users/adamamiller/software/git_repos/pandas')
import pandas as pd
os.chdir(cwd)

import glob

import numpy as np
from scipy import optimize as op
from scipy.stats import median_abs_deviation
import matplotlib.pyplot as plt

from astropy.time import Time
from astropy.coordinates import EarthLocation, SkyCoord, AltAz
import astropy.units as u

import gc
import pkg_resources



def read_ipac_fps(fps_file):
    """Read IPAC fps file, manipulate results into pandas DataFrame
    
    The IPAC forced photometry service (fps; 
    https://ztfweb.ipac.caltech.edu/cgi-bin/requestForcedPhotometry.cgi)
    produces fixed-position, PSF-model flux measurements from ZTF 
    subtraction images. This function reads fps results, an IPAC table, 
    calculates the airmass to determine the 'cloudy' parameter, then 
    returns a pandas DataFrame and an array (fcqfid) of the unique 
    field+chip+quadrant+filter combo for each observation.
    
    Parameters
    ----------
    fps_file : str
        File path to an IPAC table text file with fps results
    
    Returns
    -------
    fp_det : (pandas) DataFrame
        Data frame with fps flux measurements and updated bit mask 
    
    """
    
    ipac_version = pd.read_csv(fps_file, skiprows=1, nrows=1).columns[1][2]
    if ipac_version == '2':
        fp = pd.read_csv(fps_file, 
                         delim_whitespace=True, comment='#', skiprows=70,
                         names=['ipac_index', 'field', 'ccdid', 'qid',
                                'filter', 'pid', 'infobitssci', 'airmass', 
                                'moonalt', 'moonillf', 'moonra', 'moondec', 
                                'sciinpseeing', 'scibckgnd', 'scisigpix', 
                                'scimaglim',  'zpmaginpsci', 'zpmaginpsciunc', 
                                'zpmaginpscirms', 'clrcoeff', 'clrcoeffunc', 
                                'ncalmatches',  'exptime', 'diffmaglim', 
                                'zpdiff',  'programid', 'obsdate', 'jd',  
                                'scifilename', 'diffilename', 'rfid',  
                                'refmaglim', 'refbckgnd', 'refsigpix', 
                                'zpref', 'refcreated', 'refjdstart',
                                'refjdend', 'reffilename', 'forcediffimflux', 
                                'forcediffimfluxunc', 'forcediffimsnr', 
                                'forcediffimchisq', 'forcediffimfluxap', 
                                'forcediffimfluxuncap', 'forcediffimsnrap',
                                'aperturecorr',  'dnearestrefsrc', 
                                'nearestrefmag', 'nearestrefmagunc',  
                                'nearestrefchi', 'nearestrefsharp', 
                                'procstatus'])
        
    elif ipac_version == '3':
        fp = pd.read_csv(fps_file, 
                         delim_whitespace=True, comment='#', skiprows=57,      
                         names=['index', 'field', 'ccdid', 'qid', 'filter',  
                                'pid', 'infobitssci', 'sciinpseeing', 
                                'scibckgnd', 'scisigpix',  'zpmaginpsci', 
                                'zpmaginpsciunc', 'zpmaginpscirms',  
                                'clrcoeff', 'clrcoeffunc', 'ncalmatches',  
                                'exptime', 'adpctdif1', 'adpctdif2', 
                                'diffmaglim',  'zpdiff', 'programid', 'jd', 
                                'rfid',  'forcediffimflux', 
                                'forcediffimfluxunc', 'forcediffimsnr', 
                                'forcediffimchisq',  'forcediffimfluxap', 
                                'forcediffimfluxuncap', 'forcediffimsnrap',  
                                'aperturecorr', 'dnearestrefsrc', 
                                'nearestrefmag', 'nearestrefmagunc',  
                                'nearestrefchi', 'nearestrefsharp',  
                                'refjdstart', 'refjdend', 'procstatus'])
        
        palomar = EarthLocation.of_site('Palomar')
        ha = pd.read_csv(fps_file, skiprows=3, nrows=2, 
                         delim_whitespace=True, 
                         names=['dum', 'dum1', 'dum2', 'coords', 'dum4'])
        ra, dec = ha.coords.values
        targ = SkyCoord(ra*u.deg, dec*u.deg)
        time = Time(fp.jd.values, format='jd')
        targaltaz = targ.transform_to(AltAz(obstime=time,location=palomar))
        airmass = targaltaz.secz.value
        fp['airmass'] = airmass
        
    
    union_list = ['field', 'ccdid', 'qid', 'filter', 
                  'pid', 'infobitssci', 'airmass',
                  'sciinpseeing', 'scibckgnd', 'scisigpix', 
                  'zpmaginpsci', 'zpmaginpsciunc', 'zpmaginpscirms', 
                  'clrcoeff', 'clrcoeffunc', 'ncalmatches', 
                  'exptime', 'diffmaglim', 'zpdiff', 'programid', 'jd', 
                  'rfid', 'forcediffimflux', 'forcediffimfluxunc', 
                  'forcediffimsnr', 'forcediffimchisq', 'forcediffimfluxap', 
                  'forcediffimfluxuncap', 'forcediffimsnrap', 'aperturecorr', 
                  'dnearestrefsrc', 'nearestrefmag', 'nearestrefmagunc', 
                  'nearestrefchi', 'nearestrefsharp', 
                  'refjdstart', 'refjdend', 'procstatus']
        
    fp_det = fp[union_list].dropna(subset=['forcediffimflux']).copy()

    
    cloudy = np.zeros_like(fp_det.infobitssci.values)
    
    rcid_rel_path = 'cal_data/zp_thresholds_quadID.txt'
    rcid_file = pkg_resources.resource_stream(__name__, rcid_rel_path)
    rcid_df = pd.read_csv(rcid_file, delim_whitespace=True, 
                          names=['zp_rcid_g', 'zp_rcid_r', 'zp_rcid_i'],
                          comment='#')

    rcid = (fp_det.ccdid.values-1)*4 + fp_det.qid.values
    zp_rcid_g = rcid_df.zp_rcid_g.iloc[rcid]
    zp_rcid_r = rcid_df.zp_rcid_r.iloc[rcid]
    zp_rcid_i = rcid_df.zp_rcid_i.iloc[rcid]
    flag_clouds = np.where( ( (fp_det['filter'].values == 'ZTF_g') &
                              ( (fp_det.zpmaginpsci.values > 26.7 - 
                                 0.2*fp_det.airmass.values
                                 ) | 
                                 (fp_det.zpmaginpscirms.values > 0.06) | 
                                 (fp_det.ncalmatches.values < 80) | 
                                 (fp_det.zpmaginpsci.values < zp_rcid_g - 
                                  0.2*fp_det.airmass.values
                                 )
                              )
                            ) | 
                            ( (fp_det['filter'].values == 'ZTF_r') &
                              ( (fp_det.zpmaginpsci.values > 26.65 - 
                                 0.15*fp_det.airmass.values
                                ) |  
                                (fp_det.zpmaginpscirms.values > 0.05) | 
                                (fp_det.ncalmatches.values < 120) | 
                                (fp_det.zpmaginpsci.values < zp_rcid_r - 
                                 0.15*fp_det.airmass.values
                                )
                              )
                            ) |
                            ( (fp_det['filter'].values == 'ZTF_i') &
                              ( (fp_det.zpmaginpsci.values > 26.0 - 
                                 0.07*fp_det.airmass.values
                                 ) |  
                                 (fp_det.zpmaginpscirms.values > 0.06) | 
                                 (fp_det.ncalmatches.values < 100) | 
                                 (fp_det.zpmaginpsci.values < zp_rcid_i - 
                                  0.07*fp_det.airmass.values
                                  )
                              )
                            )
                          )
    cloudy[flag_clouds] = int(2**25)
    fp_det['infobitssci'] = fp_det.infobitssci.values + cloudy    
    
    
    obs_jd = Time(fp_det.jd.values, format='jd')
    
    fp_det = fp_det.set_index(pd.to_datetime(obs_jd.datetime)) 
    
    fid_dict = {'ZTF_g': 1, 'ZTF_r': 2, 'ZTF_i': 3}
    fcqfid = np.array(fp_det.field.values*10000 + 
                      fp_det.ccdid.values*100 + 
                      fp_det.qid.values*10 + 
                      np.array([fid_dict[x] for x in fp_det['filter'].values])
                     )
    fp_det['fcqfid'] = fcqfid

    return fp_det

