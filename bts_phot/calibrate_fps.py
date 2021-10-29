import pandas as pd
from packaging import version

if version.parse(pd.__version__) < version.parse("1.3.0"):
    raise AssertionError("pandas version must be >= 1.3.0")

import glob, gc, os, pkg_resources
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.pyplot import Figure

from scipy import optimize as op
from scipy.stats import median_abs_deviation

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation, SkyCoord, AltAz


def read_ipac_fps(fps_file, rcid_path='cal_data/zp_thresholds_quadID.txt'):
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

    rcid_path: str (default = 'cal_data/zp_thresholds_quadID.txt')
    
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
    read_opts = {
        'delim_whitespace': True,
        'names': ['zp_rcid_g', 'zp_rcid_r', 'zp_rcid_i'],
        'comment': '#'
    }
    if os.path.exists(rcid_path):
        with open(rcid_path, "r") as f:
            rcid_df = pd.read_csv(f, **read_opts)
    else:
        rcid_df = pd.read_csv(
            pkg_resources.resource_stream(__name__, rcid_path),
            **read_opts
        )

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

def get_baseline(fps_file, window="10D", 
                 write_lc=False, 
                 make_plot=False, 
                 save_fig=False,
                 talk_plot=False,
                 save_path='default',
                 rcid_path='cal_data/zp_thresholds_quadID.txt'):
    
    """
    calculate the baseline region for an IPAC fps light curve and produce 
    flux calibrated photometry
    
    Parameters
    ----------
    fps_file : str or list of str
        File path to an IPAC table text file with fps results. If a list is 
        provided forced photometry is read from each file and the results 
        appended before producing the final output.
    
    window : int, offset, or BaseIndexer subclass (optional, default = '10D')
        Size of the moving window. This is the number of observations used 
        for calculating the rolling median. 
        If its an offset then this will be the time period of each window. 
        Each window will be a variable sized based on the observations 
        included in the time-period. This is only valid for datetimelike 
        indexes.
    
    write_lc : bool or DataFrame (optional, default = 'False')
        If True, the resulting calibrated photometry is written to a text file
        If a DataFrame is provided, it will be updated with the resulting 
        calibrated photometry
    
    make_plot : bool or matplotlib Figure (optional, default = 'False')
        If True, a plot of the calibrated light curve is created
        If a matplotlib figure is provided, it will be used for plotting 
        the calibrated lightcurve
    
    save_fig : bool (optional, default = 'False')
        If True, the resulting light curve plot is saved as a png
    
    talk_plot : bool (optional, default = 'False')
        If True, the saved plot has a transparent background and white 
        axes/axis labels for projection on a dark background
    
    save_path : str (optional, default = 'default')
        Path for writing light curve figures and ascii files 

    rcid_path: str (default = 'cal_data/zp_thresholds_quadID.txt')
    
    Returns
    -------
    fcqfid_dict : dict
        Dictionary with summary statistics for the baseline region for each 
        fcqfid
                              
    """

    if type(fps_file) == str:
        ztf_name = fps_file.split('forcedphotometry_')[1].split('_')[0]
        fp_df = read_ipac_fps(fps_file, rcid_path)
        if save_path == 'default':
            save_path = fps_file.split('forced')[0]
    elif type(fps_file) == list:
        ztf_names = [n.split('forcedphotometry_')[1].split('_')[0] for 
                     n in fps_file]
        if ztf_names.count(ztf_names[0]) == len(ztf_names):
            ztf_name = ztf_names[0]
            fp_df = read_ipac_fps(fps_file[0], rcid_path)
            for n in fps_file[1:]:
                extra_df = read_ipac_fps(n, rcid_path)
                fp_df = fp_df.append(extra_df)
            if save_path == 'default':
                save_path = fps_file[-1].split('forced')[0]
        else:
            raise AssertionError("""Photometric observations can only be 
                                    combined for a single SN""")
            
    unique_fid = np.unique(fp_df.fcqfid.values).astype(int)
    
    fcqfid_dict = {}
    t_peak_list = []
    
    for ufid in unique_fid:
        fcqfid_dict[str(ufid)] = {}
        bad_obs = np.zeros_like(fp_df.ccdid.values)
        bad_obs[np.where((fp_df.infobitssci.values > 0) | 
                         (fp_df.scisigpix.values > 25) | 
                         (fp_df.sciinpseeing.values > 5)
                        )] = 1        

        this_fcqfid = np.where((fp_df.fcqfid.values == ufid) & 
                               (bad_obs == 0))
        
        if ((ufid % 10 == 3) or (len(this_fcqfid[0]) < 2)):
            continue
        else:
            fcqf_df = fp_df.iloc[this_fcqfid].copy()
            flux_series = fcqf_df.forcediffimflux
            roll_med = flux_series.rolling(window, 
                                           center=True).median().values
            t_max = fcqf_df.jd.values[np.argmax(roll_med)]
            f_max_snr = fcqf_df.forcediffimsnr.values[np.argmax(roll_med)]
            flux_max = np.max(roll_med)
            flux_scatt = median_abs_deviation(fcqf_df.forcediffimflux.values, 
                                              scale='normal')
            peak_snr = flux_max/flux_scatt
            if ((peak_snr > 5) and (ufid < 10000000)):
                fcqfid_dict[str(ufid)]['det_sn'] = True
                fcqfid_dict[str(ufid)]['t_max'] = t_max
                t_peak_list.append(t_max)
            else:
                fcqfid_dict[str(ufid)]['det_sn'] = False
                
    
    if len(t_peak_list) > 0:
        t_peak = np.mean(t_peak_list)
        if len(t_peak_list) > 1 and np.std(t_peak_list, ddof=1) > 10:
            print('Warning! Large scatter in time of maximum')
        fcqfid_dict['t_peak'] = t_peak
        around_max = np.where((fp_df.jd.values - t_peak > - 10) & 
                              (fp_df.jd.values - t_peak < 10))
        if len(around_max[0]) > 0:
            diff_flux_around_max = fp_df.forcediffimflux.values[around_max]
            mag_min = np.nanmin(fp_df.zpdiff.values[around_max] -
                                2.5*np.log10(diff_flux_around_max))
            #calculate time when SN signal is "gone" via Co56 decay at z ~ 0.09
            t_faded = t_peak + (22.5 - mag_min)/0.009
        else:
            t_faded = t_peak + 611 # catch strange cases where t_gmax != t_rmax
        for ufid in unique_fid:
            if ufid % 10 == 4:
                continue
            else:                
                this_fcqfid = np.where(fp_df.fcqfid.values == ufid)
                fcqf_df = fp_df.iloc[this_fcqfid].copy()
                
                # measure the baseline pre-peak
                pre_bl = np.where((t_peak - fcqf_df.jd.values > 100) & 
                                        (fcqf_df.infobitssci.values == 0) & 
                                        (fcqf_df.scisigpix.values < 25) & 
                                        (fcqf_df.sciinpseeing.values < 5)
                                       )
                if len(pre_bl[0]) > 1:
                    base_jd = fcqf_df.jd.values[pre_bl]
                    base_flux = fcqf_df.forcediffimflux.values[pre_bl]
                    base_flux_unc = fcqf_df.forcediffimfluxunc.values[pre_bl]
                    mask = np.where(np.abs( (base_flux - np.median(base_flux)
                                            )/base_flux_unc
                                          ) <= 5)
                    if len(mask[0]) > 1:
                        Cmean = np.average(base_flux[mask], 
                                           weights=1/base_flux_unc[mask]**2)
                        sum_diff_sq = np.sum( ( (base_flux[mask] - Cmean) / 
                                                (base_flux_unc[mask])
                                              )**2
                                            )
                        chi = 1/(len(mask[0])-1)*sum_diff_sq
                        fcqfid_dict[str(ufid)]['C_pre'] = Cmean
                        fcqfid_dict[str(ufid)]['chi_pre'] = chi
                    fcqfid_dict[str(ufid)]['N_pre_peak'] = len(mask[0])
                else:
                    fcqfid_dict[str(ufid)]['N_pre_peak'] = 0

                # measure the baseline post-peak
                post_bl = np.where((fcqf_df.jd.values > t_faded) & 
                                   (fcqf_df.infobitssci.values == 0) & 
                                   (fcqf_df.scisigpix.values < 25) & 
                                   (fcqf_df.sciinpseeing.values < 5)
                                  )

                if len(post_bl[0]) > 1:
                    base_jd = fcqf_df.jd.values[post_bl]
                    base_flux = fcqf_df.forcediffimflux.values[post_bl]
                    base_flux_unc = fcqf_df.forcediffimfluxunc.values[post_bl]

                    mask = np.where(np.abs( (base_flux - np.median(base_flux)
                                            )/base_flux_unc
                                          ) <= 5)
                    if len(mask[0]) > 1:
                        Cmean = np.average(base_flux[mask], 
                                           weights=1/base_flux_unc[mask]**2)
                        sum_diff_sq = np.sum( ( (base_flux[mask] - Cmean) / 
                                                (base_flux_unc[mask])
                                              )**2
                                            )
                        chi = 1/(len(mask[0])-1)*sum_diff_sq
                        fcqfid_dict[str(ufid)]['C_post'] = Cmean
                        fcqfid_dict[str(ufid)]['chi_post'] = chi
                    fcqfid_dict[str(ufid)]['N_post_peak'] = len(mask[0])
                else:
                    fcqfid_dict[str(ufid)]['N_post_peak'] = 0

    if make_plot is not False or write_lc is not False:
        fnu_microJy = -999.*np.ones_like(fp_df.forcediffimflux.values)
        fnu_microJy_unc = -999.*np.ones_like(fp_df.forcediffimflux.values)
        n_base_obs = np.zeros_like(fp_df.forcediffimflux.values).astype(int)
        which_base = np.zeros_like(fp_df.forcediffimflux.values).astype(int)
        C_baseline = np.zeros_like(fp_df.forcediffimflux.values).astype(int)
        sys_sigma = np.zeros_like(fp_df.forcediffimflux.values).astype(int)
        
        bad_obs = np.zeros_like(fp_df.ccdid.values)
        bad_obs[np.where((fp_df.infobitssci.values > 0) | 
                         (fp_df.scisigpix.values > 25) | 
                         (fp_df.sciinpseeing.values > 5)
                        )] = 1        
        
        for key in fcqfid_dict:
            if (key != 't_peak' and 
                'N_pre_peak' in fcqfid_dict[key].keys()):
                ufid = int(key)
                this_fcqfid = np.where(fp_df.fcqfid.values == ufid)

                if ( (fcqfid_dict[key]['N_pre_peak'] >= 25) or 
                     ( (fcqfid_dict[key]['N_pre_peak'] >= 10) and 
                       (fcqfid_dict[key]['N_post_peak'] < 25)
                     )
                   ):
                    baseline = fcqfid_dict[key]['C_pre']
                    multiplier = max(np.sqrt(fcqfid_dict[key]['chi_pre']), 1)
                    n_baseline = fcqfid_dict[key]['N_pre_peak']
                    pre_or_post = -1
                elif ( (fcqfid_dict[key]['N_post_peak'] >= 25) or 
                       ( (fcqfid_dict[key]['N_pre_peak'] < 10) and 
                         (fcqfid_dict[key]['N_post_peak'] >= 10)
                       )
                     ):
                    baseline = fcqfid_dict[key]['C_post']
                    multiplier = max(np.sqrt(fcqfid_dict[key]['chi_post']), 1)
                    n_baseline = fcqfid_dict[key]['N_post_peak']
                    pre_or_post = 1
                else:
                    n_base_obs[this_fcqfid] = fcqfid_dict[key]['N_pre_peak']
                    which_base[this_fcqfid] = -1
                    continue

                flux_dn = fp_df.forcediffimflux.values[this_fcqfid] - baseline 
                unc_fcqfid = fp_df.forcediffimfluxunc.values[this_fcqfid]
                flux_dn_unc = unc_fcqfid * multiplier
                zp_fcqfid = fp_df.zpdiff.values[this_fcqfid]
                fnu_microJy[this_fcqfid] = flux_dn*10**(29 - 
                                                        48.6/2.5 - 
                                                        0.4*zp_fcqfid)
                fnu_microJy_unc[this_fcqfid] = flux_dn_unc*10**(29 - 
                                                                48.6/2.5 - 
                                                                0.4*zp_fcqfid)
                n_base_obs[this_fcqfid] = n_baseline
                which_base[this_fcqfid] = pre_or_post
                C_baseline[this_fcqfid] = baseline
                sys_sigma[this_fcqfid] = multiplier

    
    if write_lc is not False:

        if isinstance(write_lc, pd.DataFrame):
            write_df = write_lc
            write_df['jd'] = fp_df.jd.values
        else:
            write_df = pd.DataFrame(fp_df.jd.values, columns=['jd'])
        write_df['fnu_microJy'] = fnu_microJy
        write_df['fnu_microJy_unc'] = fnu_microJy_unc
        write_df['passband'] = fp_df['filter'].values
        write_df['programid'] = fp_df.programid.values
        write_df['fcqfid'] = fp_df.fcqfid.values
        write_df['N_baseline'] = n_base_obs
        write_df['pre_or_post'] = which_base
        write_df['zpdiff'] = fp_df.zpdiff.values
        write_df['C'] = C_baseline
        write_df['sys_unc_factor'] = sys_sigma
        write_df['poor_conditions'] = bad_obs
        gr_obs = np.where(write_df.fcqfid.values % 10 != 4)
        if not isinstance(write_lc, pd.DataFrame):
            fname = save_path + ztf_name + '_fnu.csv'
            write_df.iloc[gr_obs].to_csv(fname, index=False)                
                
    if make_plot is not False:

        color_dict = {1: 'MediumAquaMarine', 2: 'Crimson', 3: 'Goldenrod'}
        nplots = 0
        jdstart = 2458119.5
        for key in fcqfid_dict:
            if key != 't_peak' and 'N_pre_peak' in fcqfid_dict[key].keys():
                if ((fcqfid_dict[key]['N_pre_peak'] > 9) or 
                    (fcqfid_dict[key]['N_post_peak'] > 9)
                   ):
                    nplots += 1

        if nplots > 0:
            if make_plot is True:
                fig = plt.figure(figsize=(8,nplots*3 + 0.5))
            else:
                fig = make_plot 
            axes = fig.subplots(nplots, 1, sharex=True)
            plot_num = 0
            for key in fcqfid_dict:
                if ((key != 't_peak') and 
                    ('N_pre_peak' in fcqfid_dict[key].keys()) and
                    (fcqfid_dict[key]['N_pre_peak'] >= 10 or 
                     fcqfid_dict[key]['N_post_peak'] >= 10)
                   ):
                    ufid = int(key)
                    this_fcqfid_good = np.where((fp_df.fcqfid.values == ufid) &
                                                (bad_obs == 0))

                    plot_jd = fp_df.jd.values[this_fcqfid_good] - jdstart
                    plot_flux = fnu_microJy[this_fcqfid_good]
                    plot_flux_unc = fnu_microJy_unc[this_fcqfid_good]

                    ax = axes[plot_num] if nplots > 1 else axes
                    ax.errorbar(plot_jd, plot_flux, plot_flux_unc,
                                            fmt='o', 
                                            mec=color_dict[ufid % 10], 
                                            ecolor=color_dict[ufid % 10],
                                            mfc='None')
                    ax.axvline(x = t_peak - jdstart, color = '0.5', ls = '--')
                    ax.axhline(y = 0, color = '0.5',
                               ls = (0, [8,1.5,1,1.5]), lw = 0.5, alpha=0.75)
                    ax.axvspan(0, t_peak - jdstart - 100, 
                               color='Cornsilk', alpha=0.6, lw=0)
                    ax.axvspan(t_faded - jdstart, 1e6, 
                               color='Cornsilk', alpha=0.6, lw=0)
                    ax.set_ylabel(r'flux ($\mu$Jy)', fontsize = 14)
                    ax.set_xlim(np.min(fp_df.jd.values - jdstart)-10, 
                                np.max(fp_df.jd.values - jdstart)+10)
                    if talk_plot:
                        ax.tick_params(axis='both', colors='white')
                        for spine in ['top','bottom','left','right']:
                            ax.spines[spine].set_color('white')
                        ax.yaxis.label.set_color('white')
                        ax.xaxis.label.set_color('white')
                    plot_num += 1
                        
            ax.set_xlabel('Time (JD - 2018 Jan 01)', fontsize = 14)

            fig.tight_layout()
            if save_fig:
                pname = save_path + ztf_name + '_fnu.png'                
                fig.savefig(pname)
                if not isinstance(make_plot, Figure):
                    plt.close(fig)
                    plt.close('all')
                    del(fig)
                    gc.collect()

    return fcqfid_dict
