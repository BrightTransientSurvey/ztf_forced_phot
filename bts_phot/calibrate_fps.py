import gc, pkg_resources
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.pyplot import Figure

from scipy.stats import median_abs_deviation

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation, SkyCoord, AltAz

from supersmoother import SuperSmoother

pkg_resources.require("pandas>=1.3")


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
                                'scimaglim', 'zpmaginpsci', 'zpmaginpsciunc',
                                'zpmaginpscirms', 'clrcoeff', 'clrcoeffunc',
                                'ncalmatches', 'exptime', 'diffmaglim',
                                'zpdiff', 'programid', 'obsdate', 'jd',
                                'scifilename', 'diffilename', 'rfid',
                                'refmaglim', 'refbckgnd', 'refsigpix',
                                'zpref', 'refcreated', 'refjdstart',
                                'refjdend', 'reffilename', 'forcediffimflux',
                                'forcediffimfluxunc', 'forcediffimsnr',
                                'forcediffimchisq', 'forcediffimfluxap',
                                'forcediffimfluxuncap', 'forcediffimsnrap',
                                'aperturecorr', 'dnearestrefsrc',
                                'nearestrefmag', 'nearestrefmagunc',
                                'nearestrefchi', 'nearestrefsharp',
                                'procstatus'])

    elif int(ipac_version) >= 3:
        fp = pd.read_csv(fps_file,
                         delim_whitespace=True, comment='#', header=0,
                         names=['sindex', 'field', 'ccdid', 'qid', 'filter',
                                'pid', 'infobitssci', 'sciinpseeing',
                                'scibckgnd', 'scisigpix', 'zpmaginpsci',
                                'zpmaginpsciunc', 'zpmaginpscirms',
                                'clrcoeff', 'clrcoeffunc', 'ncalmatches',
                                'exptime', 'adpctdif1', 'adpctdif2',
                                'diffmaglim', 'zpdiff', 'programid', 'jd',
                                'rfid', 'forcediffimflux',
                                'forcediffimfluxunc', 'forcediffimsnr',
                                'forcediffimchisq', 'forcediffimfluxap',
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
        targaltaz = targ.transform_to(AltAz(obstime=time, location=palomar))
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

    rcid_df = pd.read_csv(
        pkg_resources.resource_stream(__name__, 
        'cal_data/zp_thresholds_quadID.txt'),
        **read_opts
    )

    rcid = (fp_det.ccdid.values-1)*4 + fp_det.qid.values-1
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

def get_baseline(fps_file, window="14D",
                 write_lc=False,
                 make_plot=False,
                 roll_med_plot = False,
                 save_fig=False,
                 talk_plot=False,
                 save_path='default', 
                 deprecated=False):

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
        
    roll_med_plot: bool (optional, default = 'False')
        If True, plot the rolling median of the calibrated photometry on top of 
        the calibrated light curve.

    save_fig : bool (optional, default = 'False')
        If True, the resulting light curve plot is saved as a png

    talk_plot : bool (optional, default = 'False')
        If True, the saved plot has a transparent background and white
        axes/axis labels for projection on a dark background

    save_path : str (optional, default = 'default')
        Path for writing light curve figures and ascii files
    
    deprecated : bool (optional, default = 'False')
        Return "old" outputs from the original version of the software circa 
        01 Jan 2022. This will be removed eventually.

    Returns
    -------
    fcqfid_dict : dict
        Dictionary with summary statistics for the baseline region for each
        fcqfid
    """

# # Summary of flags >> to be removed later once flags are locked
# flag | bit value | meaning
# 1    | 1    | large scatter in max
# 2    | 2    | possible emission in pre-SN baseline
# 3    | 4    | possible emission in post-SN baseline
# 4    | 8    | outliers in the baseline region
# 5    | 16   | large scatter in baseline region
# 6    | 32   | N_baseline < 10
# 7    | 64   | scisigpix > 25
# 8    | 128  | sciinpseeing > 5
# 9    | 256  | uncertainties underestimated
# 10   | 512  | infobits > 0
# 11   | 1024 | no flux detected in forced phot
    
    bad_obs_fl = 512
    
    if isinstance(fps_file, str):
        if save_fig is True or write_lc is not False or make_plot is True:
            ztf_name = 'ZTF' + fps_file.split('ZTF')[-1][0:9]
            print(ztf_name)
        fp_df = read_ipac_fps(fps_file)
        if save_path == 'default':
            save_path = fps_file[0:-len(fps_file.split('/')[-1])]
    elif isinstance(fps_file, list):
        ztf_names = [n.split('forcedphotometry_')[1].split('_')[0] for
                     n in fps_file]
        if ztf_names.count(ztf_names[0]) == len(ztf_names):
            ztf_name = ztf_names[0]
            fp_df = read_ipac_fps(fps_file[0])
            for n in fps_file[1:]:
                extra_df = read_ipac_fps(n)
                fp_df = fp_df.append(extra_df)
            if save_path == 'default':
                save_path = fps_file[-1].split('forced')[0]
        else:
            raise AssertionError("""Photometric observations can only be
                                    combined for a single SN""")
    fp_df.drop_duplicates(['pid', 'forcediffimflux', 'forcediffimfluxunc'],
                          inplace=True)
    
    unique_fid = np.unique(fp_df.fcqfid.values).astype(int)

    fcqfid_dict = {}
    t_peak_list = []

    fp_df['flags'] = np.zeros(len(fp_df)).astype(int)
    fp_df.loc[fp_df.scisigpix.values > 25, 'flags'] += 64
    fp_df.loc[fp_df.sciinpseeing.values > 5, 'flags'] += 128
    fp_df.loc[fp_df.infobitssci.values > 0, 'flags'] += 512
    fp_df.loc[fp_df.forcediffimfluxunc == -99999, 'flags'] += 1024
    bad_obs = np.where(fp_df['flags'].values >= bad_obs_fl, 1, 0)

    for ufid in unique_fid:
        fcqfid_dict[str(ufid)] = {}
        fcqfid_dict[str(ufid)]['N_bl'] = 0
        # bad_obs = np.zeros_like(fp_df.ccdid.values)
        # if deprecated:
        #     bad_obs[np.where((fp_df.infobitssci.values > 0) |
        #                      (fp_df.scisigpix.values > 25) |
        #                      (fp_df.sciinpseeing.values > 5)
        #                     )] = 1
        # else:
        #     # AW: added additional filter
        #     bad_obs[np.where((fp_df.infobitssci.values > 0) |
        #                      (fp_df.forcediffimfluxunc == -99999)
        #                     )] = 1

        this_fcqfid = np.where((fp_df.fcqfid.values == ufid) &
                               (fp_df['flags'].values < bad_obs_fl))

        if ((ufid % 10 == 3) or (len(this_fcqfid[0]) < 2)):
            continue
        else:
            fcqf_df = fp_df.iloc[this_fcqfid].copy()
            flux_series = fcqf_df.forcediffimflux
            
            roll_med = flux_series.rolling(window,
                                           center=True).median().values
            roll_peak = np.argmax(roll_med)
            t_max = fcqf_df.jd.values[roll_peak]
            flux_max = np.max(roll_med)
            flux_scatt = median_abs_deviation(fcqf_df.forcediffimflux.values,
                                              scale='normal')
            max_over_scatt = flux_max/flux_scatt
            peak_snr = flux_max/fcqf_df.forcediffimfluxunc.values[roll_peak]
            if (max_over_scatt > 5 and peak_snr > 5):
                fcqfid_dict[str(ufid)]['det_sn'] = True
                fcqfid_dict[str(ufid)]['t_fcqfid_max'] = t_max
                if ufid < 10000000:
                    t_peak_list.append(t_max)
            else:
                fcqfid_dict[str(ufid)]['det_sn'] = False


    if len(t_peak_list) > 0:
        t_peak = np.mean(t_peak_list)
        for ufid in unique_fid:
            fcqfid_dict[str(ufid)]['t_peak'] = t_peak
            if len(t_peak_list) > 1:
                fcqfid_dict[str(ufid)]['peak_scatter'] = np.std(t_peak_list, ddof=1)
                if np.std(t_peak_list, ddof=1) > 50:
                    print('Warning! Large scatter in time of maximum')
                    fcqfid_dict[str(ufid)]['peak_warning'] = 'large scat in t_peak'
                    fp_df['flags'] += 1 
        
        around_max = np.where((fp_df.jd.values - t_peak > - 10) &
                              (fp_df.jd.values - t_peak < 10) & 
                              (fp_df.forcediffimflux.values > 0))
        if len(around_max[0]) > 0:
            diff_flux_around_max = fp_df.forcediffimflux.values[around_max]
            mag_min = np.nanmin(fp_df.zpdiff.values[around_max] -
                                2.5*np.log10(diff_flux_around_max))
            #calculate time when SN signal is "gone" via Co56 decay at z ~ 0.09
            t_faded = t_peak + (22.5 - mag_min)/0.009
        else:
            t_faded = t_peak + 611 # catch strange cases where t_gmax != t_rmax
        for ufid in unique_fid:
            if ufid % 10 == 4: # check if this ever happens - should not be a filter id = 4
                continue
            else:
                this_fcqfid = np.where(fp_df.fcqfid.values == ufid)
                fcqf_df = fp_df.iloc[this_fcqfid].copy()

                # measure the baseline
                bl = np.where(((t_peak - fcqf_df.jd.values > 100) |
                              (fcqf_df.jd.values > t_faded)) &
                              (fcqf_df['flags'].values < bad_obs_fl)
                             )
                
                if len(bl[0]) > 1:
                    base_flux = fcqf_df.forcediffimflux.values[bl]
                    base_flux_unc = fcqf_df.forcediffimfluxunc.values[bl]
                    mask = np.where(np.abs( (base_flux - np.median(base_flux)
                                            )/base_flux_unc) <= 5)
                    outl = np.where(np.abs( (base_flux - np.median(base_flux)
                                            )/base_flux_unc) > 5)
                    # little kludgy but avoids warning
                    outl_bool = np.zeros(len(fp_df)).astype(bool)          
                    outl_bool[this_fcqfid[0][bl[0]][outl]] = True
                    if sum(outl[0]) > 0:           
                        fp_df.loc[outl_bool, 'flags'] += 8
                    if len(mask[0]) > 2:
                        non_det = base_flux[mask]
                        Cmean = np.average(non_det,
                                           weights=1/base_flux_unc[mask]**2)
                        sum_diff_sq = np.sum( ( (non_det - Cmean) /
                                                (base_flux_unc[mask])
                                              )**2
                                            )
                        chi = 1/(len(mask[0])-1)*sum_diff_sq
                        mean = np.mean(non_det)
                        median = np.median(non_det)
                        # using bootstrapping to place unc on median
                        medians = np.empty(1000)
                        for i in range(1000): 
                            medians[i] = np.median(np.random.choice(non_det, 
                                         size = len(non_det), replace = True))
                        median_unc = np.diff(np.percentile(medians, (16,84)))[0]/2
                        
                        low_limit = np.percentile(non_det, 10)
                        high_limit = np.percentile(non_det, 90)
                        trim = np.where((non_det > low_limit) &
                                        (non_det < high_limit)
                                        )
                        trim_mean = np.mean(non_det[trim])
                        # is this warning needed?
                        if len(non_det[trim]) == 0:
                            print(f'{ztf_name} {ufid} bl empty')
                        scatter = np.diff(np.percentile(non_det, (16,84)))[0]/2
                        
                        fcqfid_dict[str(ufid)]['C_bl'] = Cmean
                        fcqfid_dict[str(ufid)]['chi_bl'] = chi
                        fcqfid_dict[str(ufid)]['mean_bl'] = mean
                        fcqfid_dict[str(ufid)]['median_bl'] = median
                        fcqfid_dict[str(ufid)]['median_unc_bl'] = median_unc
                        fcqfid_dict[str(ufid)]['trim_mean_bl'] = trim_mean
                        fcqfid_dict[str(ufid)]['scatter_bl'] = scatter
                    fcqfid_dict[str(ufid)]['N_bl'] = len(mask[0])
                    if fcqfid_dict[str(ufid)]['N_bl'] < 10:
                        fp_df.loc[fp_df.fcqfid == ufid, 'flags'] += 32
                else:
                    fcqfid_dict[str(ufid)]['N_bl'] = 0
                    fp_df.loc[fp_df.fcqfid == ufid, 'flags'] += 32
                
                # measure the baseline pre-peak
                if deprecated:
                    pre_bl = np.where((t_peak - fcqf_df.jd.values > 100) &
                                            (fcqf_df.infobitssci.values == 0) &
                                            (fcqf_df.scisigpix.values < 25) &
                                            (fcqf_df.sciinpseeing.values < 5)
                                           )
                else:
                    pre_bl = np.where((t_peak - fcqf_df.jd.values > 100) &
                                      (fcqf_df['flags'].values < bad_obs_fl)
                                      )
                if len(pre_bl[0]) > 1:
                    base_flux = fcqf_df.forcediffimflux.values[pre_bl]
                    base_flux_unc = fcqf_df.forcediffimfluxunc.values[pre_bl]
                    mask = np.where(np.abs( (base_flux - np.median(base_flux)
                                            )/base_flux_unc
                                          ) <= 5)
                    if len(mask[0]) > 2:
                        non_det = base_flux[mask]
                        Cmean = np.average(non_det,
                                           weights=1/base_flux_unc[mask]**2)
                        sum_diff_sq = np.sum( ( (non_det - Cmean) /
                                                (base_flux_unc[mask])
                                              )**2
                                            )
                        chi = 1/(len(mask[0])-1)*sum_diff_sq
                        mean = np.mean(non_det)
                        median = np.median(non_det)
                        
                        medians = np.empty(1000)
                        for i in range(1000): 
                            medians[i] = np.median(np.random.choice(non_det, 
                                         size = len(non_det), replace = True))
                        median_unc = np.diff(np.percentile(medians, (16,84)))[0]/2
                        
                        low_limit = np.percentile(non_det, 10)
                        high_limit = np.percentile(non_det, 90)
                        trim = np.where((non_det > low_limit) &
                                        (non_det < high_limit)
                                        )
                        trim_mean = np.mean(non_det[trim])
                        if len(non_det[trim]) == 0:
                            print(f'{ztf_name} {ufid} pre empty')
                        scatter = np.diff(np.percentile(non_det, (16,84)))[0]/2
                        
                        fcqfid_dict[str(ufid)]['C_pre'] = Cmean
                        fcqfid_dict[str(ufid)]['chi_pre'] = chi
                        fcqfid_dict[str(ufid)]['mean_pre'] = mean
                        fcqfid_dict[str(ufid)]['median_pre'] = median
                        fcqfid_dict[str(ufid)]['median_unc_pre'] = median_unc
                        fcqfid_dict[str(ufid)]['trim_mean_pre'] = trim_mean
                        fcqfid_dict[str(ufid)]['scatter_pre'] = scatter
                    fcqfid_dict[str(ufid)]['N_pre_peak'] = len(mask[0])
                else:
                    fcqfid_dict[str(ufid)]['N_pre_peak'] = 0

                # measure the baseline post-peak
                if deprecated:
                    post_bl = np.where((fcqf_df.jd.values > t_faded) &
                                       (fcqf_df.infobitssci.values == 0) &
                                       (fcqf_df.scisigpix.values < 25) &
                                       (fcqf_df.sciinpseeing.values < 5)
                                      )
                else:
                    post_bl = np.where((fcqf_df.jd.values > t_faded) &
                                       (fcqf_df['flags'].values < bad_obs_fl)
                                      )
                    
                if len(post_bl[0]) > 1:
                    base_flux = fcqf_df.forcediffimflux.values[post_bl]
                    base_flux_unc = fcqf_df.forcediffimfluxunc.values[post_bl]

                    mask = np.where(np.abs( (base_flux - np.median(base_flux)
                                            )/base_flux_unc
                                          ) <= 5)
                    if len(mask[0]) > 2:
                        non_det = base_flux[mask]
                        Cmean = np.average(non_det,
                                           weights=1/base_flux_unc[mask]**2)
                        sum_diff_sq = np.sum( ( (base_flux[mask] - Cmean) /
                                                (base_flux_unc[mask])
                                              )**2
                                            )
                        chi = 1/(len(mask[0])-1)*sum_diff_sq
                        mean = np.mean(non_det)
                        median = np.median(non_det)
                        medians = np.empty(1000)
                        for i in range(1000): 
                            medians[i] = np.median(np.random.choice(non_det, 
                                         size = len(non_det), replace = True))
                        median_unc = 0.5*np.diff(np.percentile(medians, 
                                                           (16,84)))[0]
                        low_limit = np.percentile(non_det, 10)
                        high_limit = np.percentile(non_det, 90)
                        trim = np.where((non_det > low_limit) &
                                        (non_det < high_limit)
                                        )
                        trim_mean = np.mean(non_det[trim])
                        if len(non_det[trim]) == 0:
                            print(f'{ztf_name} {ufid} post empty')
                        scatter = np.diff(np.percentile(non_det, (16,84)))[0]/2
                        
                        fcqfid_dict[str(ufid)]['C_post'] = Cmean
                        fcqfid_dict[str(ufid)]['chi_post'] = chi
                        fcqfid_dict[str(ufid)]['mean_post'] = mean
                        fcqfid_dict[str(ufid)]['median_post'] = median
                        fcqfid_dict[str(ufid)]['median_unc_post'] = median_unc
                        fcqfid_dict[str(ufid)]['trim_mean_post'] = trim_mean
                        fcqfid_dict[str(ufid)]['scatter_post'] = scatter
                    fcqfid_dict[str(ufid)]['N_post_peak'] = len(mask[0])
                else:
                    fcqfid_dict[str(ufid)]['N_post_peak'] = 0

    if make_plot is not False or write_lc is not False:
        fnu_microJy = -999.*np.ones_like(fp_df.forcediffimflux.values)
        fnu_microJy_unc = -999.*np.ones_like(fp_df.forcediffimflux.values)
        
        sys_sigma = np.zeros_like(fp_df.forcediffimflux.values)
        if deprecated:
            n_base_obs = np.zeros_like(fp_df.forcediffimflux.values).astype(int)
            which_base = np.zeros_like(fp_df.forcediffimflux.values).astype(int)
            C_baseline = np.zeros_like(fp_df.forcediffimflux.values)

        for key in fcqfid_dict:
            if fcqfid_dict[key]['N_bl'] > 1:
                ufid = int(key)
                this_fcqfid = np.where(fp_df.fcqfid.values == ufid)
                if deprecated:
                    sys_unc = max(fcqfid_dict[key]['chi_pre']**0.5, 1)
                else:
                    good_fcqfid = np.where((fp_df.fcqfid.values == ufid) & 
                                           (fp_df['flags'].values < bad_obs_fl))
                    chi_ser = fp_df.forcediffimchisq.iloc[good_fcqfid].copy()
                    good_diffl = fp_df.forcediffimflux.iloc[good_fcqfid].values
                    this_diffl = fp_df.forcediffimflux.iloc[this_fcqfid].values
                    med_chi = np.median(chi_ser.values)
                    if med_chi < 1.5:
                        sys_unc = med_chi**0.5 * np.ones_like(this_diffl)
                    elif len(good_fcqfid[0]) > 0: 
                        try:
                            model = SuperSmoother()
                            model.fit(fp_df.forcediffimflux.iloc[good_fcqfid], 
                                      fp_df.forcediffimchisq.iloc[good_fcqfid])
                            # find the smoothed fit to the data
                            # interpolate to non-good obs
                            yfit = model.predict(np.sort(good_diffl)) 
                            chi_interp = np.interp(this_diffl,
                                                np.sort(good_diffl), 
                                                yfit)
                            chi_interp = np.where(chi_interp < 1, 1, 
                                                  chi_interp)
                            sys_unc = chi_interp**0.5
                            # is this warning necessary?
                            if np.isnan(sys_unc).any():
                                print(f'{ztf_name} {ufid} bad interp')
                        except:                            
                            chi_ser.index = pd.to_datetime(good_diffl, 
                                                           unit='s', 
                                                           origin='unix')
                            chi_ser.sort_index(inplace=True)
                            runmed = chi_ser.rolling('1000s', 
                                                     center=True).median()
                            sys_unc = np.interp(this_diffl, 
                                                np.sort(good_diffl), 
                                                runmed.values)**0.5
                    else:
                        continue
                sys_unc = np.where(sys_unc < 1, 1, sys_unc)

                if deprecated:
                    if ( (fcqfid_dict[key]['N_pre_peak'] >= 25) or
                         ( (fcqfid_dict[key]['N_pre_peak'] >= 10) and
                           (fcqfid_dict[key]['N_post_peak'] < 25)
                         )
                       ):
                        baseline = fcqfid_dict[key]['C_pre']                        
                        n_baseline = fcqfid_dict[key]['N_pre_peak']
                        pre_or_post = -1
                    elif ( (fcqfid_dict[key]['N_post_peak'] >= 25) or
                           ( (fcqfid_dict[key]['N_pre_peak'] < 10) and
                             (fcqfid_dict[key]['N_post_peak'] >= 10)
                           )
                         ):
                        baseline = fcqfid_dict[key]['C_post']
                        n_baseline = fcqfid_dict[key]['N_post_peak']
                        pre_or_post = 1
                    else:
                        n_base_obs[this_fcqfid] = fcqfid_dict[key]['N_pre_peak']
                        which_base[this_fcqfid] = -1
                        continue
                else:
                    # determine if there is emission in pre/post SN baseline
                    if 'C_bl' not in fcqfid_dict[key].keys():
                        print('{ztf_name}, {key} no C_bl'.format(ztf_name = 
                                                                 ztf_name, 
                                                                 key = key))
                        continue
                    baseline = fcqfid_dict[key]['median_bl']
                    baseline_unc = fcqfid_dict[key]['median_unc_bl']
                    base_scat = fcqfid_dict[key]['scatter_bl']
                    fcqfid_dict[key]['which_bl'] = 'pre+post SN'
                    
                    good_df = fp_df.iloc[good_fcqfid].copy()
                    flux_series = good_df.forcediffimflux
                    roll_med = flux_series.rolling("14D",
                                                   center=True).median().values
                    this_fl = fp_df['flags'].iloc[this_fcqfid].values
                    scale = sys_unc[np.where(this_fl < bad_obs_fl)]
                    scale_unc = scale * good_df.forcediffimfluxunc.values

                    pre_bl = np.where(t_peak - good_df.jd.values > 100)
                    pre_em = (roll_med[pre_bl] - baseline)/scale_unc[pre_bl]
                    
                    post_bl = np.where(good_df.jd.values > t_faded)
                    post_em = (roll_med[post_bl] - baseline)/scale_unc[post_bl]

                    # test scatter after systematic correction
                    bl = np.where(((t_peak - good_df.jd.values > 100) |
                                  (good_df.jd.values > t_faded)) & 
                                  (good_df['flags'].values & 8 == 0)
                                 )
                    chi2nu = np.sum(((good_df.iloc[bl].forcediffimflux.values -
                                      baseline) /
                                     scale_unc[bl])**2)/len(good_df.iloc[bl])
                    fcqfid_dict[key]['chi2nu'] = chi2nu
                    if chi2nu > 2:
                        print('Warning! scaled unc are underestimated')
                        print(f'{ztf_name} {key} has chi2nu = {chi2nu:.3f}')
                        fp_df.loc[fp_df.fcqfid == ufid, 'flags'] += 256

                    #### So many nested ifs - cover your eye's Guido
                    if len(pre_em) < 10 and len(post_em) < 10:
                        continue
                    
                    pre_rise_em = False
                    if len(pre_em) == 0:
                        fcqfid_dict[key]['which_bl'] = 'post SN'
                    if len(pre_em) >= 1:
                        if len(pre_bl[0]) < 5 or sum(pre_em >= 7) >= 2:
                            if ((len(post_em) >= 10) and 
                                (fcqfid_dict[key]['N_post_peak'] > 2)):
                                baseline = fcqfid_dict[key]['median_post']
                                baseline_unc = fcqfid_dict[key]['median_unc_post']
                                base_scat = fcqfid_dict[key]['scatter_post']
                                fcqfid_dict[key]['which_bl'] = 'post SN'
                            if sum(pre_em >= 7) >= 2:
                                print(f'Warning {ztf_name} {ufid} pre-SN')
                                fcqfid_dict[key]['Warning'] = 'pre-SN emission'
                                pre_rise_em = True
                                fp_df.loc[fp_df.fcqfid == ufid, 'flags'] += 2
                                
                    post_rise_em = False
                    if len(post_em) == 0:
                        fcqfid_dict[key]['which_bl'] = 'pre SN'
                    if len(post_em) >= 1:
                        if len(post_em) < 5 or sum(post_em >= 7) >= 2:
                            if ((len(pre_em) >= 10) and
                                (fcqfid_dict[key]['N_pre_peak'] > 2)):
                                baseline = fcqfid_dict[key]['median_pre']
                                baseline_unc = fcqfid_dict[key]['median_unc_pre']
                                base_scat = fcqfid_dict[key]['scatter_pre']
                                fcqfid_dict[key]['which_bl'] = 'pre SN'
                            if sum(post_em >= 7) >= 2:
                                print(f'Warning {ztf_name} {ufid} post-SN')
                                fcqfid_dict[key]['Warning'] = 'post-SN emission'
                                post_rise_em = True
                                fp_df.loc[fp_df.fcqfid == ufid, 'flags'] += 4
                                
                    if pre_rise_em + post_rise_em == 2:
                        baseline = fcqfid_dict[key]['median_bl']
                        baseline_unc = fcqfid_dict[key]['median_unc_bl']
                        base_scat = fcqfid_dict[key]['scatter_bl']
                        print(f'Warning {ztf_name} {ufid} bad baseline')
                        fcqfid_dict[key]['Warning'] = 'bad baseline'
                        fcqfid_dict[key]['which_bl'] = 'pre+post SN'
                        # fp_df.loc[fp_df.fcqfid == ufid, 'flags'] += 512
                if base_scat > 100:
                    fp_df.loc[fp_df.fcqfid == ufid, 'flags'] += 16
                flux_dn = fp_df.forcediffimflux.values[this_fcqfid] - baseline
                unc_fcqfid = fp_df.forcediffimfluxunc.values[this_fcqfid]

                # AW: fixed bug, add deprecated version
                if deprecated:
                    flux_dn_unc = np.sqrt(unc_fcqfid**2) * sys_unc
                else:
                    flux_dn_unc = np.sqrt(unc_fcqfid**2 + baseline_unc**2) * sys_unc

                zp_fcqfid = fp_df.zpdiff.values[this_fcqfid]
                fnu_microJy[this_fcqfid] = flux_dn*10**(29 -
                                                        48.6/2.5 -
                                                        0.4*zp_fcqfid)
                
                fnu_microJy_unc[this_fcqfid] = flux_dn_unc*10**(29 -
                                                                48.6/2.5 -
                                                                0.4*zp_fcqfid)  
                if deprecated:
                    n_base_obs[this_fcqfid] = n_baseline
                    which_base[this_fcqfid] = pre_or_post
                    C_baseline[this_fcqfid] = baseline
                else:
                    sys_sigma[this_fcqfid] = sys_unc

                    
    if write_lc is not False:

        if isinstance(write_lc, pd.DataFrame):
            write_df = write_lc
            write_df['jd'] = fp_df.jd.values
            write_df['forcediffimflux'] = fp_df.forcediffimflux.values
            write_df['forcediffimchisq'] =  fp_df.forcediffimchisq.values
        else:
            write_df = pd.DataFrame(fp_df.jd.values, columns=['jd'])
        write_df['fnu_microJy'] = fnu_microJy
        write_df['fnu_microJy_unc'] = fnu_microJy_unc
        write_df['passband'] = fp_df['filter'].values
        write_df['programid'] = fp_df.programid.values
        write_df['fcqfid'] = fp_df.fcqfid.values
        write_df['zpdiff'] = fp_df.zpdiff.values
        write_df['sys_unc_factor'] = sys_sigma
        write_df['flags'] = fp_df['flags'].values
        if deprecated:
            write_df['C'] = C_baseline
            write_df['N_baseline'] = n_base_obs
            write_df['pre_or_post'] = which_base

        gr_obs = np.where(write_df.fcqfid.values % 10 != 4)
        if not isinstance(write_lc, pd.DataFrame):
            fname = save_path + ztf_name + '_fnu.csv'
            write_df.iloc[gr_obs].to_csv(fname, index=False)
            
            storehdf = pd.HDFStore(save_path + ztf_name + '_fnu.h5')
            storehdf.put('light_curve', write_df)
            storehdf.get_storer('light_curve').attrs.metadata = fcqfid_dict
            storehdf.close()

    if make_plot is not False:

        color_dict = {1: 'MediumAquaMarine', 2: 'Crimson', 3: 'Goldenrod'}
        nplots = 0
        jdstart = 2458119.5
        for key in fcqfid_dict:
            if fcqfid_dict[key]['N_bl'] > 1:
                this_fcqfid_good = np.where((fp_df.fcqfid.values == int(key)) & 
                                            (bad_obs == 0))
                plot_flux = fnu_microJy[this_fcqfid_good]

                if not ((plot_flux == -999).sum() == len(plot_flux)):
                    if ((fcqfid_dict[key]['N_pre_peak'] > 2) or
                        (fcqfid_dict[key]['N_post_peak'] > 2)
                       ):
                        nplots += 1

        if nplots > 0:
            fig = plt.figure() if make_plot is True else make_plot
            fig.set_size_inches(8, nplots * 3 + 0.5)
            axes = fig.subplots(nplots, 1, sharex=True)
            plot_num = 0
            for key in fcqfid_dict:
                if fcqfid_dict[key]['N_bl'] > 1:
                    ufid = int(key)
                    this_fcqfid_good = np.where((fp_df.fcqfid.values == ufid) & 
                                                (bad_obs == 0))

                    plot_jd = fp_df.jd.values[this_fcqfid_good] - jdstart
                    plot_flux = fnu_microJy[this_fcqfid_good]
                    plot_flux_unc = fnu_microJy_unc[this_fcqfid_good]
                    
                    if (plot_flux == -999).sum() == len(plot_flux):
                        continue

                    ax = axes[plot_num] if nplots > 1 else axes
                    ax.errorbar(plot_jd, plot_flux, plot_flux_unc,
                                            fmt='o',
                                            mec=color_dict[ufid % 10],
                                            ecolor=color_dict[ufid % 10],
                                            mfc='None')
                    
                    if roll_med_plot == True:
                        jd_time = Time(plot_jd + jdstart, format='jd')
                        f_ser = pd.Series(plot_flux, 
                                          index = 
                                          pd.to_datetime(jd_time.datetime))
                        plot_roll = f_ser.rolling(window, 
                                                  center=True).median().values
                        ax.plot(plot_jd, plot_roll, 
                                color = 'lightgrey', zorder = 2)
                    
                    ax.axvline(x = t_peak - jdstart, color = '0.5', ls = '--')
                    ax.axhline(y = 0, color = '0.5',
                               ls = (0, [8, 1.5, 1, 1.5]), lw = 0.5, alpha=0.75)
                    ax.axvspan(0, t_peak - jdstart - 100,
                               color='Cornsilk', alpha=0.6, lw=0)
                    ax.axvspan(t_faded - jdstart, 1e6,
                               color='Cornsilk', alpha=0.6, lw=0)
                    ax.set_ylabel(r'flux ($\mu$Jy)', fontsize = 14)
                    ax.set_xlim(np.min(fp_df.jd.values - jdstart)-10,
                                np.max(fp_df.jd.values - jdstart)+10)
                    ax.set_title(f"{ztf_name}, {ufid}")
                    
                    if talk_plot:
                        ax.tick_params(axis='both', colors='white')
                        for spine in ['top', 'bottom', 'left', 'right']:
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
