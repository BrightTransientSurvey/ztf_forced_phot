import numpy as np
import pymc as pm
from pymc import Continuous, Model, Normal, Slice, sample, traceplot
from pymc.distributions import Interpolated
import arviz as az
import threading
from multiprocessing.pool import ThreadPool
import pandas as pd
import os
import matplotlib.pyplot as plt
from collections import defaultdict
from scipy.stats import median_abs_deviation

import argparse

core_count = 8


def calc_sn_exp_rise(t, A, B, t0, trise, offset):
    '''
    Calculate the rising portion of the SN model from Villar+19
    
    Parameters
    ----------
    t : array-like
        Times at which the model is evaluated
    
    A - float
        Constant of proportionality
    
    B - float
        constant of proportionality

    t0 - float
        reference time for exponential rise
    
    trise - float
        exponential decay factor
    
    offset - float
        scalar offset relative to 0 flux
    
    Returns
    -------
    f : array-like
        Eqn (1) from Villar+19 for rising flux of SN evaluated at
        times t
    
    References
    ----------
    (1) V. A. Villar, E. Berger, G. Miller et al. (2019) 
    Astrophysical Journal, 884, 83; doi: 10.3847/1538-4357/ab418c
    '''
    f = ((A + B*(t - t0))/(1 + np.exp(-(t - t0)/trise))) + offset
    return f
    
    
def calc_sn_exp_decline(t, A, B, t0, gamma, trise, tfall, offset):
    '''
    Calculate the declining portion of the SN model from Villar+19
    
    Parameters
    ----------
    t : array-like
        Times at which the model is evaluated
    
    A - float
        Constant of proportionality
    
    B - float
        constant of proportionality

    t0 - float
        reference time for exponential rise
    
    gamma - float
        time of plateau onset minus t0
    
    trise - float
        exponential decay factor
    
    tfall - float
        exponential decay factor controlling the decline
    
    offset - float
        scalar offset relative to 0 flux
    
    Returns
    -------
    f : array-like
        Eqn (1) from Villar+19 for rising flux of SN evaluated at
        times t
    
    References
    ----------
    (1) V. A. Villar, E. Berger, G. Miller et al. (2019) 
    Astrophysical Journal, 884, 83; doi: 10.3847/1538-4357/ab418c
    '''
    f = offset + ( (A + B*gamma) *
                   np.exp(-(t-(gamma+t0))/tfall) / 
                   (1 + np.exp(-(t-t0)/trise)))
    return f

def fit_gr(sn, lc_path='', out_path=''):
    '''
    Fit parametric model from Villar+19 to ZTF light curve
    
    Parameters
    ----------
    sn : string
        ZTF name of the SN to be fit by the model
    
    lc_path : string (optional, default = '')
        File path to the csv file with the ZTF light curve
    
    out_path : string (optional, default = '')
        File path to write output files (MCMC chains and summary statistics)
    '''
    
    lc_df = pd.read_csv(f"{lc_path}/{sn}_fnu.csv")
    for pb in lc_df.passband.unique():
        filt = pb[-1]
    
        # fit the g band model
        pb_lc = np.where((lc_df["passband"] == f'ZTF_{filt}') & 
                         (lc_df["flags"] == 0) & 
                         # # CHECK IF THIS IS WHAT SHOULD BE DONE!!!
                         # (lc_df["programid"] == 1) &
                         (lc_df["fnu_microJy"] != -999)
                        )
        if len(pb_lc[0]) < 10:
            continue
        jd0 = 2458119.5 # 2018 Jan 01
        time_axis = (lc_df['jd'].values)[pb_lc] - jd0

        Y_observed = ((lc_df['fnu_microJy']).values)[pb_lc]
        # Y_observed_filter.append(light_data)

        Y_unc = ((lc_df['fnu_microJy_unc']).values)[pb_lc]
    
        with pm.Model() as model:
            x_grid = np.logspace(-3, np.log10(60), 500)
            logistic_pdf = lambda x, tau, x0: 1/(1 + np.exp(-(x - x0)/tau)) 
            trise = pm.Interpolated('trise', x_grid, 
                                    logistic_pdf(x_grid, 5e-3, 0.1))
            tfall = pm.Uniform('tfall', lower = 1, upper = 300)
    
            flux_scatter = np.std(Y_observed, ddof=2)
            sigma_est = np.sqrt(np.mean(Y_unc**2))
    
            Amp_Guess = np.max(Y_observed)
            Amplitude = pm.TruncatedNormal('Amplitude', 
                                           mu=Amp_Guess,
                                           sigma=Amp_Guess/10, lower=0)
    
            Beta = pm.Uniform('Beta', 
                              lower = -np.max(Y_observed)/150, 
                              upper = 0)
    
            t_max = time_axis[np.argmax(Y_observed)] + 10
            t_min = np.min(time_axis) - 50
            t0 = pm.Uniform('t0', 
                            lower = time_axis[np.argmax(Y_observed)]-25,
                            upper = time_axis[np.argmax(Y_observed)]+25)
            
            scalar = pm.TruncatedNormal('scalar', 
                                        mu = 0, sigma = sigma_est, 
                                        lower=-2*sigma_est, 
                                        upper=2*sigma_est)

            no_p = pm.Normal.dist(mu = 5, sigma = 5)
            yes_p = pm.Normal.dist(mu = 60, sigma = 30)

            gamma = pm.NormalMixture("gamma", 
                                     w=[2/3, 1/3], 
                                     mu=[5, 60], 
                                     sigma=[5,30])

            # Expected value of outcome
            mu_rise = calc_sn_exp_rise(time_axis, 
                                       Amplitude, 
                                       Beta, 
                                       t0, 
                                       trise, 
                                       scalar)
            mu_fall = calc_sn_exp_decline(time_axis, 
                                          Amplitude, 
                                          Beta, 
                                          t0, 
                                          gamma, 
                                          trise, 
                                          tfall, 
                                          scalar)

            mu_switch = pm.math.switch(gamma+t0 >= time_axis, 
                                       mu_rise, mu_fall)

            # Likelihood (sampling distribution) of observations
            Y_obs = pm.Normal('Y_obs', mu=mu_switch, 
                              sigma=Y_unc, observed=Y_observed)

            data = pm.sample(1000, tune = 15000, cores = 4, 
                             return_inferencedata=True, 
                             progressbar = True, 
                             target_accept=0.99)
        data.to_netcdf(f"{out_path}/{sn}_{filt}.nc")
        out_summary = f"{out_path}/{sn}_{filt}.csv"
        az.summary(data,
                   stat_focus='median').to_csv(out_summary)

def plot_posterior_draws(sn, lc_path='', out_path='', save_fig=True):
    '''
    Plot posterior draws of model from Villar+19 to ZTF light curve
    
    Parameters
    ----------
    sn : string
        ZTF name of the SN to be fit by the model
    
    lc_path : string (optional, default = '')
        File path to the csv file with the ZTF light curve
    
    out_path : string (optional, default = '')
        File path to write output files (MCMC chains and summary statistics)
    
    save_fig : boolean (optional, default = True)
        Boolean flag indicating whether or not to save the plot as a png
    '''
    color_dict = {'ZTF_g': "MediumAquaMarine", 
                  'ZTF_r': "Crimson",
                  'ZTF_i': "GoldenRod"}

    lc_df = pd.read_csv(f"{lc_path}/{sn}_fnu.csv")
    fig, ax = plt.subplots(figsize=(10,4))
    for pb in lc_df.passband.unique():
        # this_pb = np.where(lc_df.passband == pb)
        filt = pb[-1]
        this_pb = np.where((lc_df["passband"] == f'ZTF_{filt}') & 
                         (lc_df["flags"] == 0) & 
                         # # CHECK IF THIS IS WHAT SHOULD BE DONE!!!
                         # (lc_df["programid"] == 1) &
                         (lc_df["fnu_microJy"] != -999)
                        )
        if len(this_pb[0]) < 10:
            continue
        
        jd0 = 2458119.5 # 2018 Jan 01
        
        chains = az.from_netcdf(f"{out_path}/{sn}_{filt}.nc")
        ax.errorbar(lc_df.jd.values[this_pb] - jd0, 
                    lc_df.fnu_microJy.values[this_pb], 
                    lc_df.fnu_microJy_unc.values[this_pb], 
                    fmt='o', color=color_dict[pb])    

        # max posterior plot
        pi_max_index = np.argmax(chains.sample_stats.lp.values.flatten())
        pi_max_t0 = chains.posterior.t0.values.flatten()[pi_max_index]
        pi_max_amp = chains.posterior.Amplitude.values.flatten()[pi_max_index]
        pi_max_beta = chains.posterior.Beta.values.flatten()[pi_max_index]
        pi_max_gamma = chains.posterior.gamma.values.flatten()[pi_max_index]
        pi_max_trise = chains.posterior.trise.values.flatten()[pi_max_index]
        pi_max_tfall = chains.posterior.tfall.values.flatten()[pi_max_index]
        pi_max_scalar = chains.posterior.scalar.values.flatten()[pi_max_index]

        t_grid_rise = np.linspace(pi_max_t0 - 150,
                                  pi_max_t0+pi_max_gamma,
                                  num = 10000)
        ax.plot(t_grid_rise, 
                calc_sn_exp_rise(t_grid_rise, 
                                 pi_max_amp, 
                                 pi_max_beta, 
                                 pi_max_t0, 
                                 pi_max_trise, 
                                 pi_max_scalar), 
                color = color_dict[pb])
        t_grid_decline = np.linspace(pi_max_t0+pi_max_gamma,
                                     np.max(lc_df.jd.values[this_pb]) - jd0,
                                     num = 10000)
        ax.plot(t_grid_decline, 
                calc_sn_exp_decline(t_grid_decline, 
                                    pi_max_amp, 
                                    pi_max_beta, 
                                    pi_max_t0, 
                                    pi_max_gamma, 
                                    pi_max_trise, 
                                    pi_max_tfall, 
                                    pi_max_scalar), 
                color = color_dict[pb])


        # posterior samples
        n_samples = len(chains.posterior.t0.values.flatten())
        rand_idx = np.random.choice(range(n_samples), 
                                    10, replace=False)
        pi_t0 = chains.posterior.t0.values.flatten()[rand_idx]
        pi_amp = chains.posterior.Amplitude.values.flatten()[rand_idx]
        pi_beta = chains.posterior.Beta.values.flatten()[rand_idx]
        pi_gamma = chains.posterior.gamma.values.flatten()[rand_idx]
        pi_trise = chains.posterior.trise.values.flatten()[rand_idx]
        pi_tfall = chains.posterior.tfall.values.flatten()[rand_idx]
        pi_scalar = chains.posterior.scalar.values.flatten()[rand_idx]

        t_grid_rise = np.linspace(pi_t0 - 150,
                                  pi_t0+pi_gamma,
                                  num = 10000)
        ax.plot(t_grid_rise, 
                calc_sn_exp_rise(t_grid_rise, 
                                 pi_amp, 
                                 pi_beta, 
                                 pi_t0, 
                                 pi_trise, 
                                 pi_scalar), 
                color = color_dict[pb], ls='--', lw=0.6, alpha=0.3)
        t_grid_decline = np.linspace(pi_t0+pi_gamma,
                                     np.max(lc_df.jd.values[this_pb]) - jd0,
                                     num = 10000)
        ax.plot(t_grid_decline, 
                calc_sn_exp_decline(t_grid_decline, 
                                    pi_amp, 
                                    pi_beta, 
                                    pi_t0, 
                                    pi_gamma, 
                                    pi_trise, 
                                    pi_tfall, 
                                    pi_scalar), 
                color = color_dict[pb], ls='--', lw=0.6, alpha=0.3)
    
        x_max = np.min([pi_max_t0 + pi_max_gamma + 10*pi_max_tfall, 
                       np.max(lc_df.jd.values[this_pb]) - jd0 + 10])
        ax.set_xlim(pi_max_t0 - 75, x_max)
        ax.set_ylim(-3*median_abs_deviation(lc_df.fnu_microJy.values[this_pb]), 
                    1.2*np.percentile(lc_df.fnu_microJy.values, 99.5))
        ax.set_xlabel('Time (JD - 2018 Jan 01)',fontsize=14)
        ax.set_ylabel(r'Flux ($\mu$Jy)',fontsize=14)
        ax.tick_params(axis='both', which='major', labelsize=12)
        fig.subplots_adjust(left=0.0.8,bottom=0.13,right=0.99, top=0.99)
    if save_fig:
        fig.savefig(f"{lc_path}/{sn}_posterior.png", 
                    dpi = 600, transparent=True)
        

def main():
    
    # Initialize argument parser
    parser = argparse.ArgumentParser(prog='fit_villar.py <sn>',
                                     description='Run Villar+19 light curve fits on a BTS fps lc.')
    # Necessary arguments
    parser.add_argument('sn', type=str, nargs='?', default=None,
                   help='ZTF transient name')
    # Optional arguments
    parser.add_argument('lc_path', type=str, nargs='?', default=None,
                   help='path to the processed fps lc file from Miller+24')
    parser.add_argument('out_path', type=str, nargs='?', default=None,
                   help='path for output MCMC chains')
    
    
    try:
        args = parser.parse_args()
        
        run = True
    except Exception:
        run = False
        
    if run:
        fit_gr(**vars(args))
    
if __name__ == "__main__":
    main()