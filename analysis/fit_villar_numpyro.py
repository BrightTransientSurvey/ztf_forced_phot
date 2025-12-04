import jax
import jax.numpy as jnp
import numpyro
from numpyro import distributions as dist, infer

import matplotlib.pyplot as plt
import numpy as np
import arviz as az
import pandas as pd
from scipy.stats import median_abs_deviation, truncnorm, rv_continuous, norm
from numpyro.infer.initialization import init_to_mean, init_to_uniform, init_to_value, init_to_feasible, init_to_median, init_to_sample

import argparse
import glob
import os



numpyro.set_host_device_count(4)


def y_model_villar19(t, A, B, t0, gamma, trise, tfall, offset):
    """
    Calculate the SN model from Villar+19. Incompatible with ztf priors from Superphot+

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
        Eqn (1) from Villar+19 for flux of SN evaluated at
        times t

    References
    ----------
    (1) V. A. Villar, E. Berger, G. Miller et al. (2019)
    Astrophysical Journal, 884, 83; doi: 10.3847/1538-4357/ab418c
    """
    f = jnp.where(gamma + t0 >= t,
                  ((A + B * (t - t0)) / (1 + jnp.exp(-(t - t0) / trise))) + offset,
                  offset + ((A + B * gamma) *
                            jnp.exp(-(t - (gamma + t0)) / tfall) /
                            (1 + jnp.exp(-(t - t0) / trise))))
    return f

def y_model_superphot(t, A, B, t0, gamma, trise, tfall, offset):
    """
    Calculate the SN model from de Soto+24

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
        Eqn (1) from deSoto+24 for flux of SN evaluated at
        times t

    References
    ----------
    (1) K. M. de Soto, A. Villar, E. Berger et al. (2024)
    Astrophysical Journal, 974, 169; doi: 10.3847/1538-4357/ad6a4f
    """
    f = jnp.where(t - t0 < gamma,
                  
                 ((A * (1 - (B * (t - t0)))) / (1 + jnp.exp(-(t - t0) / trise))) + offset,
                 
                  offset + (
                      (A * (1 - (B * gamma))) *
                            jnp.exp((gamma - (t - t0)) / tfall)
                      /
                            (1 + jnp.exp(-(t - t0) / trise))
                           ))
    return f

def load_lc_df(sn, lc_path, min_num_obs=2):#10): #TODO: REVERT THIS
    """
    Load ZTF light curve (LC) into a dataframe, with quality cuts

    Parameters
    ----------
    sn : string
        ZTF name of the SN to be fit by the model

    lc_path : string (optional, default = '')
        File path to folder containing processed fps lc files from Miller+24

    min_num_obs : int (optional, default = 10)
        If there are fewer than min_num_obs observations for any given filter,
        exclude from output dataframe

    Returns
    -------
    lc_df_clean : pandas dataframe
        Dataframe of LC obs. with flags==0 and valid flux measurements. Entire
        passbands may be excluded if there were fewer than min_num_obs valid points
    """

    # light curve data frame for all observations
    lc_df = pd.read_csv(f"{lc_path}/{sn}_fnu.csv")
    
    # make quality cuts
    keep_ind = []
    for pb in lc_df.passband.unique(): # loop over passbands
        pb_lc = np.where((lc_df["passband"] == pb) &
                         (lc_df["flags"] == 0) &
                         # # CHECK IF THIS IS WHAT SHOULD BE DONE!!!
                         # (lc_df["programid"] == 1) &
                         (lc_df["fnu_microJy"] != -999)
                         )
        if len(pb_lc[0]) < min_num_obs:
            # if fewer than min_num_obs pts in this filter, don't include in df
            continue
        else:
            keep_ind.append(pb_lc[0])
    if len(keep_ind) == 0:
        print(f"{sn} has 0 data points passing quality cuts; skipping this source")
    lc_df_clean = lc_df.iloc[np.concatenate(keep_ind)]
    
    return lc_df_clean

def villar_fit_constraint(x):
    """
    Apply constraints described in Appendix B of de Soto+24
    1. ensure piecewise transition happens after plateau
    2. at t - t0 = gamma, enforce negative (or 0) derivative of first piecewise portion 
    3. at t - t0 = gamma, ensure derivative of second piecewise portion is more negative than that of first portion
    4. (via 1-3) model flux is always positive

    Parameters
    ----------
    x : list
        list containing beta, gamma, trise, tfall

    Returns
    -------
     : tuple
        expressions (B2) and (B3) from de Soto+24

    References
    ----------
    (1) K. M. de Soto, A. Villar, E. Berger et al. (2024)
    Astrophysical Journal, 974, 169; doi: 10.3847/1538-4357/ad6a4f
    (2) https://github.com/VTDA-Group/superphot-plus
    """

    #amplitude, beta, gamma, t0, trise, tfall = x
    
    # return ( 
    #     jnp.maximum(gamma * beta - 1., 0.) + ##<-
    #     jnp.maximum(jnp.exp(-gamma/trise)*(tfall/trise-1) - 1., 0.) + ##<-
    #     jnp.maximum(beta*tfall - 1. + beta*gamma, 0.) ##<- same as eq B3 (ensure deriv. of second piecewise portion is more negative than that of first at t-t0=gamma)
    #     #jnp.maximum(jnp.exp(-gamma / trise) * (1.0 / beta - trise - gamma) - trise, 0.) ##<- same as eq B2, and it's not getting used!! (enforce transition happens after plateau)
    # ) <- TRANSLATED newly defined on superphot+ github 2025


    # return (
    #     jnp.maximum(x[2] * x[1] - 1., 0.) +
    #     jnp.maximum(jnp.exp(-x[2]/x[4])*(x[5]/x[4]-1) - 1., 0.) +
    #     jnp.maximum(x[1]*x[5] - 1. + x[1]*x[2], 0.)
    #     #jnp.maximum(jnp.exp(-x[2] / x[4]) * (1.0 / x[1] - x[4] - x[2]) - x[4], 0.)
    # ) <- newly defined on superphot+ github 2025

    # -------------------------------
    beta, gamma, trise, tau_fall = x

    return (
        jnp.maximum(gamma - (1.0 - beta * tau_fall) / beta, 0.) +  #<- same as eq B3
        jnp.maximum(jnp.exp(-gamma / trise) * (1.0/beta - trise - gamma) - trise, 0.)
    ) # <- as once defined on superphot+ 2024

def lc_model_superphot(t_obs, Y_unc, Y_obs=None, prior_vars='old'):
    """
    Model for ZTF light curves, to be handed to `fit_gr_numpyro` -> `numpyro.infer.NUTS`. 
    Uses priors based on superphot_plus/samplers/numpyro_sampler.py

    Parameters
    ----------
    t_obs : array-like
        Time values at which the model is evaluated

    Y_unc : array-like
        Flux uncertainties at the observed times t_obs

    Y_obs : array-like (optional, default = None)
        Flux observations at times t_obs

    References
    ----------
    (1) https://github.com/VTDA-Group/superphot-plus
    """

    # define time of maximum, t0 prior depends on it
    tFmax = jnp.array(t_obs)[jnp.argmax(Y_obs)]

    # define priors variables from superphot_plus/surveys/ztf.yaml (reference, works for r-band)
    old_prior_vars = {'log_amplitude': (0.0957, 0.0575,-0.3, 0.5), 
            'beta': (0.00833, 0.00385, 0, 0.03),
            't0': (tFmax - 17.878, 9.916, tFmax - 100, tFmax + 30), 
            'log_gamma': (1.4258, 0.3079, 0, 3.5),
            'log_trise': (0.6664, 0.4250, -2.0, 4.0),
            'log_tfall': (1.5261, 0.3037, 0, 4), 
            'log_extra_sigma': (-1.6629, 0.3378, -3, -0.8)}
    # define priors variables from RECENTLY CHANGED, WIDER superphot_plus/surveys/ztf.yaml (reference, works for r-band)
    new_prior_vars = {'log_amplitude': (0.0957, 0.15,-0.3, 0.5), 
            'beta': (0.00833, 0.012, -0.01, 0.03), 
            't0': (tFmax - 17.878, 30, tFmax - 100, tFmax + 30), 
            'log_gamma': (1.4258, 0.9, 0, 3.5),
            'log_trise': (0.6664, 1.2, -2.0, 4.0),
            'log_tfall': (1.5261, 0.9, 0, 4), 
            'log_extra_sigma': (-1.6629, 0.9, -3, -0.8)}
    # choose which set to use
    if prior_vars == 'old':
        prior_vars = old_prior_vars
        print('using old superphot+ priors (narrow)')
    elif prior_vars == 'new':
        prior_vars = new_prior_vars
        print('using new superphot+ priors (wide)')
    else:
        print('undefined prior variables, quitting')
        return

    # define and un-log10- the priors where applicable
    amplitude = 10 ** numpyro.sample("log_amplitude", dist.TruncatedNormal(loc=prior_vars["log_amplitude"][0],
                                                                    scale=prior_vars["log_amplitude"][1],
                                                                    low=prior_vars["log_amplitude"][2], 
                                                                    high=prior_vars["log_amplitude"][3], validate_args=False))

    beta = numpyro.sample("beta", dist.TruncatedNormal(loc=prior_vars["beta"][0],
                                                        scale=prior_vars["beta"][1],
                                                        low=prior_vars["beta"][2],
                                                        high=prior_vars["beta"][3], validate_args=False))

    gamma = 10 ** numpyro.sample("log_gamma", dist.TruncatedNormal(loc=prior_vars["log_gamma"][0], 
                                                                   scale=prior_vars["log_gamma"][1],
                                                                   low=prior_vars["log_gamma"][2], 
                                                                   high=prior_vars["log_gamma"][3], validate_args=False))
    
    t0 = numpyro.sample("t0", dist.TruncatedNormal(loc=prior_vars["t0"][0],
                                                    scale=prior_vars["t0"][1],
                                                    low=prior_vars["t0"][2], 
                                                    high=prior_vars["t0"][3], validate_args=False))

    trise = 10 ** numpyro.sample("log_trise", dist.TruncatedNormal(loc=prior_vars["log_trise"][0], 
                                                                    scale=prior_vars["log_trise"][1],
                                                                    low=prior_vars["log_trise"][2], 
                                                                    high=prior_vars["log_trise"][3],
                                                                      validate_args=False))

    tfall = 10 ** numpyro.sample("log_tfall", dist.TruncatedNormal(loc=prior_vars["log_tfall"][0],
                                                                    scale=prior_vars["log_tfall"][1],
                                                                    low=prior_vars["log_tfall"][2], 
                                                                    high=prior_vars["log_tfall"][3], validate_args=False))

    extra_sigma = 10 ** numpyro.sample("log_extra_sigma", dist.TruncatedNormal(loc=prior_vars["log_extra_sigma"][0],
                                                                                scale=prior_vars["log_extra_sigma"][1],
                                                                                low=prior_vars["log_extra_sigma"][2], 
                                                                                high=prior_vars["log_extra_sigma"][3], validate_args=False))
    
    # # scalar term was in Villar+19, although not in superphot+ models
    # sigma_est = jnp.sqrt(jnp.mean(Y_unc ** 2))
    # scalar = numpyro.sample("scalar", dist.TruncatedNormal(loc=0, scale=sigma_est,
    #                                                       low=-2 * sigma_est, high=2 * sigma_est))

    # incorporate constraint from de Soto+24 Appendix B, this adds an arbitrary log prob. factor
    #constraint = villar_fit_constraint([beta, gamma, trise, tfall])
    #numpyro.factor(
    #    "vf_constraint",
    #    -1000. * jnp.max(constraint)
    #)

    # define uncertainties (add extra_sigma to observed uncertainties in quad.)
    sigma_tot = jnp.sqrt(Y_unc**2 + extra_sigma**2)

    # expected value of outcome
    mu_switch = y_model_superphot(t_obs, amplitude, beta, t0, gamma, trise, tfall, 
                                 0)#,scalar)

    # Sample!
    numpyro.sample("y",
                   dist.Normal(mu_switch, sigma_tot),
                   obs=Y_obs)

def lc_model_villar19_miller(t_obs, Y_unc, Y_obs=None):
    """
    model for ZTF light curves, with priors loosely based on Villar+19, and Adam's pymc scripts

    Parameters
    ----------
    t_obs : array-like
        Time values at which the model is evaluated

    Y_unc : array-like
        Flux uncertainties at the observed times t_obs

    Y_obs : array-like (optional, default = None)
        Flux observations at times t_obs
    """

    # Define priors based on Villar+19
    trise = numpyro.sample("trise", dist.Uniform(low=0.01, high=50))  # dist.continuous.Uniform?
    #trise = numpyro.sample("trise", dist.TruncatedNormal(loc=6, scale=2, low=0))
    tfall = numpyro.sample("tfall", dist.Uniform(low=1, high=300))

    Amp_Guess = jnp.max(Y_obs)
    amplitude = numpyro.sample("amplitude", dist.TruncatedNormal(
        loc=Amp_Guess,
        scale=Amp_Guess / 10, low=0))

    beta = numpyro.sample("beta", dist.Uniform(low=-jnp.max(Y_obs) / 150, high=0))
    #uniformbeta_low = -jnp.max(Y_obs) / 150
    #beta = numpyro.sample("beta", dist.Normal(0, uniformbeta_low * (-0.3)))

    t0 = numpyro.sample("t0", dist.Uniform(
        low=jnp.array(t_obs)[jnp.argmax(Y_obs)] - 25,
        high=jnp.array(t_obs)[jnp.argmax(Y_obs)] + 25))
    #t0_uniform_prior_center = jnp.array(t_obs)[jnp.argmax(Y_obs)]
    #t0 = numpyro.sample("t0", dist.Normal(t0_uniform_prior_center - 10, 10))

    sigma_est = jnp.sqrt(jnp.mean(Y_unc ** 2))
    scalar = numpyro.sample("scalar", dist.TruncatedNormal(loc=0, scale=sigma_est,
                                                           low=-2 * sigma_est,
                                                           high=2 * sigma_est))

    ## gamma's prior is a normal mixture - not straightforward in numpyro
    # Define the weights for the mixture components
    weights = jnp.array([2 / 3, 1 / 3])  # Should sum to 1
    # Define means and standard deviations for the normal components
    means, stds = jnp.array([5, 60]), jnp.array([5, 30])
    # Define the normal distributions
    components = [dist.Normal(mu, sigma) for mu, sigma in zip(means, stds)]
    # Define the mixture distribution
    mixture = dist.Mixture(
        dist.Categorical(probs=weights),  # Mixture weights
        components  # List of component distributions
    )
    gamma = numpyro.sample("gamma", mixture)


    # Expected value of outcome
    mu_switch = y_model_villar19(t_obs, amplitude, beta, t0, gamma, trise, tfall, scalar)

    # Sample!
    numpyro.sample("y",
                   dist.Normal(mu_switch, Y_unc),
                   obs=Y_obs)

def fit_gr_numpyro(sn, lc_path, out_path, num_warmup=15000, num_samples=1000, num_chains=4, init_strat='uniform', model=lc_model_superphot, randomkey=0, jd0 = 2458119.5):
    """
    Fit parametric model from [Villar+19 OR de Soto+24/Superphot+] to all light curves (in available filters) from a single object `sn`.

    Parameters
    ----------
    sn : string
        ZTF name of the SN to be fit by the model

    lc_path : string (optional, default = '')
        File path to folder containing processed fps lc files from Miller+24

    out_path : string (optional, default = '')
        File path to write output files (MCMC chains and summary statistics)

    num_warmup : int (optional, default = 15000)
        Number of warmup steps to run MCMC chains for

    num_samples : int (optional, default = 1000)
        Number of samples to run MCMC chains for

    num_chains : int (optional, default = 4)
        Number of chains to run MCMC, make sure num_chains =< jax.local_device_count()
        ( you can change device count using: numpyro.set_host_device_count(#) )

    init_strat : str (optional, default = 'uniform')
        Initialization strategy for numpyro, can be 'uniform', 'mean', or 'value'.
        See https://num.pyro.ai/en/stable/utilities.html#init-strategy for details.

    model : function (optional, default = lc_model_superphot)
        numpyro modeling function from Villar19+ (lc_model_villar19_miller) OR non-hierarchical superphot+ (lc_model_superphot)
    """

    # load light curve dataframe with some quality cuts
    lc_df = load_lc_df(sn=sn, lc_path=lc_path)

    for pb in lc_df.passband.unique(): # loop over available passbands
        # preprocessing data
        filt = pb[-1]
        ## fit the model with the current filter
        lc_df_thisfilt = lc_df[(lc_df["passband"] == f'ZTF_{filt}')]
        ## shift obs. times relative to 2018 Jan 01
        t_obs_shifted = (lc_df_thisfilt['jd'].values) - jd0
        ## scale flux and flux uncertainty down by max. observed flux
        Y_obs_orig = ((lc_df_thisfilt['fnu_microJy']).values)
        Y_obs_max = jnp.max(Y_obs_orig)
        print(f'scaling data by dividing by Y_observed_max: {Y_obs_max}. keep in mind when visualizing posterior draws!!!')
        np.save(f'{out_path}/{sn}_{filt}_maxflux_{Y_obs_max}.npy', Y_obs_max) # save Y_obs_max, it's important for later!!!

        Y_obs_scaled = Y_obs_orig / Y_obs_max
        Y_unc_scaled = ((lc_df_thisfilt['fnu_microJy_unc']).values) / Y_obs_max

        # initalization strategy
        ## map input "init_strat" string to define strategy
        init_strat_dict = {'uniform' : init_to_uniform,
                       'mean': init_to_mean,
                       'feasible': init_to_feasible,
                       'median': init_to_median,
                       'sample': init_to_sample,
                       'value': "need init. values"}
        this_init_strat = init_strat_dict[init_strat]
        ## specify initalization values (hardcoded for now, should only be used for testing)
        if (this_init_strat == "need init. values") & (model == lc_model_superphot):
            print('warning: init_to_value method is hardcoded')
            # amp	beta	t0	gamma	trise	tfall
            #truths = [1.196678, 0.013459, 1382.580604, 60.483029, 5.535795, 22.128788]
            truths = [0.116, 0.012, 1.548, 1409.233, 0.415, 1.519, -1.66] # LOGGED APPROPRIATELY, in order of table in diagnos. plot
            ## note t0 and log_gamma listed in swapped order above
            this_init_strat = init_to_value(values = {"log_amplitude": truths[0], "beta": truths[1], 
                                                    "t0": truths[3], "log_gamma": truths[2],
                                                    "log_trise": truths[4], 
                                                    "log_tfall": truths[5], "scalar": truths[6]})
        elif init_strat not in ['uniform', 'mean', 'value', 'feasible', 'median', 'sample']:
            print('could not define initialization strategy based on inputs!')
            return

        # numpyro setup
        sampler = infer.MCMC(
            infer.NUTS(model, init_strategy = this_init_strat),
            num_warmup=num_warmup, num_samples=num_samples, num_chains=num_chains,
            progress_bar=True, chain_method="parallel")
        ## collecting warmup sample chains (for testing)
        sampler.warmup(jax.random.PRNGKey(randomkey), t_obs_shifted, Y_unc_scaled, Y_obs=Y_obs_scaled, collect_warmup=True)
        warmup_samples = sampler.get_samples()
        np.save(f'{out_path}/{sn}_{filt}_warmupsamples.npy', warmup_samples)

        ## draw samples from the posterior
        sampler.run(jax.random.PRNGKey(randomkey), t_obs_shifted, Y_unc_scaled, Y_obs=Y_obs_scaled)
        ## collecting sampling chains (for testing)
        #samples = sampler.get_samples()
        #np.save(f'{out_path}/{sn}_{filt}_samples.npy', samples)

        ## save results
        data = az.from_numpyro(sampler)
        data.to_netcdf(f"{out_path}/{sn}_{filt}_numpyro.nc")
        out_summary = f"{out_path}/{sn}_{filt}_numpyro.csv"
        az.summary(data, stat_focus='median').to_csv(out_summary)

# --- --- functions for visualizing results from Villar19+ fits (does not work for Superphot+ fits)

class MixtureModel(rv_continuous):
    def __init__(self, submodels, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.submodels = submodels

    def _pdf(self, x):
        pdf = self.submodels[0].pdf(x)
        for submodel in self.submodels[1:]:
            pdf += submodel.pdf(x)
        pdf /= len(self.submodels)
        return pdf

    def rvs(self, size):
        submodel_choices = np.random.randint(len(self.submodels), size=size)
        submodel_samples = [submodel.rvs(size=size) for submodel in self.submodels]
        rvs = np.choose(submodel_choices, submodel_samples)
        return rvs

def posterior_from_chains(chains_filt, y_obs_filt, y_unc_filt):
    """
    for one filter, calculate posterior (at every step in every chain)
    """
    # posterior = prior * likelihood
    ## grab likelihoods from chains, define priors from scratch
    
    # define log likelihood for the (4) chains
    summed = np.sum(chains_filt.log_likelihood.y.values, axis=2)

    # define all ~gaussian priors (ignore uniform)

    ## scalar (truncnorm)
    sigma_est = np.sqrt(np.mean(y_unc_filt**2))
    myclip_low = -2* sigma_est
    myclip_high = 2* sigma_est
    low, high = (myclip_low - 0) / sigma_est, (myclip_high - 0) / sigma_est
    rv = truncnorm(low, high, loc=0, scale=sigma_est)
    ln_scalar = np.log(rv.pdf(chains_filt.posterior.scalar))

    ## Amplitude (truncnorm)
    amp_guess = np.max(y_obs_filt)
    myclip_low_amp = 0
    myclip_high_amp = 3*amp_guess
    low_amp, high_amp = (myclip_low_amp - amp_guess) / (amp_guess/10), (myclip_high_amp - amp_guess) / (amp_guess/10)
    rv_2 = truncnorm(low_amp, high_amp, loc=amp_guess, scale=amp_guess/10)
    ln_amp = np.log(rv_2.pdf(chains_filt.posterior.Amplitude))

    ## Gamma (normal mixture)
    mixture_gaussian_model = MixtureModel([norm(5, 5), norm(5, 5), norm(60, 30)]) # cheat the 2/3 and 1/3 weights by mixing 3 normals
    ln_gamma = np.log(mixture_gaussian_model.pdf(chains_filt.posterior.gamma))

    # posterior calc 
    ln_posterior = summed + ln_scalar + ln_amp + ln_gamma
    return ln_posterior

def plot_posterior_draws_numpyro(sn, lc_path='', out_path='', save_fig=True, jd0 = 2458119.5):
    """
    Plot posterior draws of model from Villar+19 to ZTF light curve

    Parameters
    ----------
    sn : string
        ZTF name of the SN to be fit by the model

    lc_path : string (optional, default = '')
        File path to folder containing processed fps lc files from Miller+24

    out_path : string (optional, default = '')
        File path to files saved by fit_gr()  (MCMC chains and summary statistics)

    save_fig : boolean (optional, default = True)
        Boolean flag indicating whether or not to save the plot as a png
    """
    color_dict = {'ZTF_g': "MediumAquaMarine",
                  'ZTF_r': "Crimson",
                  'ZTF_i': "GoldenRod"}

    lc_df = load_lc_df(sn, lc_path)
    fig, ax = plt.subplots(figsize=(10, 4))
    xlims, ylims = [], []
    for pb in lc_df.passband.unique():
        filt = [-1]

        lc_df_thisfilt = lc_df[(lc_df["passband"] == f'ZTF_{filt}')]

        try:
            chains = az.from_netcdf(f"{out_path}/{sn}_{filt}_numpyro.nc")
        except:
            print(f'Unable to find chains for {sn} in {pb}, skipping this filter.')
            continue
        ax.errorbar(lc_df_thisfilt.jd.values - jd0,
                    lc_df_thisfilt.fnu_microJy.values,
                    lc_df_thisfilt.fnu_microJy_unc.values,
                    fmt='o', color=color_dict[pb])

        # max posterior plot
        post = posterior_from_chains(chains,
                                     lc_df_thisfilt['fnu_microJy'],
                                     lc_df_thisfilt['fnu_microJy_unc'])
        pi_max_index = np.argmax(post.flatten())
        pi_max_t0 = chains.posterior.t0.values.flatten()[pi_max_index]
        pi_max_amp = chains.posterior.amplitude.values.flatten()[pi_max_index]
        pi_max_beta = chains.posterior.beta.values.flatten()[pi_max_index]
        pi_max_gamma = chains.posterior.gamma.values.flatten()[pi_max_index]
        pi_max_trise = chains.posterior.trise.values.flatten()[pi_max_index]
        pi_max_tfall = chains.posterior.tfall.values.flatten()[pi_max_index]
        pi_max_scalar = chains.posterior.scalar.values.flatten()[pi_max_index]

        t_grid = jnp.linspace(pi_max_t0 - 150,
                              jnp.max(lc_df_thisfilt.jd.values) - jd0,
                              num=20000)
        ax.plot(t_grid,
                y_model_villar19(t_grid,
                                 pi_max_amp,
                                 pi_max_beta,
                                 pi_max_t0,
                                 pi_max_gamma,
                                 pi_max_trise,
                                 pi_max_tfall,
                                 pi_max_scalar),
                color=color_dict[pb])

        # posterior samples
        n_samples = len(chains.posterior.t0.values.flatten())
        rand_idx = np.random.choice(range(n_samples),
                                    10, replace=False)
        pi_t0 = chains.posterior.t0.values.flatten()[rand_idx]
        pi_amp = chains.posterior.amplitude.values.flatten()[rand_idx]
        pi_beta = chains.posterior.beta.values.flatten()[rand_idx]
        pi_gamma = chains.posterior.gamma.values.flatten()[rand_idx]
        pi_trise = chains.posterior.trise.values.flatten()[rand_idx]
        pi_tfall = chains.posterior.tfall.values.flatten()[rand_idx]
        pi_scalar = chains.posterior.scalar.values.flatten()[rand_idx]

        t_grid = jnp.linspace(pi_t0 - 150,
                              np.max(lc_df_thisfilt.jd.values) - jd0,
                              num=20000)
        ax.plot(t_grid,
                y_model_villar19(t_grid,
                                 pi_amp,
                                 pi_beta,
                                 pi_t0,
                                 pi_gamma,
                                 pi_trise,
                                 pi_tfall,
                                 pi_scalar),
                color=color_dict[pb], ls='--', lw=0.6, alpha=0.3)

        
        x_max = np.min([pi_max_t0 + pi_max_gamma + 10 * pi_max_tfall,
                        np.max(lc_df_thisfilt.jd.values) - jd0 + 10])
        #if filt != 'i': ## make sure to NOT set axis limits based on i-band
        # save axis limits based on this filter
        xlims.append((pi_max_t0 - 75, x_max))
        ylims.append((-3 * median_abs_deviation(lc_df_thisfilt.fnu_microJy.values),
                        1.5 * np.percentile(lc_df_thisfilt.fnu_microJy.values, 99.5)))

    # choose min and max of all lims
    ax.set_xlim(np.min([i[0] for i in xlims]), np.max([i[1] for i in xlims]))
    ax.set_ylim(np.min([i[0] for i in ylims]), np.max([i[1] for i in ylims]))

    # labels and save figure
    ax.set_xlabel('Time (JD - 2018 Jan 01)', fontsize=14)
    ax.set_ylabel(r'Flux ($\mu$Jy)', fontsize=14)
    ax.tick_params(axis='both', which='major', labelsize=12)
    if save_fig:
        fig.savefig(f"{lc_path}/{sn}_posterior_numpyro.png",
                    dpi=600, transparent=True)

# --- DIAGNOSTIC PLOTS ---
def lnl_from_sim(fake_lc, truth, max_flux, jd0 = 2458119.5):
    """
    knowing the true model param. values (from simulation), calculate lnl with sim. data and sim. truth model
    returns: a single value!
    """
    
    # simulated SN data
    data = fake_lc['fnu_microJy']
    N = len(data)
    #max_flux = np.max(fake_lc_smooth['fnu_microJy'])
    unc = fake_lc['fnu_microJy_unc']

    # model using the "true" params from which SN was simulated
    model = max_flux * y_model_superphot(fake_lc['jd'].values-jd0, truth['amp'][0], 
                                   truth['beta'][0], truth['t0'][0], 
                                   truth['gamma'][0], truth['trise'][0], truth['tfall'][0], 0)

    # log likelihood calc.
    lnl = (-1 * (N/2) * np.log(2*np.pi)) - (np.sum(np.log(unc))) -1 * np.sum((data - model)**2 / (2 * unc**2))
    return lnl

def lnl_from_chains(fake_lc, chains, max_flux, jd0 = 2458119.5):
    """
    making this because we think the log_likelihood var is incorrect. calc manually!!
    returns: list of len = nchains, each entry is a list of len = nsteps (aka one lnl value for every step in each chain)
    """

    # simulated SN data
    data = fake_lc['fnu_microJy']
    N = len(data)
    unc = fake_lc['fnu_microJy_unc']
    

    # access parameter values across nchains and over nsteps
    nchains = chains.posterior.log_amplitude.values.shape[0]
    nsteps = chains.posterior.log_amplitude.values.shape[1]

    lnl_chains = []

    for i in range(nchains): # loop over chains
        amp_step = 10 ** chains.posterior.log_amplitude.values[i]
        beta_step = chains.posterior.beta.values[i]
        gamma_step = 10 ** chains.posterior.log_gamma.values[i]
        t0_step = chains.posterior.t0.values[i]
        tfall_step = 10 ** chains.posterior.log_tfall.values[i]
        trise_step = 10 ** chains.posterior.log_trise.values[i]
        extra_sigma_step = 10 ** chains.posterior.log_extra_sigma.values[i]

        lnl_chain = []
        for j in range(nsteps):
            # model at every step of the chains
            model = max_flux * y_model_superphot(fake_lc['jd'].values-jd0, amp_step[j], 
                                    beta_step[j], t0_step[j],
                                    gamma_step[j], trise_step[j], tfall_step[j], 0)
            
            
            sigma_tot = np.sqrt(unc**2 + extra_sigma_step[j]**2)

            # log likelihood calc.
            lnl_step = (-1 * (N/2) * np.log(2*np.pi)) - (np.sum(np.log(sigma_tot))) -1 * np.sum((data - model)**2 / (2 * sigma_tot**2))
            lnl_chain.append(lnl_step)
        lnl_chains.append(lnl_chain)
    return lnl_chains

def get_truncated_normal(mean=0, sd=1, low=0, upp=10):
    '''redefines inputs, see note at: 
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.truncnorm.html'''
    return truncnorm(
        (low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)

def posterior_from_chains_superphot(data, chains, max_flux, prior_vers = 'old', jd0 = 2458119.5):
    """
    for one filter, calculate posterior (at every step in every chain)
    """

    y_obs_filt, y_unc_filt, t_obs = data['fnu_microJy'], data['fnu_microJy_unc'], data['jd'] - jd0
    tFmax = np.array(t_obs)[np.argmax(y_obs_filt)]

    old_prior_vars = {'log_amplitude': (0.0957, 0.0575,-0.3, 0.5), 
            'Beta': (0.00833, 0.00385, 0, 0.03),
            't0': (tFmax - 17.878, 9.916, tFmax - 100, tFmax + 30), 
            'log_gamma': (1.4258, 0.3079, 0, 3.5),
            'log_trise': (0.6664, 0.4250, -2.0, 4.0),
            'log_tfall': (1.5261, 0.3037, 0, 4), 
            'log_extra_sigma': (-1.6629, 0.3378, -3, -0.8)}

    # posterior = prior * likelihood
    ## define priors from scratch, grab likelihoods from chains
    
    # define log likelihood for the (4) chains
    ##summed = np.sum(chains_filt.log_likelihood.y.values, axis=2)
    lnl_all_steps = lnl_from_chains(data, chains, max_flux)
    print(np.array(lnl_all_steps).shape)
    #lnl_summed = np.sum(np.array(lnl_all_steps), axis=1)
    #print(lnl_summed.shape)

    
    # define all ~gaussian priors

    if prior_vers == 'old':
        p = old_prior_vars
        ## amplitude ( log truncnorm )

        log10_amp = get_truncated_normal(mean=p['log_amplitude'][0], sd=p['log_amplitude'][1], low=p['log_amplitude'][2], upp=p['log_amplitude'][3])
        ln_amp = np.log((log10_amp.pdf(chains.posterior.log_amplitude)))

        ## beta ( truncnorm )
        beta = get_truncated_normal(mean=p['Beta'][0], sd=p['Beta'][1], low=p['Beta'][2], upp=p['Beta'][3])
        ln_beta = np.log(beta.pdf(chains.posterior.beta))

        ## gamma ( log truncnorm )
        log10_gamma = get_truncated_normal(mean=p['log_gamma'][0], sd=p['log_gamma'][1], low=p['log_gamma'][2], upp=p['log_gamma'][3])
        ln_gamma = np.log((log10_gamma.pdf(chains.posterior.log_gamma)))
        
        ## t0 ( truncnorm )
        t0 = get_truncated_normal(mean=p['t0'][0], sd=p['t0'][1], low=p['t0'][2], upp=p['t0'][3])
        ln_t0 = np.log(t0.pdf(chains.posterior.t0))

        ## trise ( log truncnorm )
        log10_trise = get_truncated_normal(mean=p['log_trise'][0], sd=p['log_trise'][1], low=p['log_trise'][2], upp=p['log_trise'][3])
        ln_trise = np.log((log10_trise.pdf(chains.posterior.log_trise)))

        ## tfall ( log truncnorm )
        log10_tfall = get_truncated_normal(mean=p['log_tfall'][0], sd=p['log_tfall'][1], low=p['log_tfall'][2], upp=p['log_tfall'][3])
        ln_tfall = np.log((log10_tfall.pdf(chains.posterior.log_tfall)))

        ## extra_sigma ( log truncnorm )
        log10_extrasigma = get_truncated_normal(mean=p['log_extra_sigma'][0], sd=p['log_extra_sigma'][1], low=p['log_extra_sigma'][2], upp=p['log_extra_sigma'][3])
        ln_extrasigma = np.log(log10_extrasigma.pdf(chains.posterior.log_extra_sigma))

        # posterior calc 
        ln_posterior = lnl_all_steps + ln_amp + ln_beta + ln_gamma + ln_t0 + ln_trise + ln_tfall + ln_extrasigma
        
        return ln_posterior

    else:
        print('undefined prior variables')
        return

def posterior_from_sim_superphot(data, truths, max_flux, prior_vers = 'old', jd0 = 2458119.5):
    """
    for one filter, calculate posterior (at every step in every chain)
    """

    y_obs_filt, y_unc_filt, t_obs = data['fnu_microJy'], data['fnu_microJy_unc'], data['jd'] - jd0
    tFmax = np.array(t_obs)[np.argmax(y_obs_filt)]

    old_prior_vars = {'log_amplitude': (0.0957, 0.0575,-0.3, 0.5), 
            'Beta': (0.00833, 0.00385, 0, 0.03),
            't0': (tFmax - 17.878, 9.916, tFmax - 100, tFmax + 30), 
            'log_gamma': (1.4258, 0.3079, 0, 3.5),
            'log_trise': (0.6664, 0.4250, -2.0, 4.0),
            'log_tfall': (1.5261, 0.3037, 0, 4), 
            'log_extra_sigma': (-1.6629, 0.3378, -3, -0.8)}

    # define priors variables from RECENTLY CHANGED, WIDER superphot_plus/surveys/ztf.yaml (reference, works for r-band)
    new_prior_vars = {'log_amplitude': (0.0957, 0.15,-0.3, 0.5), 
            'beta': (0.00833, 0.012, -0.01, 0.03), 
            't0': (tFmax - 17.878, 30, tFmax - 100, tFmax + 30), 
            'log_gamma': (1.4258, 0.9, 0, 3.5),
            'log_trise': (0.6664, 1.2, -2.0, 4.0),
            'log_tfall': (1.5261, 0.9, 0, 4), 
            'log_extra_sigma': (-1.6629, 0.9, -3, -0.8)}

    # posterior = prior * likelihood
    ## define priors from scratch, grab likelihoods from chains
    
    # define log likelihood for the (4) chains
    ##summed = np.sum(chains_filt.log_likelihood.y.values, axis=2)
    lnl_summed = lnl_from_sim(data, truths, max_flux)

    
    # define all ~gaussian priors

    if prior_vers == 'old':
        p = old_prior_vars
    elif prior_vers == 'new':
        p = new_prior_vars
    else:
        print('undefined prior variables')
        return
    
    ## amplitude ( log truncnorm )

    log10_amp = get_truncated_normal(mean=p['log_amplitude'][0], sd=p['log_amplitude'][1], low=p['log_amplitude'][2], upp=p['log_amplitude'][3])
    ln_amp = np.log(log10_amp.pdf(np.log10(truths['amp'][0])))


    ## beta ( truncnorm )
    beta = get_truncated_normal(mean=p['Beta'][0], sd=p['Beta'][1], low=p['Beta'][2], upp=p['Beta'][3])
    ln_beta = np.log(beta.pdf(truths['beta'][0]))

    ## gamma ( log truncnorm )
    log10_gamma = get_truncated_normal(mean=p['log_gamma'][0], sd=p['log_gamma'][1], low=p['log_gamma'][2], upp=p['log_gamma'][3])
    ln_gamma = np.log(log10_gamma.pdf(np.log10(truths['gamma'])))
    
    ## t0 ( truncnorm )
    t0 = get_truncated_normal(mean=p['t0'][0], sd=p['t0'][1], low=p['t0'][2], upp=p['t0'][3])
    ln_t0 = np.log(t0.pdf(truths['t0'][0]))

    ## trise ( log truncnorm )
    log10_trise = get_truncated_normal(mean=p['log_trise'][0], sd=p['log_trise'][1], low=p['log_trise'][2], upp=p['log_trise'][3])
    ln_trise = np.log(log10_trise.pdf(np.log10(truths['trise'][0])))

    ## tfall ( log truncnorm )
    log10_tfall = get_truncated_normal(mean=p['log_tfall'][0], sd=p['log_tfall'][1], low=p['log_tfall'][2], upp=p['log_tfall'][3])
    ln_tfall = np.log(log10_tfall.pdf(np.log10(truths['tfall'][0])))

    ## extra_sigma ( log truncnorm )
    log10_extrasigma = get_truncated_normal(mean=p['log_extra_sigma'][0], sd=p['log_extra_sigma'][1], low=p['log_extra_sigma'][2], upp=p['log_extra_sigma'][3])
    ln_extrasigma = np.log(log10_extrasigma.pdf(p['log_extra_sigma'][0])) # TODO: check

    # posterior calc 
    ln_posterior = lnl_summed + ln_amp + ln_beta + ln_gamma + ln_t0 + ln_trise + ln_tfall + ln_extrasigma

    
    return ln_posterior

    

def superphot_sim_fits(fake_lc, chains, max_flux, n_idx = 30, jd0 = 2458119.5):

    ## posterior samples
    n_samples = len(chains.posterior.t0.values.flatten())
    rand_idx = np.random.choice(range(n_samples),
                                n_idx, replace=False)
    pi_t0 = chains.posterior.t0.values.flatten()[rand_idx]
    pi_amp =  max_flux * 10**chains.posterior.log_amplitude.values.flatten()[rand_idx]
    pi_beta = chains.posterior.beta.values.flatten()[rand_idx]
    pi_gamma = 10**chains.posterior.log_gamma.values.flatten()[rand_idx]
    pi_trise = 10**chains.posterior.log_trise.values.flatten()[rand_idx]
    pi_tfall = 10**chains.posterior.log_tfall.values.flatten()[rand_idx]

    ## overplot posterior draws
    t_grid = jnp.linspace(pi_t0[0] - 150,
                            np.max(fake_lc.jd.values) - jd0,
                            num=1000)
    fit_arrays = []
    for e in range(n_idx):
        f_model = y_model_superphot(t_grid, pi_amp[e], pi_beta[e], pi_t0[e], pi_gamma[e],
                                        pi_trise[e], pi_tfall[e], 0)
        fit_arrays.append(f_model)

    return t_grid, fit_arrays

def steps_from_npy(path, nchains, nsteps, param_names = ['log_amplitude', 'beta', 'log_gamma', 't0', 'log_trise', 'log_tfall', 'log_extra_sigma']):
    loaded = np.load(path, allow_pickle=True).item()
    if loaded['beta'].shape[0]/nsteps != nchains:
        print('error loading warmup: something is wrong with input nchains and/or nsteps')
        return
    else:
        chainsteps = []
        for param in param_names:
            chainsteps.append(loaded[param].reshape(nchains,nsteps))
    return chainsteps

def params_table_truth_init(warmup_path, nchains, nsteps_warmup, truth_path=None):
    if truth_path is not None:
        # true values, make base df
        truth = pd.read_csv(truth_path)
        tab_logs = pd.DataFrame(np.array([np.log10(truth['amp']), truth['beta'], np.log10(truth['gamma']), 
                                        truth['t0'], np.log10(truth['trise']), np.log10(truth['tfall']), np.array([0])]).T, 
                                        columns=['log_amp', 'beta', 'log_gamma', 't0', 'log_trise', 'log_tfall', 'log_extra_sigma'])
        index_names = ['truth']
    else:
        tab_logs = pd.DataFrame(columns=['log_amp', 'beta', 'log_gamma', 't0', 'log_trise', 'log_tfall', 'log_extra_sigma'])
        index_names = []
    pd.set_option('display.precision', 3)
    pd.set_option('display.width', 90)
    
    
    # add rows w inital values from warmup chains
    warmup_steps = steps_from_npy(warmup_path, nchains, nsteps_warmup)
    for c in range(nchains):
        tab_logs.loc[-1] = [float(array[c][0]) for array in warmup_steps]
        tab_logs.index = tab_logs.index + 1  # shifting index
        index_names.append(f'c{c}_init') # name this row
    tab_logs.index = index_names # assign row names
    return tab_logs

def summary_plot(sn_name, lc_folder_path, fit_folder_path, band, is_sim, init_method, nsteps_warmup, fmax_inject = None, xlim_for_draws=None, save_fig=True, jd0 = 2458119.5):

    #- run and fit a simulation
    #- as output w plots for a single run;
    #- “truth” values of simulation
    #- exactly where we are initializing (store warmup)
    #    - (clearly something strange, e.g. trise at 2*)
    #- plot of every individual chain for every parameter (see whether or not params are where we expect)
    #- plots of log likelihood
    #    - one i calculate and numpyro calculates
    #- plots of the posterior (log likelihood + log prior)
    #- plot of some fits

    # define all paths based on inputs
    fitname = f'{fit_folder_path}/{sn_name}_{band}_numpyro.nc'
    chains = az.from_netcdf(fitname) #### many things come from this
    # access parameter values across nchains and over nsteps
    nchains = chains.posterior.log_amplitude.values.shape[0]
    nsteps = chains.posterior.log_amplitude.values.shape[1]
    fake_lc_path = f'{lc_folder_path}/{sn_name}_fnu.csv'
    truth_path = None
    if is_sim:
        truth_path = f'{lc_folder_path}/{sn_name}_truth.csv'
        truth = pd.read_csv(truth_path)
        #fmax_inject = 328.73734846133704 #FMax_realsn
        if fmax_inject is None:
            fmax_inject = 328.73734846133704 #fmax of 'ZTF18abxoqkd', model SN for fakeSN{0-4}
            print('warning; no `fmax_inject` provided; what SN maxflux did you use to make these fake SNe?')
        fmax_eject = float(
            glob.glob(f'{fit_folder_path}/{sn_name}_{band}_maxflux*')[0].split('maxflux_')[1].split('.npy')[0])
    else:
        fmax = np.max(fake_lc['fnu_microJy']) # TODO: THIS WILL CAUSE ISSUES...
        print('ERRORRRR  i need fmax_inject AND fmax_eject')
    warmupsteppath = f'{fit_folder_path}/{sn_name}_{band}_warmupsamples.npy'
    fake_lc_path_smooth = f'{lc_folder_path}/{sn_name}_fnu_smooth.csv' # this will eventually interact w band
    # open dfs, chains using paths
    fake_lc_smooth = pd.read_csv(fake_lc_path_smooth)
    fake_lc = pd.read_csv(fake_lc_path) 

    ## markers and color palettes (needed for viz. of different chains)
    markers = [':', '-', ':', '-', ':', '-', ':', '-']
    chain_lws = [2, 0.9, 2, 0.9, 2, 0.9, 2, 0.9]
    chain_alphas = [0.9, 0.9, 0.9, 0.9 ,0.7, 0.7, 0.7, 0.7]
    chain_cols = ['#9CA3DB', '#9CA3DB', '#677DB7', '#677DB7', 'red', 'brown', 'cornflowerblue', 'black','goldenrod', 'hotpink', 'green']
    truth_col = "#B52E2C" 
    post_draws_col = '#191308'
    color_dict = {'g': "MediumAquaMarine", 'r': "Crimson", 'i': "GoldenRod"}


    ## skeleton of subplots
    # fig, axes = plt.subplot_mosaic("AABC;AADE;FFGH;IIJK;LLMN", figsize=(14,8)) # format without SNR plot
    fig, axes = plt.subplot_mosaic("AABC;AADE;SSGH;FXJK;ILMN", figsize=(13,9)) # format with SNR plot
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = 'Avenir'


    ## remove underlying axes to place text instead
    for ax_name in ["B", "C", "D", "X"]:
        axes[ax_name].remove()
    textstart_x = 0.52
    fig.text(textstart_x, 0.96, f'fit: {fitname}')
    fig.text(textstart_x, 0.77, f'init. method: {init_method}')
    if is_sim:
        ptable = params_table_truth_init(warmupsteppath, nchains, nsteps_warmup, truth_path)
        fig.text(textstart_x, 0.82, f'\n{ptable}')
        lnl_from_truth = lnl_from_sim(fake_lc, truth, max_flux = fmax_inject) #TODO: THIS IS CRUCIAL
        fig.text(0.3, 0.35, f'lnl from truth: {lnl_from_truth:.3f}')
        post_from_sim = posterior_from_sim_superphot(fake_lc, truth, max_flux=fmax_inject)[0]
        fig.text(0.3, 0.3, f'posterior from truth: {post_from_sim:.3f}')
    else:
        ptable = params_table_truth_init(warmupsteppath, nchains, nsteps_warmup)
        fig.text(textstart_x, 0.82, f'\n{ptable}')

    ## axis lims and labels
    axes["S"].set_xlabel('time (JD - 2018 Jan 01)')
    if xlim_for_draws is not None:
        axes["A"].set_xlim(xlim_for_draws[0], xlim_for_draws[1])
        axes["S"].set_xlim(xlim_for_draws[0], xlim_for_draws[1])
    for ax_name in ['L', 'M', 'N', 'I', 'L']:
        axes[ax_name].set_xlabel('step #')
    for ax_name in ['E', 'F', 'G', 'H', 'J', 'K', 'S', 'X']:
        axes[ax_name].set_xticks([])
        axes[ax_name].set_xticklabels([])
    axes["A"].set_ylabel(r'flux ($\mu$Jy)')
    axes["S"].set_ylabel('SNR')
    axes["F"].set_ylabel('lnl (numpyro)')
    axes["I"].set_ylabel('lnl (saarah)')
    axes["L"].set_ylabel('posterior (saarah)') 


    ## lnl and posterior plots
    lnl_numpyro = np.sum(chains.log_likelihood.y.values, axis=2)
    lnl_saarah = lnl_from_chains(fake_lc, chains, max_flux = fmax_eject)
    post_saarah = posterior_from_chains_superphot(fake_lc, chains, max_flux = fmax_eject)
    for i in range(nchains):
        axes["F"].plot(np.arange(len(lnl_numpyro[i])), lnl_numpyro[i],alpha=chain_alphas[i],lw=chain_lws[i],ls=markers[i], c=chain_cols[i])
        axes["I"].plot(np.arange(len(lnl_saarah[i])), lnl_saarah[i], alpha=chain_alphas[i],lw=chain_lws[i],ls=markers[i], c=chain_cols[i], label='chain '+str(i))
        axes["L"].plot(np.arange(len(post_saarah[i])), post_saarah[i], alpha=chain_alphas[i],lw=chain_lws[i],ls=markers[i], c=chain_cols[i])

    ## parameter step plots
    param_names = ["log_amp", "beta", "log_gamma", "t0", "log_trise", "log_tfall", "log_extra_sigma"]
    param_vals = [chains.posterior.log_amplitude, chains.posterior.beta, chains.posterior.log_gamma, chains.posterior.t0, chains.posterior.log_trise, chains.posterior.log_tfall, chains.posterior.log_extra_sigma]
    for n, ax_name in enumerate(["G", "H", "J", "K", "M", "N", "E"]):
        axes[ax_name].set_ylabel(param_names[n])
        for i in range(nchains):
            axes[ax_name].plot(np.arange(nsteps), param_vals[n][i].values, lw=chain_lws[i], ls=markers[i], c=chain_cols[i], alpha=chain_alphas[i])
        if is_sim:
            if param_names[n] == 'log_extra_sigma':
                continue
            elif param_names[n] not in truth.columns:
                tru_val = np.log10(truth[param_names[n].split("_")[1]].values[0])
                axes[ax_name].axhline(tru_val, c=truth_col, lw=5)
            else:
                tru_val = truth[param_names[n]].values[0]
                if param_names[n] == 'beta':
                    axes[ax_name].axhline(tru_val, c=truth_col, lw=5, label='true value')
                else:
                    axes[ax_name].axhline(tru_val, c=truth_col, lw=5)



    ## one legend for all chain-based plots
    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    fig.legend(lines, labels, loc=(0.901,0.815))


    ## posterior draws and data
    if is_sim:
        axes["A"].plot(fake_lc_smooth['jd']-jd0, fake_lc_smooth['fnu_microJy'], c=truth_col, ls='-', lw=5, label='simulated model')
        ## posterior draws and data
        axes["A"].errorbar(fake_lc['jd']- jd0, fake_lc['fnu_microJy'], yerr=fake_lc['fnu_microJy_unc'], fmt='o', color=color_dict[band], 
                    label='noisy simulated data')
    else:
        axes["A"].errorbar(fake_lc['jd']- jd0, fake_lc['fnu_microJy'], yerr=fake_lc['fnu_microJy_unc'], fmt='o', color=color_dict[band], 
                    label='data')
    fits_to_plot = superphot_sim_fits(fake_lc, chains, fmax_eject)
    for p in range(len(fits_to_plot[1])):
        if p == 0:
            axes["A"].plot(fits_to_plot[0], fits_to_plot[1][p], color=post_draws_col, ls='--', lw=1.2, alpha=0.3, label='fit model')
        else:
            axes["A"].plot(fits_to_plot[0], fits_to_plot[1][p], color=post_draws_col, ls='--', lw=1.2, alpha=0.3)
    axes["A"].legend()

    ## SNR vs. time plot (baseline investigation)
    SNR = (fake_lc.fnu_microJy - fake_lc.fnu_microJy_unc) / (np.sqrt(fake_lc.fnu_microJy_unc))
    axes["S"].plot(fake_lc.jd - jd0, SNR, 's', ms=1, color='k')
    snrplot_yspan = axes["S"].get_ylim()
    axes["S"].axhspan(snrplot_yspan[1], 5, alpha=0.4, label = 'SNR > 5', color='darkseagreen')
    axes["S"].axhspan(0, snrplot_yspan[0], alpha=0.4, label = 'SNR < 0', color='slategray')
    axes["S"].legend()


    plt.tight_layout()
    plt.subplots_adjust(hspace=0.20)
    plt.subplots_adjust(wspace=0.30)

    if save_fig:
            fig.savefig(f"{fit_folder_path}/{sn_name}_fitsummary.png",
                        dpi=900, transparent=False)

    return


# --- ---
def mkdir_and_numpyro(sn, out_path, nwarmup, nsample, nchain, init_strat, model, randomkey, lc_path='bts_lcs/'):
    bool_out = os.path.isdir(out_path)
    if not bool_out:
        print('time to make a new path')
        os.mkdir(out_path)
        print('made a new path')
    
    fit_gr_numpyro(sn, lc_path, out_path, nwarmup, nsample, nchain, 
                   init_strat, model, randomkey=randomkey)
      
def main():

    # Initialize argument parser
    parser = argparse.ArgumentParser(prog='fit_villar_numpyro.py <sn>',
                                     description='Run Villar+19 light curve fits on a BTS fps lc. with numpyro')
    # Necessary arguments
    parser.add_argument('sn', type=str, nargs='?', default=None,
                   help='ZTF transient name')
    # Optional arguments
    parser.add_argument('lc_path', type=str, nargs='?', default=None,
                   help='path to folder containing processed fps lc file from Miller+24')
    parser.add_argument('out_path', type=str, nargs='?', default=None,
                   help='path for output MCMC chains')
    parser.add_argument('num_warmup', type=int, nargs='?', default=None,
                        help='number of warmup steps for MCMC')
    parser.add_argument('num_samples', type=int, nargs='?', default=None,
                        help='number of samples for MCMC')
    parser.add_argument('num_chains', type=int, nargs='?', default=None,
                        help='number of chains for MCMC')
    parser.add_argument('init_strat', type=str, nargs='?', default=None,
                        help='initialization strategy for MCMC')

    try:
        args = parser.parse_args()

        run = True
    except Exception:
        run = False

    if run:
        fit_gr_numpyro(**vars(args))

if __name__ == "__main__":
    main()