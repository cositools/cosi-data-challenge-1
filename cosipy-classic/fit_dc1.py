import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import matplotlib as mpl

from tqdm.autonotebook import tqdm
from IPython.display import HTML

import warnings
warnings.filterwarnings('ignore')

import pandas as pd
from shapely.geometry import Polygon

import pickle

# spectral fits
from priors_dc1 import *
import time
import emcee
import corner
import scipy.optimize as op



class fit():
    """
    Fitting class that includes the dataset to be analysed, the pointings, the response, and a background model.
    :option: bg_only:   Default = False: performs a background only fit when True
    :option: priors:    Default: uninformative (normal) priors for the sky;
    :option: verbose:   Default = False: more verbose output
    """

    def __init__(self,
                 dataset,          # COSIpy dataset
                 pointings,        # POINTINGS object
                 response,         # SkyResponse object for certain source position
                 background,       # BG object including cuts and tracer
                 bg_only=False,    # True performs BG only fit
                 priors=None,      # can set priors for sky components
                 verbose=False,    # more verbose output
                 reduced=True):    # use full data space
    
        # init objects
        self.dataset = dataset
        self.pointings = pointings
        self.response = response
        self.background = background

        self.count_data = self.reduce_dataset_CDS()  # construct data set with reduced CDS (ignore bins that are always zero)
        
        self.bg_only = bg_only
        self.priors = priors
        self.verbose = verbose

        self.data_per_energy_bin = self.make_dictionaries_for_stan() # make dictionaries for fit with Stan
   
    
    def reduce_dataset_CDS(self):
        
        count_data_reduced = [] # init list of data sets per bin, will be irregularly-shaped because CDS population depends on energy

        for i in range(self.dataset.energies.n_energy_bins):
            
            yp_tmp = self.dataset.binned_data[:,i,:,:].reshape(self.dataset.times.n_ph,self.background.bg_model.shape[2]*self.background.bg_model.shape[3])  # reshape count data to reduce CDS if possible, this combines the 3 CDS angles into a 1D array for all times at the chosen energy
    
            yp_tmp = yp_tmp[:,self.background.calc_this[i]] # reshape count data grid the same way as background and choose only non-zero indices

            count_data_reduced.append(yp_tmp) # append reduced data set
            
        return count_data_reduced

    
    
    def make_dictionaries_for_stan(self):
        """
        Create dictionaries that can be read in by the Stan model to fit the data
        """

        all_dicts = [] # init dictionaries for each energy bin

        # loop over energy bins
        for i in range(self.dataset.energies.n_energy_bins):
            Np, Nrsp = self.background.bg_model_reduced[i].shape         # initialise sizes of arrays
            N = Np*Nrsp                                                  # total number of data points

            Nsky = 1  # right now, only one sky model allowed to be fitted

            if np.any(self.priors == None): # standard priors scaling the initial flux as in response calculation

                mu_flux_scl = np.array([0.])    # mean, prior centroids for sky, we use 10 because we are ignorant; this has to be an array because it could be more than one
                sigma_flux_scl = np.array([1e8]) # std, same for the width (so, easily 0 but also high values possible)

                # priors for background model components
                mu_Abg = 1.       # for the moment set to a useful value if bg model is ~normalised to data, initially normalised to 1, so mean would be 1, variance very large (uninformative)
                sigma_Abg = 1e4   # same
                
            else:

                # set priors yourself
                mu_flux_scl = np.array([self.priors[0,0]])
                sigma_flux_scl = np.array([self.priors[0,1]])

                mu_Abg = self.priors[1,0]       # for the moment set to a useful value if bg model is ~normalised to data
                sigma_Abg = self.priors[1,1]   # same

            # dictionary for data set and prior
            data2D = dict(N = Nrsp,                                                      # number of CDS bins
                          Np = Np,                                                       # number of observations
                          Nsky = Nsky,                                                   # number of sky models (now: 1)
                          Ncuts = self.background.Ncuts,                                 # number of background cuts / renormalisations
                          bg_cuts = self.background.bg_cuts,                             # bg cuts at
                          bg_idx_arr = self.background.idx_arr,                          # bg cut indices
                          y = self.count_data[i].ravel().astype(int),                    # data
                          bg_model = self.background.bg_model_reduced[i],                # background model 
                          conv_sky = self.response.sky_response[i].reshape(Nsky,Np,Nrsp),# this has to be reshaped because it could be more than one
                          mu_flux = mu_flux_scl,                                         # priors for sky (mean)
                          sigma_flux = sigma_flux_scl,                                   # std
                          mu_Abg = mu_Abg,                                               # BG mean
                          sigma_Abg = sigma_Abg)                                         # std

            all_dicts.append(data2D) # append dictionary for energy bin
            
        return all_dicts
          
            
    def fit(self,iters=1000,pars=['flux','Abg']):
        """
        Fitting COSIpy fit object of a data set with pointing definition, background model, and (now only) point source response.
        Fitting background only is only possible when object is initialised with bg_only = True.
        :option: par    Parameters to save in fit object: pars=['flux','Abg','model_tot','model_bg','model_sky','ppc'], default pars=['flux','Abg'],
                        i.e. no models will be saved, only the fitted parameters.
        :option: iters  Number of iterations for fit, default 1000.
        Saves fitting results in .fit_pars, including all posterior distributions.
        Creates .diff_rate and .diff_rate_err (1sigma uncertainty) that includes the differential count rate in units of ph/s/keV for all energy bins
        """

        # init arrays
        self.fit_pars = []
        self.diff_rate = np.zeros(self.dataset.energies.n_energy_bins)
        self.diff_rate_err = np.zeros((self.dataset.energies.n_energy_bins,2))

        self.systematics = np.array([9.61348192, 4.5787978 , 1.15368693, 1.22003529, 1, 1, 1, 1, 1, 1]) # standard systematics (from sky-only fits)
        
        for i in range(self.dataset.energies.n_energy_bins):
            print("working on energy bin " + str(i))
            init = np.array([np.sum(self.data_per_energy_bin[i]['y'])*0.05,0.99]) # initial guess for emcee fitting (1 bg parameter only)
            init_var = init*1e-4

            # here standard emcee workflow
            ndim, nwalkers = 2, 10
            pos = [init + np.random.randn(ndim)*init_var for i in range(nwalkers)]
                    
            sampler = emcee.EnsembleSampler(nwalkers,
                                            ndim,
                                            COSImodfit,
                                            args = (self.data_per_energy_bin[i]['y'], #data
                                                    self.data_per_energy_bin[i]['conv_sky'].ravel(), # sky model
                                                    self.data_per_energy_bin[i]['bg_model'].ravel())) # bg model)
                    
            _ = sampler.run_mcmc(pos, iters, progress=True)

            # extract samples
            samples = sampler.get_chain()
            samplesf = sampler.flatchain
            n_par = 2
                    
            n_samples = iters
            n_walkers = nwalkers
                    
            burnin = int(0.5*n_samples)
                        
            spec_params = np.zeros((n_par,7))

            for p in range(n_par):
                mean_val   = np.mean(samples[burnin:,:,p])
                std_val    = np.std(samples[burnin:,:,p])
                median_val = np.median(samples[burnin:,:,p])
                ub1_val    = np.percentile(samples[burnin:,:,p],50+68.3/2)
                lb1_val    = np.percentile(samples[burnin:,:,p],50-68.3/2)
                ub3_val    = np.percentile(samples[burnin:,:,p],50+99.73/2)
                lb3_val    = np.percentile(samples[burnin:,:,p],50-99.73/2)
                spec_params[p,:] = [mean_val,std_val,lb3_val,lb1_val,median_val,ub1_val,ub3_val]

            self.fit_pars.append(spec_params) # append fitting results

            norm_tmp = self.dataset.energies.energy_bin_wid[i]*2/self.response.flux_norm # calculate rates

            self.diff_rate[i] = self.fit_pars[i][0,4]/norm_tmp # median value is representative for spectrum

            # 1sigma error bars
            self.diff_rate_err[i,1] = self.fit_pars[i][0,5]/norm_tmp - self.diff_rate[i] # upper boundary uncertainty
            self.diff_rate_err[i,0] = self.fit_pars[i][0,3]/norm_tmp - self.diff_rate[i] # lower boundary uncertainty


    def plot_extracted_spectrum(self,save_file,ul=None,with_systematics=False,col1='black',col2='gray'):

        if with_systematics == True:
            syst_scl = self.systematics
        else:
            syst_scl = np.ones(len(self.dataset.energies.energy_bin_cen))
        
        plt.figure(figsize=(10.24,7.68))

        snr = self.diff_rate/np.nanmax(np.abs(self.diff_rate_err.T)*syst_scl,axis=0)
        if ul == None:
            plt.errorbar(self.dataset.energies.energy_bin_cen,
                         self.diff_rate,
                         xerr=self.dataset.energies.energy_bin_wid,
                         yerr=np.abs(self.diff_rate_err.T)*syst_scl[None,:],
                         fmt='o',linewidth=3,color=col1,markersize=10,
                         label=r'COSI Data Fit ($1\sigma$)')
        else:
            # mark data points which have SNR < ul with downward arrows
            # snr = self.diff_rate/np.nanmax(np.abs(self.diff_rate_err.T)*syst_scl,axis=0)
            ul_idx = np.where(snr < ul)[0]
            gl_idx = np.where(snr >= ul)[0]
            plt.errorbar(self.dataset.energies.energy_bin_cen[gl_idx],
                           self.diff_rate[gl_idx],
                           xerr=self.dataset.energies.energy_bin_wid[gl_idx],
                           yerr=np.abs(self.diff_rate_err.T)[:,gl_idx]*syst_scl[None,gl_idx],
                           fmt='o',linewidth=3,color=col1,markersize=10,
                           label=r'COSI Data Fit ($1\sigma$)')
            plt.errorbar(self.dataset.energies.energy_bin_cen[ul_idx],
                           np.max(np.abs(self.diff_rate_err.T),axis=0)[ul_idx]*ul*syst_scl[ul_idx],
                           xerr=self.dataset.energies.energy_bin_wid[ul_idx],
                           yerr=np.zeros(len(ul_idx)),
                           fmt='v',linewidth=3,color=col2,markersize=15,
                           label=r'COSI Data Fit ({0}$\sigma$ upper limit)'.format(ul))

        plt.xlabel('Energy [keV]')
        plt.ylabel('Counts [cnts/keV]')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlim(150,5000)
        plt.legend()

        # Write data:
        save_path = os.getcwd()
        save_file = os.path.join(save_path,save_file)
        d = {"ebin_center[keV]":self.dataset.energies.energy_bin_cen,\
                        "xerr[keV]":self.dataset.energies.energy_bin_wid,\
                        "rate[ct/keV]":self.diff_rate,\
                        "rate_err_low[ct/keV]":np.abs(self.diff_rate_err.T)[0],\
                        "rate_err_high[ct/keV]":np.abs(self.diff_rate_err.T)[1],\
                        "SNR":snr}
        df = pd.DataFrame(data=d)
        df.to_csv(save_file,float_format='%10.5e',index=False,sep="\t",\
                columns=["ebin_center[keV]","xerr[keV]","rate[ct/keV]",\
                "rate_err_low[ct/keV]","rate_err_high[ct/keV]","SNR"])

def COSImodfit(theta, data, sky_model, background_model, eval=False):
    """
    Returns:
    cash-stat for model fit
    
    one-parameter background model used only

    Parameters:
    :param theta:     Array of to-be-fitted parameters 
    :param other: ...
    """

    total_model = theta[0]*sky_model + theta[1]*background_model
    
    stat = -2*np.sum(total_model - data*(np.log(total_model)))
    if np.isnan(stat) | (theta[0] < 0) | (theta[1] < 0):
        return -np.inf
    else:
        return stat


