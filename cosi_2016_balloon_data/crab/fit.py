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
#from COSIpy import FISBEL
#from COSIpy import dataset
#from COSIpy import GreatCircle
#from COSIpy import angular_distance

import pickle
import pystan

# spectral fits
from priors import *
from spectral_shapes import *
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

        # construct data set with reduced CDS (ignore bins that are always zero)
        self.count_data = self.reduce_dataset_CDS()
        
        self.bg_only = bg_only
        self.priors = priors
        self.verbose = verbose

        # make dictionaries for fit with Stan
        # load corresponding Stan model
        if not self.bg_only:
            self.data_per_energy_bin = self.make_dictionaries_for_stan()
            self.load_stan_model()
        else:
            if self.verbose:
                print('This will be a background only fit.')
            self.data_per_energy_bin = self.make_dictionaries_for_stan_bg_only()
            self.load_stan_model_bg_only()
   
    
    
    def reduce_dataset_CDS(self):
        
        # init list of data sets per bin
        # will be irregularly-shaped because CDS population depends on energy
        count_data_reduced = []

        # loop ove renergies
        for i in range(self.dataset.energies.n_energy_bins):
            
            # reshape count data to reduce CDS if possible
            # this combines the 3 CDS angles into a 1D array for all times at the chosen energy
            yp_tmp = self.dataset.binned_data[:,i,:,:].reshape(self.dataset.times.n_ph,#self.dataset.times.n_time_bins,
                                                          self.background.bg_model.shape[2]*self.background.bg_model.shape[3])
    
            # reshape count data grid the same way as backgrorund and choose only non-zero indices
            yp_tmp = yp_tmp[:,self.background.calc_this[i]]

            # append reduced data set
            count_data_reduced.append(yp_tmp)
            
        return count_data_reduced

    
    
    def make_dictionaries_for_stan(self):
        """
        Create dictionaries that can be read in by the Stan model to fit the data
        """

        # init dictionaries for each energy bini
        all_dicts = []

        # loop over energy bins
        for i in range(self.dataset.energies.n_energy_bins):
            Np, Nrsp = self.background.bg_model_reduced[i].shape         # initialise sizes of arrays
            N = Np*Nrsp                                                  # total number of data points

            # right now, only one sky model allowed to be fitted
            Nsky = 1

            # standard priors scaling the initial flux as in response calculation
            if np.any(self.priors == None):

                # mean
                mu_flux_scl = np.array([0.])    # prior centroids for sky, we use 10 because we are ignorant;
                # this has to be an array because it could be more than one
                # std
                sigma_flux_scl = np.array([1e8]) # same for the width (so, easily 0 but also high values possible)

                # priors for backgrorund model components
                # initially normalised to 1, so mean would be 1, variance very large (uninformative)
                mu_Abg = 1.       # for the moment set to a useful value if bg model is ~normalised to data
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

            # append dictionary for energy bin
            all_dicts.append(data2D)
            
        return all_dicts
    
    
    def make_dictionaries_for_stan_bg_only(self):
        """
        Same as above just for background-only fit
        """
        # init dictionary per energy bin
        all_dicts = []

        # loop over energies
        for i in range(self.dataset.energies.n_energy_bins):
            
            Np, Nrsp = self.background.bg_model_reduced[i].shape         # initialise sizes of arrays
            N = Np*Nrsp                                                  # total number of data points

            # standard priors scaling the initial flux as in response calculation
            if np.any(self.priors == None):

                mu_Abg = 1.       # for the moment set to a useful value if bg model is ~normalised to data
                sigma_Abg = 1e4   # same

            else:

                # set priors yourself
                mu_Abg = self.priors[1,0]       # for the moment set to a useful value if bg model is ~normalised to data
                sigma_Abg = self.priors[1,1]   # same

            
            # dictionary for data set and prior
            data2D = dict(N = Nrsp,
                          Np = Np,
                          Ncuts = self.background.Ncuts,
                          bg_cuts = self.background.bg_cuts,
                          bg_idx_arr = self.background.idx_arr,
                          y = self.count_data[i].ravel().astype(int),
                          bg_model = self.background.bg_model_reduced[i],
                          mu_Abg = mu_Abg,
                          sigma_Abg = sigma_Abg)
            all_dicts.append(data2D)
            
        return all_dicts
    
    
    def load_stan_model(self):
        """
        Loading the Stan model COSImodfit.stan.
        Compiles it if not already done.
        """
        try:
            #read COSImodefit.pkl (if already compiled)
            self.model = pickle.load(open('COSImodfit.pkl', 'rb'))
            
        except:
            print('Model not yet compiled, doing that now (might take a while).')
            ## compile model (if not yet compiled):
            self.model = pystan.StanModel('COSImodfit.stan')

            ## save it to the file 'filename.pkl' for later use
            with open('COSImodfit.pkl', 'wb') as f:
                pickle.dump(self.model, f)
    
    
    def load_stan_model_bg_only(self):
        """
        Loading Stan model for background only.
        Compiles it of not already done.
        """
        try:
            #read COSImodfit_BGonly.pkl (if already compiled)
            self.model = pickle.load(open('COSImodfit_BGonly.pkl', 'rb'))
            
        except:
            print('Model not yet compiled, doing that now (might take a while).')
            ## compile model (if not yet compiled):
            self.model = pystan.StanModel('COSImodfit_BGonly.stan')

            ## save it to the file 'filename.pkl' for later use
            with open('COSImodfit_BGonly.pkl', 'wb') as f:
                pickle.dump(self.model, f)
                
    
    def MAP_solution(self,guess=1.0,method='LBFGS',scipy=False,ebins='all'):
        """
        Performs optimisation of the joint posterior distribution and returns Median A-Posteriori (MAP) point estimate.
        Returns no error bars and serves as quick cross check or for calls to likelihood ratio tests.
        Creates array .diff_rate_map that includes the differential count rate in units of ph/s/keV for each energy bin.
        Saves all fit results and quality in .fit_pars_map.
        """
        # init arrays to save info
        self.fit_pars_map = []
        self.diff_rate_map = np.zeros(self.dataset.energies.n_energy_bins)

        if scipy == True:
            if not self.bg_only:
                nll = lambda *args: -COSImodfit(*args)
            else:
                nll = lambda *args: -COSImodfit_bgonly(*args)

        if np.any(ebins == 'all'):
            ebins = np.arange(self.dataset.energies.n_energy_bins)

        guess_arr = np.array([130351.9, 1018422.0, 2676585.3, 3723391.2, 4756865.1, 5539038.8, 5843760.6, 5051422.9, 3110379.8, 2139293.1])
        
        # loop over energy bins
        for i in tqdm(range(self.dataset.energies.n_energy_bins),desc='Loop over energy bins:'):

            if ebins[i] == i:
                #print('energy bin', i, ':', np.sum(self.data_per_energy_bin[i]['y']), np.sum(self.data_per_energy_bin[i]['bg_model'].ravel()), np.sum(self.data_per_energy_bin[i]['conv_sky'].ravel()))
                
                if self.verbose:
                    print('Start optimising energy bin '+str(i+1)+'/'+str(self.dataset.energies.n_energy_bins)+'...')
                    print('\nEnergy range: '+str(self.dataset.energies.energy_bin_min[i])+'-'+str(self.dataset.energies.energy_bin_max[i])+' keV ...')
                
                # not necessary in general, but good for quicklook MAP estimate
                init = {}
                if not self.bg_only:
                    init['flux'] = np.array([np.sum(self.data_per_energy_bin[i]['y'])*0.05])
                    init['Abg'] = np.repeat(1.01,self.background.Ncuts)
                else:
                    init['Abg'] = np.repeat(1.01,self.background.Ncuts)


                if scipy == False:
                    # optimising the model
                    res = self.model.optimizing(data=self.data_per_energy_bin[i],
                                                verbose=False,init=init,as_vector=False,algorithm=method)#,tol_rel_grad=1e4)
                else:
                    if not self.bg_only:
                        res = op.minimize(nll,
                                          [np.sum(self.data_per_energy_bin[i]['y'])*0.015,1.],
                                          args=(self.data_per_energy_bin[i]['y'], #data
                                                self.data_per_energy_bin[i]['conv_sky'].ravel(), # sky model
                                                self.data_per_energy_bin[i]['bg_model'].ravel()), # bg model
                                          options={'gtol': 1e-4})
                        print(res.x, res.x[0]*np.sum(self.data_per_energy_bin[i]['conv_sky'].ravel()), res.x[1]*np.sum(self.data_per_energy_bin[i]['bg_model'].ravel()))
                    else:
                        res = op.minimize(nll,
                                          [np.sum(self.data_per_energy_bin[i]['y'])*0.015,1.],
                                          args=(self.data_per_energy_bin[i]['y'], #data
                                                self.data_per_energy_bin[i]['bg_model'].ravel()), # bg model
                                          options={'gtol': 1e-4})
                    

                # append the result
                self.fit_pars_map.append(res)

                # calculate the flux
                if not self.bg_only:
                    norm_tmp = self.dataset.energies.energy_bin_wid[i]*2/self.response.flux_norm*self.dataset.times.total_time
                    if scipy == False:
                        #self.diff_rate_map[i] = self.fit_pars_map[i]['par']['flux']/norm_tmp
                        self.diff_rate_map[i] = res['par']['flux']/norm_tmp
                    else:
                        #self.diff_rate_map[i] = self.fit_pars_map[i].x[0]/norm_tmp
                        self.diff_rate_map[i] = res.x[0]/norm_tmp
                        print(res.x[0]/norm_tmp)
                        
            else:

                if self.verbose:
                    print('Skipping energy range: '+str(self.dataset.energies.energy_bin_min[i])+'-'+str(self.dataset.energies.energy_bin_max[i])+' keV ...')
                # append a zero result
                #self.fit_pars_map.append(0)
            
            
    def fit(self,iters=1000,pars=['flux','Abg'],use_emcee=False):
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
        self.diff_rate_err2 = np.zeros((self.dataset.energies.n_energy_bins,2))


        # standard systematics (from sky-only fits)
        self.systematics = np.array([9.61348192, 4.5787978 , 1.15368693, 1.22003529, 1, 1, 1, 1, 1, 1])
        #0.75795775, 0.44679336, 0.7485999 , 0.82748034, 0.7139853 , 0.16266929])
        
        # loop over energy bins
        for i in tqdm(range(self.dataset.energies.n_energy_bins),desc='Loop over energy bins:'):

            if self.verbose:
                print('###################################################################')
                print('\nStart fitting energy bin '+str(i+1)+'/'+str(self.dataset.energies.n_energy_bins)+'...')
                print('\nEnergy range: '+str(self.dataset.energies.energy_bin_min[i])+'-'+str(self.dataset.energies.energy_bin_max[i])+' keV ...')
                

            # fit including sky
            if not self.bg_only:

                if not use_emcee:
                
                    # sample full posterior
                    fit = self.model.sampling(data=self.data_per_energy_bin[i],
                                              chains=1,iter=iters,n_jobs=-1,verbose=False,
                                              pars=pars,control={'adapt_delta':0.9})

                    if self.verbose:
                        print('Summary for energy bin '+str(i+1)+'/'+str(self.dataset.energies.n_energy_bins)+':\n')
                        print(fit.stansummary(['flux','Abg']))

                    # append fitting results
                    self.fit_pars.append(fit)

                    # calculate rates
                    norm_tmp = self.dataset.energies.energy_bin_wid[i]*2/self.response.flux_norm*self.dataset.times.total_time
                
                    # median value es representative for spectrum
                    self.diff_rate[i] = np.percentile(self.fit_pars[i]['flux'],50)/norm_tmp

                    # 1sigma error bars
                    # upper boundary uncertainty
                    self.diff_rate_err[i,1] = np.percentile(self.fit_pars[i]['flux'],50+68.3/2)/norm_tmp - self.diff_rate[i]
                    # lower boundary uncertainty
                    self.diff_rate_err[i,0] = np.percentile(self.fit_pars[i]['flux'],50-68.3/2)/norm_tmp - self.diff_rate[i]

                    # 2sigma error bars
                    # upper boundary uncertainty
                    self.diff_rate_err2[i,1] = np.percentile(self.fit_pars[i]['flux'],50+95.4/2)/norm_tmp - self.diff_rate[i]
                    # lower boundary uncertainty
                    self.diff_rate_err2[i,0] = np.percentile(self.fit_pars[i]['flux'],50-95.4/2)/norm_tmp - self.diff_rate[i]

                    
                else:

                    # emcee fitting (1 bg parameter only)
                    # initial guess
                    init = np.array([np.sum(self.data_per_energy_bin[i]['y'])*0.05,
                                     0.99])

                    #print('init',init)
                    
                    init_var = init*1e-4

                    # here standard emcee workflow
                    ndim, nwalkers = 2, 10
                    pos = [init + np.random.randn(ndim)*init_var for i in range(nwalkers)]

                    #print('pos',pos)
                    
                    start = time.time()
                    
                    sampler = emcee.EnsembleSampler(nwalkers,
                                                    ndim,
                                                    COSImodfit,
                                                    args = (self.data_per_energy_bin[i]['y'], #data
                                                            self.data_per_energy_bin[i]['conv_sky'].ravel(), # sky model
                                                            self.data_per_energy_bin[i]['bg_model'].ravel())) # bg model)
                    
                    _ = sampler.run_mcmc(pos, iters, progress=True)

                    end = time.time()

                    # extract samples
                    samples = sampler.get_chain()
                    samplesf = sampler.flatchain
                    n_par = 2
                    
                    n_samples = iters#*nwalkers
                    n_walkers = nwalkers
                    
                    burnin = int(0.5*n_samples)
                    
                    ttime = end - start
                    print("Processing took {0:.1f} seconds".format(ttime))

                    if self.verbose:

                        # output here
                        print('\n')
                        print('Results:\n')

                        
                    spec_params = np.zeros((n_par,7))

                    if self.verbose:
                    
                        # formatting the table
                        row_format ='{:>10}' * 8

                        # first table row
                        print(row_format.format(*['Parameter','mean','std','0.15','15.85','50.00','84.15','99.85']))

                    for p in range(n_par):
                        mean_val   = np.mean(samples[burnin:,:,p])
                        std_val    = np.std(samples[burnin:,:,p])
                        median_val = np.median(samples[burnin:,:,p])
                        ub1_val    = np.percentile(samples[burnin:,:,p],50+68.3/2)
                        lb1_val    = np.percentile(samples[burnin:,:,p],50-68.3/2)
                        ub3_val    = np.percentile(samples[burnin:,:,p],50+99.73/2)
                        lb3_val    = np.percentile(samples[burnin:,:,p],50-99.73/2)
                        spec_params[p,:] = [mean_val,std_val,lb3_val,lb1_val,median_val,ub1_val,ub3_val]

                        
                        if self.verbose:
                        
                            print(row_format.format(str(p)+':',
                                                    str('%1.2e' % mean_val),
                                                    str('%1.2e' % std_val),
                                                    str('%1.2e' % lb3_val),
                                                    str('%1.2e' % lb1_val),
                                                    str('%1.2e' % median_val),
                                                    str('%1.2e' % ub1_val),
                                                    str('%1.2e' % ub3_val)))

                    # append fitting results
                    self.fit_pars.append(spec_params)

                    # calculate rates
                    #norm_tmp = self.dataset.energies.energy_bin_wid[i]*2/self.response.flux_norm*self.dataset.times.total_time
                    norm_tmp = self.dataset.energies.energy_bin_wid[i]*2/self.response.flux_norm

                    # median value es representative for spectrum
                    self.diff_rate[i] = self.fit_pars[i][0,4]/norm_tmp

                    # 1sigma error bars
                    # upper boundary uncertainty
                    self.diff_rate_err[i,1] = self.fit_pars[i][0,5]/norm_tmp - self.diff_rate[i]
                    # lower boundary uncertainty
                    self.diff_rate_err[i,0] = self.fit_pars[i][0,3]/norm_tmp - self.diff_rate[i]


                    
            # BG-only fit
            else:

                # sample full posterior
                fit = self.model.sampling(data=self.data_per_energy_bin[i],
                                          chains=1,iter=iters,n_jobs=-1,verbose=False,
                                          pars=['Abg','model_tot','model_bg'])

                if self.verbose:
                    print('Summary for energy bin '+str(i+1)+'/'+str(self.dataset.energies.n_energy_bins)+':\n')
                    print(fit.stansummary(['Abg']))
                    
                # append fitting results
                self.fit_pars.append(fit)

            if self.verbose:
                print('###################################################################')

                
    def fit_new(self,iters=1000,pars=['flux','Abg'],use_emcee=False,time_div='tot'):
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
        self.diff_rate_err2 = np.zeros((self.dataset.energies.n_energy_bins,2))


        # standard systematics (from sky-only fits)
        self.systematics = np.array([9.61348192, 4.5787978 , 1.15368693, 1.22003529, 1, 1, 1, 1, 1, 1])
        #0.75795775, 0.44679336, 0.7485999 , 0.82748034, 0.7139853 , 0.16266929])
        
        # loop over energy bins
        for i in tqdm(range(self.dataset.energies.n_energy_bins),desc='Loop over energy bins:'):

            if self.verbose:
                print('###################################################################')
                print('\nStart fitting energy bin '+str(i+1)+'/'+str(self.dataset.energies.n_energy_bins)+'...')
                print('\nEnergy range: '+str(self.dataset.energies.energy_bin_min[i])+'-'+str(self.dataset.energies.energy_bin_max[i])+' keV ...')
                

            # fit including sky
            if not self.bg_only:

                if not use_emcee:
                
                    # sample full posterior
                    fit = self.model.sampling(data=self.data_per_energy_bin[i],
                                              chains=1,iter=iters,n_jobs=-1,verbose=False,
                                              pars=pars,control={'adapt_delta':0.9})

                    if self.verbose:
                        print('Summary for energy bin '+str(i+1)+'/'+str(self.dataset.energies.n_energy_bins)+':\n')
                        print(fit.stansummary(['flux','Abg']))

                    # append fitting results
                    self.fit_pars.append(fit)

                    # calculate rates
                    norm_tmp = self.dataset.energies.energy_bin_wid[i]*2/self.response.flux_norm*self.dataset.times.total_time
                
                    # median value es representative for spectrum
                    self.diff_rate[i] = np.percentile(self.fit_pars[i]['flux'],50)/norm_tmp

                    # 1sigma error bars
                    # upper boundary uncertainty
                    self.diff_rate_err[i,1] = np.percentile(self.fit_pars[i]['flux'],50+68.3/2)/norm_tmp - self.diff_rate[i]
                    # lower boundary uncertainty
                    self.diff_rate_err[i,0] = np.percentile(self.fit_pars[i]['flux'],50-68.3/2)/norm_tmp - self.diff_rate[i]

                    # 2sigma error bars
                    # upper boundary uncertainty
                    self.diff_rate_err2[i,1] = np.percentile(self.fit_pars[i]['flux'],50+95.4/2)/norm_tmp - self.diff_rate[i]
                    # lower boundary uncertainty
                    self.diff_rate_err2[i,0] = np.percentile(self.fit_pars[i]['flux'],50-95.4/2)/norm_tmp - self.diff_rate[i]

                    
                else:

                    # emcee fitting (1 bg parameter only)
                    # initial guess
                    init = np.array([np.sum(self.data_per_energy_bin[i]['y'])*0.05,
                                     0.99])

                    #print('init',init)
                    
                    init_var = init*1e-4

                    # here standard emcee workflow
                    ndim, nwalkers = 2, 10
                    pos = [init + np.random.randn(ndim)*init_var for i in range(nwalkers)]

                    #print('pos',pos)
                    
                    start = time.time()
                    
                    sampler = emcee.EnsembleSampler(nwalkers,
                                                    ndim,
                                                    COSImodfit,
                                                    args = (self.data_per_energy_bin[i]['y'], #data
                                                            self.data_per_energy_bin[i]['conv_sky'].ravel(), # sky model
                                                            self.data_per_energy_bin[i]['bg_model'].ravel())) # bg model)
                    
                    _ = sampler.run_mcmc(pos, iters, progress=True)

                    end = time.time()

                    # extract samples
                    samples = sampler.get_chain()
                    samplesf = sampler.flatchain
                    n_par = 2
                    
                    n_samples = iters#*nwalkers
                    n_walkers = nwalkers
                    
                    burnin = int(0.5*n_samples)
                    
                    ttime = end - start
                    print("Processing took {0:.1f} seconds".format(ttime))

                    if self.verbose:

                        # output here
                        print('\n')
                        print('Results:\n')

                        
                    spec_params = np.zeros((n_par,7))

                    if self.verbose:
                    
                        # formatting the table
                        row_format ='{:>10}' * 8

                        # first table row
                        print(row_format.format(*['Parameter','mean','std','0.15','15.85','50.00','84.15','99.85']))

                    for p in range(n_par):
                        mean_val   = np.mean(samples[burnin:,:,p])
                        std_val    = np.std(samples[burnin:,:,p])
                        median_val = np.median(samples[burnin:,:,p])
                        ub1_val    = np.percentile(samples[burnin:,:,p],50+68.3/2)
                        lb1_val    = np.percentile(samples[burnin:,:,p],50-68.3/2)
                        ub3_val    = np.percentile(samples[burnin:,:,p],50+99.73/2)
                        lb3_val    = np.percentile(samples[burnin:,:,p],50-99.73/2)
                        spec_params[p,:] = [mean_val,std_val,lb3_val,lb1_val,median_val,ub1_val,ub3_val]

                        
                        if self.verbose:
                        
                            print(row_format.format(str(p)+':',
                                                    str('%1.2e' % mean_val),
                                                    str('%1.2e' % std_val),
                                                    str('%1.2e' % lb3_val),
                                                    str('%1.2e' % lb1_val),
                                                    str('%1.2e' % median_val),
                                                    str('%1.2e' % ub1_val),
                                                    str('%1.2e' % ub3_val)))

                    # append fitting results
                    self.fit_pars.append(spec_params)

                    # calculate rates
                    if time_div == 'tot':
                        norm_tmp = self.dataset.energies.energy_bin_wid[i]*2/self.response.flux_norm*self.dataset.times.total_time
                    elif time_div == 'none':
                        norm_tmp = self.dataset.energies.energy_bin_wid[i]*2/self.response.flux_norm
                    else:
                        norm_tmp = self.dataset.energies.energy_bin_wid[i]*2/self.response.flux_norm*time_div[i] # remove [i] if not array

                    # median value es representative for spectrum
                    self.diff_rate[i] = self.fit_pars[i][0,4]/norm_tmp

                    # 1sigma error bars
                    # upper boundary uncertainty
                    self.diff_rate_err[i,1] = self.fit_pars[i][0,5]/norm_tmp - self.diff_rate[i]
                    # lower boundary uncertainty
                    self.diff_rate_err[i,0] = self.fit_pars[i][0,3]/norm_tmp - self.diff_rate[i]


                    
            # BG-only fit
            else:

                # sample full posterior
                fit = self.model.sampling(data=self.data_per_energy_bin[i],
                                          chains=1,iter=iters,n_jobs=-1,verbose=False,
                                          pars=['Abg','model_tot','model_bg'])

                if self.verbose:
                    print('Summary for energy bin '+str(i+1)+'/'+str(self.dataset.energies.n_energy_bins)+':\n')
                    print(fit.stansummary(['Abg']))
                    
                # append fitting results
                self.fit_pars.append(fit)

            if self.verbose:
                print('###################################################################')


    def fit_spectrum(self,function,prior,guess=None,iters=2000,e_select=None,with_systematics=False):
        """
        Fit the extracted spectrum with the energy redistribution matrix, including effective area,
        with a model defined by function, and a prior.

        An initial guess for the fitted parameters can be provided in guess;
        if it is not provided, the prior is used to set initial guesses

        No checks for consistency are performed (yet)

        Can select bins to fit with e_elect
        """

        print('Fitting spectrum ... with some info included here ...')
        
        self.prior = prior
        self.function = function
        self.iters = iters
        
        self.n_par = len(self.prior)
        self.n_e = self.dataset.energies.n_energy_bins
        
        # define initial guess
        if guess == None:
            # this is not always useeful ...
            init = np.array([self.prior[i][1] for i in range(self.n_par)])
        else:
            init = np.copy(guess)

        # set initial position scale (jitter)
        init_var = init*1e-4

        if np.any(e_select == None):
            self.e_select = np.arange(self.n_e)
        else:
            self.e_select = e_select
            self.n_e = len(self.e_select)

        if with_systematics == True:
            syst_scl = self.systematics
        else:
            syst_scl = np.ones(self.n_e)
            
        # here standard emcee workflow
        ndim, nwalkers = self.n_par, self.n_par*4
        pos = [init + np.random.randn(ndim)*init_var for i in range(nwalkers)]

        start = time.time()

        self.sampler = emcee.EnsembleSampler(nwalkers,
                                             ndim,
                                             ln_posterior_COSI_spectrum,
                                             args=(self.dataset.energies.energy_bin_cen[self.e_select],
                                                   self.diff_rate[self.e_select],
                                                   self.dataset.energies.energy_bin_wid[self.e_select]*2,
                                                   np.max(np.abs(self.diff_rate_err[self.e_select,:].T),axis=0)*syst_scl[self.e_select], # maximum of asymmetric error bars
                                                   self.response.rmf[self.e_select[0]:self.e_select[-1]+1,self.e_select[0]:self.e_select[-1]+1].T,
                                                   self.response.e_min[self.e_select],
                                                   self.response.e_max[self.e_select],
                                                   self.function,
                                                   self.prior))
        
        _ = self.sampler.run_mcmc(pos, iters, progress=True)
    
        end = time.time()

        # extract samples
        self.samples = self.sampler.get_chain()
        self.samplesf = self.sampler.flatchain


        self.n_samples = self.iters#*nwalkers
        self.n_walkers = nwalkers
        self.labels = [self.prior[i][-1] for i in range(self.n_par)]

        self.burnin = int(0.5*self.n_samples)
        
        ttime = end - start
        print("Processing took {0:.1f} seconds".format(ttime))


        # output here
        print('\n')
        print('Results:\n')

        self.spec_params = np.zeros((self.n_par,7))
        
        # formatting the table
        row_format ='{:>10}' * 8

        # first table row
        print(row_format.format(*['Parameter','mean','std','0.15','15.85','50.00','84.15','99.85']))
        
        for i in range(self.n_par):
            mean_val   = np.mean(self.samples[self.burnin:,:,i])
            std_val    = np.std(self.samples[self.burnin:,:,i])
            median_val = np.median(self.samples[self.burnin:,:,i])
            ub1_val    = np.percentile(self.samples[self.burnin:,:,i],50+68.3/2)
            lb1_val    = np.percentile(self.samples[self.burnin:,:,i],50-68.3/2)
            ub3_val    = np.percentile(self.samples[self.burnin:,:,i],50+99.73/2)
            lb3_val    = np.percentile(self.samples[self.burnin:,:,i],50-99.73/2)
            self.spec_params[i,:] = [mean_val,std_val,lb3_val,lb1_val,median_val,ub1_val,ub3_val]
            
            print(row_format.format(self.labels[i]+':',
                                    str('%1.2e' % mean_val),
                                    str('%1.2e' % std_val),
                                    str('%1.2e' % lb3_val),
                                    str('%1.2e' % lb1_val),
                                    str('%1.2e' % median_val),
                                    str('%1.2e' % ub1_val),
                                    str('%1.2e' % ub3_val)))

        # chi2
        self.chi2 = COSI_model_fit(self.spec_params[:,0],
                                   self.dataset.energies.energy_bin_cen[self.e_select],
                                   self.diff_rate[self.e_select],
                                   self.dataset.energies.energy_bin_wid[self.e_select]*2,
                                   np.max(np.abs(self.diff_rate_err[self.e_select,:].T),axis=0)*syst_scl[self.e_select], # maximum af asymmetric error bars
                                   self.response.rmf[self.e_select[0]:self.e_select[-1]+1,self.e_select[0]:self.e_select[-1]+1].T,
                                   self.response.e_min[self.e_select],
                                   self.response.e_max[self.e_select],
                                   self.function,
                                   eval=False)*(-2)

        # dof
        if function == 'powerlaw_boxcar':
            self.dof = self.n_e - self.n_par + 2
        else:
            self.dof = self.n_e - self.n_par

        print('\n')
        print('Chi2 (dof): {0:.1f} ({1:d})'.format(self.chi2,self.dof))
        

    def plot_posterior_chains(self):
        """
        Plot posterior as a function of iteration
        """
        plt.plot(self.sampler.lnprobability)
        plt.xlabel('Iteration')
        plt.ylabel('Log(Posterior probability)')
        plt.xscale('log')


    def plot_parameter_chains(self,truths=None):
        """
        Plot parameters as a function of iteration
        """
        
        fig, axes = plt.subplots(self.n_par, figsize=(10, self.n_par*2.5), sharex=True)

        for i in range(self.n_par):
            ax = axes[i]
            ax.plot(np.arange(self.n_samples),self.samples[:, :, i], alpha=0.1,marker='o')
            ax.set_xlim(1, self.n_samples)
            if truths != None:
                ax.plot([1,self.n_samples],[truths[i],truths[i]],color='orange')
            ax.set_ylabel(self.labels[i])
            ax.yaxis.set_label_coords(-0.1, 0.5)
            ax.set_xscale('log')
            
        axes[-1].set_xlabel("Iteration");


    def corner_plot(self,truths=None):
        """
        Corner plot of posteriors
        """

        sigma = 68.3
        fig = corner.corner(self.samplesf[self.burnin*self.n_walkers:,:],
                            labels=self.labels,
                            #truths=np.array(params0)*scl*[ff,1,1,1],
                            quantiles=(50+sigma/np.array([-2,+2]))/100.,
                            show_titles=True,
                            bins=25,
                            fill_contours=True,
                            contourf_kwargs={"cmap": plt.cm.plasma, "colors":None},
                            levels=[1-np.exp(-4.5),1-np.exp(-2.0),1-np.exp(-0.5)],
                            #range=[(1.0,2.5),(-0.2-2.2,0.2-2.2),(150,250),(2900,3100)],
                            #       (-1.3,-0.8),(0.6,1.6),(-2.6,-1.5),
                            #       (40,65),(0,3),(0,6),
                            #       (2,2.5),(0,30),(0,5)],
                            truth_color='cyan',
                            label_kwargs={"fontsize": 8},
                            max_n_ticks=4,
                            title_kwargs={"fontsize": 8})
        for ax in fig.get_axes():
            #ax.tick_params(axis='both', which='major', labelsize=14)
            #ax.tick_params(axis='both', which='minor', labelsize=12)    
            ax.tick_params(axis='both', labelsize=14)
        fig.set_size_inches(8,8)



    def calculate_model_posteriors(self,n_use=250):
        """
        calculate model posterior in model an data space
        """

        # how many samples to use for plotting and calculation of posterior model
        #n_use = 250
        # data space
        n_plot_samples = self.n_walkers*n_use
        self.y_models = np.zeros((self.n_e,n_plot_samples))
        
        # where to evaluate model
        self.x_model = np.logspace(np.log10(100),np.log10(5500),250)
        N_model = len(self.x_model)
        self.y_modelsm = np.zeros((N_model,n_plot_samples))

        last_x_samples = self.iters-n_use

        print(self.n_walkers*last_x_samples,self.n_walkers*last_x_samples+n_plot_samples)
        
        for i in tqdm(range(self.n_walkers*last_x_samples,self.n_walkers*last_x_samples+n_plot_samples),'Loop over samples:'):
            #i = p - self.n_walkers*last_x_samples
            self.y_models[:,i-self.n_walkers*last_x_samples] = COSI_model_fit(self.samplesf[i,:],
                                                                              self.dataset.energies.energy_bin_cen[self.e_select],0,
                                                                              self.dataset.energies.energy_bin_wid[self.e_select]*2,0,
                                                                              self.response.rmf[self.e_select[0]:self.e_select[-1]+1,self.e_select[0]:self.e_select[-1]+1].T,
                                                                              self.response.e_min[self.e_select],
                                                                              self.response.e_max[self.e_select],
                                                                              self.function,          # function to fit/evaluate
                                                                              eval=True) 
    
            self.y_modelsm[:,i-self.n_walkers*last_x_samples] = globals()[self.function](self.x_model,self.samplesf[i,:])


        # auxiliary arrays
        self.ee =     self.dataset.energies.energy_bin_cen[self.e_select]
        self.ee_err = self.dataset.energies.energy_bin_wid[self.e_select]

        level = 95.4
        self.tot_model95_lo = np.percentile(self.y_modelsm, 50 - 0.5*level, axis=1 )
        self.tot_model95_hi = np.percentile(self.y_modelsm, 50 + 0.5*level, axis=1 )

        level = 68.3
        self.tot_model68_lo = np.percentile(self.y_modelsm, 50 - 0.5*level, axis=1 )
        self.tot_model68_hi = np.percentile(self.y_modelsm, 50 + 0.5*level, axis=1 )
        
        level = 0
        self.tot_model_median = np.percentile(self.y_modelsm, 50 - 0.5*level, axis=1 )



    def plot_posterior_spectrum(self,model_mode='flux',ul=None,with_systematics=False):
        """
        Plot the final result of the spectral fit
        """
        
        if model_mode == 'flux':
            edx = 0.
            ylabel_text = 'Differential Flux \n'+r'[$\mathrm{ph\,cm^{-2}\,s^{-1}\,keV^{-1}}$]'
        elif model_mode == 'eflux':
            edx = 1.
            ylabel_text = r'Flux $\times$'+ ' Differential Flux \n'+r'[$\mathrm{ph\,cm^{-2}\,s^{-1}\,keV^{-1}\,keV}$]'
        elif model_mode == 'eeflux':
            edx = 2.
            ylabel_text = r'Flux$^2$ $\times$'+ ' Differential Flux \n'+r'[$\mathrm{ph\,cm^{-2}\,s^{-1}\,keV^{-1}\,keV^{2}}$]'
        else:
            edx = 0.
            ylabel_text = 'Differential Flux \n'+r'[$\mathrm{ph\,cm^{-2}\,s^{-1}\,keV^{-1}}$]'


        if with_systematics == True:
            syst_scl = self.systematics
        else:
            syst_scl = np.ones(len(self.e_select))
        
            
        fig, ax = plt.subplots(nrows=2,ncols=1,figsize=(16,10))

        if ul == None:
            ax[0].errorbar(self.dataset.energies.energy_bin_cen[self.e_select],
                           self.diff_rate[self.e_select],
                           xerr=self.dataset.energies.energy_bin_wid[self.e_select],
                           yerr=np.abs(self.diff_rate_err[self.e_select,:].T)*syst_scl[None,self.e_select],
                           fmt='o',linewidth=3,color='black',markersize=10,
                           label=r'COSI Data Fit ($1\sigma$)')
        else:
            # mark data points which have SNR < ul with downward arrows
            #snr = self.diff_rate/np.nanmax(np.abs(self.diff_rate_err.T)*syst_scl,axis=0)
            snr = self.diff_rate[self.e_select]/np.nanmax(np.abs(self.diff_rate_err[self.e_select,:].T)*syst_scl[self.e_select],axis=0)
            ul_idx = np.where(snr < ul)[0]
            gl_idx = np.where(snr >= ul)[0]
            ax[0].errorbar(self.dataset.energies.energy_bin_cen[self.e_select][gl_idx],
                           self.diff_rate[self.e_select][gl_idx],
                           xerr=self.dataset.energies.energy_bin_wid[self.e_select][gl_idx],
                           yerr=np.abs(self.diff_rate_err[self.e_select,:].T)[:,gl_idx]*syst_scl[self.e_select][None,gl_idx],
                           fmt='o',linewidth=3,color='black',markersize=10,
                           label=r'COSI Data Fit ($1\sigma$)')
            ax[0].errorbar(self.dataset.energies.energy_bin_cen[self.e_select][ul_idx],
                           np.max(np.abs(self.diff_rate_err[self.e_select,:].T),axis=0)[ul_idx]*ul*syst_scl[self.e_select][ul_idx],
                           xerr=self.dataset.energies.energy_bin_wid[self.e_select][ul_idx],
                           yerr=np.zeros(len(ul_idx)),
                           fmt='v',linewidth=3,color='gray',markersize=15,
                           label=r'COSI Data Fit ({0}$\sigma$ upper limit)'.format(ul))

            
        self.tot_fit_model = np.median(self.y_models[:,:],axis=1)

        
        
        for i in range(self.n_e):
            level = 95.4
            ax[0].fill_between([self.ee[i]-self.ee_err[i],self.ee[i]+self.ee_err[i]],
                               np.repeat(np.percentile(self.y_models[i,:], 50 - 0.5*level),2),
                               np.repeat(np.percentile(self.y_models[i,:], 50 + 0.5*level),2),
                               color='xkcd:cobalt blue',alpha=1.0,step='mid')
    
            level = 68.3
            ax[0].fill_between([self.ee[i]-self.ee_err[i],self.ee[i]+self.ee_err[i]],
                               np.repeat(np.percentile(self.y_models[i,:], 50 - 0.5*level),2),
                               np.repeat(np.percentile(self.y_models[i,:], 50 + 0.5*level),2),
                               color='xkcd:yellowish orange',alpha=1,step='mid')


            ax[0].plot([self.ee[i]-self.ee_err[i],self.ee[i]+self.ee_err[i]],
                       np.repeat(self.tot_fit_model[i],2),
                       linewidth=4,color='gold')
            if i != self.n_e-1:
                ax[0].plot([self.ee[i]+self.ee_err[i],self.ee[i+1]-self.ee_err[i+1]],
                           [self.tot_fit_model[i],self.tot_fit_model[i+1]],
                           linewidth=4,color='gold')

        ax[0].set_xlabel('')
        ax[0].set_ylabel('Count Rate \n'+r'[$\mathrm{cts\,s^{-1}\,keV^{-1}}$]',fontsize=24)

        ax[0].set_xscale('log')
        ax[0].set_yscale('log')
        
        ax[0].set_xlim(125,5500)
        ax[0].xaxis.set_ticks([])
        #ax[0].set_ylim(5e-2,5e0)

        #ax[0].yaxis.set_ticks([1e-6,1e-5,1e-4])
        #ax[0].yaxis.set_ticklabels([r'$10^{-6}$',r'$10^{-5}$',r'$10^{-4}$'],fontsize=24)
        
        ax[0].legend(title='Data Space')


        ax[1].fill_between(self.x_model,
                           self.tot_model95_lo*self.x_model**edx,
                           self.tot_model95_hi*self.x_model**edx,
                           color='xkcd:cobalt blue',alpha=1.0,label='Total Model (95th)',zorder=1000)

        ax[1].fill_between(self.x_model,
                           self.tot_model68_lo*self.x_model**edx,
                           self.tot_model68_hi*self.x_model**edx,
                           color='xkcd:yellowish orange',alpha=1.0,label='Total Model (68th)',zorder=1000)

        #true_model = powerlaw_boxcar(x_model,params0)

        ax[1].plot(self.x_model,self.tot_model_median*self.x_model**edx,color='gold',linewidth=3,zorder=1000,label='Median model')
        #ax[1].plot(self.x_model,self.true_model,color='cyan',linewidth=3,zorder=1000,label='True model',linestyle='--')

        ax[1].set_xlabel('Energy [keV]',fontsize=24)
        ax[1].set_ylabel(ylabel_text,fontsize=24)

        ax[1].set_yscale('log')
        ax[1].set_xscale('log')

        ax[1].set_xlim(125,5500)
        #ax[1].ylim(1e-8,1e-3)

        ax[1].legend(title='Model Space')

        ax[1].xaxis.set_ticks([200,300,500,1000,2000,3000])
        ax[1].xaxis.set_ticklabels([200,300,500,1000,2000,3000],fontsize=24)

        #ax[1].yaxis.set_ticks([1e-5,1e-6,1e-7])
        #ax[1].yaxis.set_ticklabels([r'$10^{-7}$',r'$10^{-6}$',r'$10^{-5}$'],fontsize=24)

        plt.subplots_adjust(hspace=0)

        return ax
        

    def TS_map(self,grid,scipy=False,ebins='all',lookup=True):
        """
        TS: still experimental
        Creating a test statistics map from optimising a grid of point source positions.
        :param: grid: 2D grid of longitude/latitude coordinates to test above a background-only fit.
        
        Output: .TS_bg_only and .TS_vals that include the absolute values of all likelihoods.
        TS: will need to include plotting routine for illustrate the results of this call
        """

        if np.any(ebins == 'all'):
            n_e_tmp = self.dataset.energies.n_energy_bins
        else:
            n_e_tmp = 0
            for e in range(len(ebins)):
                if e == ebins[e]:
                    n_e_tmp += 1

        print(n_e_tmp)
        
        
        # tested grid
        self.grid = grid

        # make a BG-only run
        bg_only_results = fit(self.dataset,
                              self.pointings,
                              self.response,
                              background=self.background,
                              bg_only=True)#,
                              #priors=np.array([[0,1e8],[1e-6,1e-4]]))

        # BG-only optimisation
        print('Fitting background only...')
        #bg_only_results.fit()
        bg_only_results.MAP_solution(scipy=scipy,ebins=ebins)

        # save BG-only result
        if scipy == False:
            #self.TS_bg_only = np.array([np.mean(bg_only_results.fit_pars[i]['lp__']) for i in range(9)]) #fit_pars_map[0]['value']
            self.TS_bg_only = np.array([bg_only_results.fit_pars_map[i]['value']*2 for i in range(n_e_tmp)])
        else:
            self.TS_bg_only = np.array([bg_only_results.fit_pars_map[i].fun for i in range(n_e_tmp)])
        
        # init array for saving results
        self.TS_vals   = np.zeros((n_e_tmp,len(self.grid[0].ravel())))
        self.flux_vals = np.zeros((n_e_tmp,len(self.grid[0].ravel())))
        self.bg_vals   = np.zeros((n_e_tmp,len(self.grid[0].ravel()),self.background.Ncuts))
        
        # loop over grid points
        for i in tqdm(range(len(self.grid[0].ravel())),desc='Loop over grid points:'):

            if self.verbose:
                print('Now fitting object at (l,b) = (%.1f,%.1f)' % (self.grid[0].ravel()[i],self.grid[1].ravel()[i]))
                #print(self.grid[0].ravel()[i],self.grid[1].ravel()[i])


            
            # calculate response for current position
            self.response.calculate_PS_response(self.dataset,
                                                self.pointings,
                                                self.grid[0].ravel()[i],self.grid[1].ravel()[i],1,
                                                background=self.background,
                                                pixel_size=self.dataset.pixel_size,
                                                lookup=lookup)

            # if something with response entries == zero und dann checken ob es ausserhalb des FoV ist

            # fitting object for current source position
            tmp_results = fit(self.dataset,
                              self.pointings,
                              self.response,
                              background=self.background)#,
                              #priors=np.array([[0,1e8],[1e-6,1e-4]]))
                              #priors=np.array([[0,1e8],[1,1e8]]))

            try:
                # fit this position
                #tmp_results.fit()
                if scipy == False:
                    tmp_results.MAP_solution(ebins=ebins)
                    # save result if not failed
                    #self.TS_vals[:,i] = np.array([np.mean(tmp_results.fit_pars[e]['lp__']) for e in range(9)])#fit_pars_map[0]['value']
                    self.TS_vals[:,i] = np.array([tmp_results.fit_pars_map[e]['value']*2 for e in range(n_e_tmp)])
                    self.flux_vals[:,i] = np.array([tmp_results.fit_pars_map[e]['par']['flux'] for e in range(n_e_tmp)])
                    if self.background.Ncuts > 1:
                        self.bg_vals[:,i,:] = np.array([tmp_results.fit_pars_map[e]['par']['Abg'] for e in range(n_e_tmp)])
                    else:
                        self.bg_vals[:,i,0] = np.array([tmp_results.fit_pars_map[e]['par']['Abg'] for e in range(n_e_tmp)])
                else:
                    tmp_results.MAP_solution(scipy=scipy,ebins=ebins)
                    self.TS_vals[:,i] = np.array([tmp_results.fit_pars_map[e].fun for e in range(n_e_tmp)])
                    self.flux_vals[:,i] = np.array([tmp_results.fit_pars_map[e].x[0] for e in range(n_e_tmp)])
                    if self.background.Ncuts > 1:
                        self.bg_vals[:,i,:] = np.array([tmp_results.fit_pars_map[e].x[1] for e in range(n_e_tmp)])
                    else:
                        self.bg_vals[:,i,0] = np.array([tmp_results.fit_pars_map[e].x[1] for e in range(n_e_tmp)])
                    
            except RuntimeError:
                # if fit failed (zeros or something else) through RuntimeError
                if self.verbose:
                    print('Something went wrong with the fit, ignoring grid point ',self.grid[0].ravel()[i],self.grid[1].ravel()[i])
                self.TS_vals[:,i] = np.nan
        
        
    def plot_TS_map_results(self,mode='all',l_src=None,b_src=None,ebins='all'):
        """
        requires TS_map to have worked once
        """
        
        # tmporary energy bin number
        if np.any(ebins == 'all'):
            n_e_tmp = self.dataset.energies.n_energy_bins
        else:
            n_e_tmp = 0
            for e in range(len(ebins)):
                if e == ebins[e]:
                    n_e_tmp += 1
        edges_tmp = self.dataset.energies.energy_bin_edges
        
        # temporary coordinate grid arrays
        l_arr   = self.grid[0][0,:]
        l_arr[l_arr>0] -= 360 # this is fudged to get continuous maps
        b_arr   = self.grid[1][:,0]
        delta_l = np.diff(self.grid[0][0,:])[0]
        delta_b = np.diff(self.grid[1][:,0])[0]

        l_arrg = np.concatenate([l_arr-delta_l/2,np.array([l_arr[-1]+delta_l/2])])
        b_arrg = np.concatenate([b_arr-delta_b/2,np.array([b_arr[-1]+delta_b/2])])

        L_ARRg, B_ARRg = np.meshgrid(l_arrg,b_arrg)
        L_ARR, B_ARR = np.meshgrid(l_arr,b_arr)

        
        # plt maps that we got from the grid search
        if mode == 'all':

            tot_plots = n_e_tmp
            nCols = 3
            nRows = tot_plots // nCols 
            nRows += tot_plots % nCols
            
            fig,ax = plt.subplots(nrows=nRows,ncols=nCols,figsize=(nCols*5+1,nRows*4+2))

            for i in range(n_e_tmp):
                mapi = np.abs(self.TS_vals.reshape((n_e_tmp,self.grid[0].shape[0],self.grid[1].shape[1]))[i,:,:]-self.TS_bg_only[i])
                pc = ax.ravel()[i].pcolormesh(L_ARRg,B_ARRg,mapi,vmin=0,cmap=plt.cm.plasma)
                ax.ravel()[i].contour(L_ARR,B_ARR,mapi,colors='black',levels=np.nanmax(mapi)*np.array([0,0.3,0.6,0.9]))
                plt.colorbar(pc,ax=ax.ravel()[i])
                if (l_src != None) & (b_src != None):
                    ax.ravel()[i].plot(l_src,b_src,marker='*',markersize=15,color='cyan',linestyle='')

                if np.any(ebins != 'all'):
                    for l in range(self.dataset.energies.n_energy_bins):
                        if l == ebins[l]:
                            ax.ravel()[i].set_title('%i-%i keV' % (edges_tmp[l],edges_tmp[l+1]))
                else:
                    ax.ravel()[i].set_title('%i-%i keV' % (edges_tmp[i],edges_tmp[i+1]))

                ax.ravel()[i].grid(alpha=0.5,color='white')
                ax.ravel()[i].set_xlabel('Gal. Lon. [deg]')
                ax.ravel()[i].set_ylabel('Gal. Lat. [deg]')
            plt.tight_layout()
                
        elif mode == 'sum':

            fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(9,8))
            
            map_sum = np.abs(np.nansum(self.TS_vals.reshape((n_e_tmp,self.grid[0].shape[0],self.grid[1].shape[1]))[:,:,:]-self.TS_bg_only[:,None,None],axis=0))
            pc = ax.pcolormesh(L_ARRg,B_ARRg,map_sum,vmin=0,cmap=plt.cm.plasma)
            ax.contour(L_ARR,B_ARR,map_sum,colors='black',levels=np.nanmax(map_sum)*np.array([0,0.3,0.6,0.9]))
            plt.colorbar(pc,ax=ax,label='TS')
            if (l_src != None) & (b_src != None):
                ax.plot(l_src,b_src,marker='*',markersize=15,color='cyan',linestyle='')
            ax.set_title('%i-%i keV' % (edges_tmp[0],edges_tmp[-1]))
            ax.grid(alpha=0.5,color='white')
            ax.set_xlabel('Gal. Lon. [deg]')
            ax.set_ylabel('Gal. Lat. [deg]')
            
            
        else:
            print('no other modes, yet ...')

        return ax


    def plot_MAP_spectrum(self):
        plt.figure(figsize=(10.24,7.68))

        # fit
        plt.errorbar(self.dataset.energies.energy_bin_cen,
                     self.diff_rate_map,
                     xerr=self.dataset.energies.energy_bin_wid,
                     yerr=np.repeat(0,len(self.dataset.energies.energy_bin_cen)),
                     fmt='o-',label='Fit (extracted spectrum; MAP solution)',
                     linewidth=3)

        plt.xlabel('Energy [keV]')
        plt.ylabel(r'Rate [$\mathrm{cts\,s^{-1}\,keV^{-1}}$]')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlim(100,5500)
        #plt.ylim(1e-7,1e-4)

        plt.legend()


    def plot_extracted_spectrum(self,ul=None,with_systematics=False,col1='black',col2='gray'):

        if with_systematics == True:
            syst_scl = self.systematics
        else:
            syst_scl = np.ones(len(self.dataset.energies.energy_bin_cen))
        
        #plt.figure(figsize=(10.24,7.68))

        # fit
        # data
        if ul == None:
            plt.errorbar(self.dataset.energies.energy_bin_cen,
                         self.diff_rate,
                         xerr=self.dataset.energies.energy_bin_wid,
                         yerr=np.abs(self.diff_rate_err.T)*syst_scl[None,:],
                         fmt='o',linewidth=3,color=col1,markersize=10,
                         label=r'COSI Data Fit ($1\sigma$)')
        else:
            # mark data points which have SNR < ul with downward arrows
            snr = self.diff_rate/np.nanmax(np.abs(self.diff_rate_err.T)*syst_scl,axis=0)
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
        plt.ylabel(r'Rate [$\mathrm{cts\,s^{-1}\,keV^{-1}}$]')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlim(100,5500)
        #plt.ylim(1e-7,1e-4)

        plt.legend()




    def plot_rmf(self):
        """
        Plot the average rmf for the source position of interest
        """

        plt.pcolormesh(self.response.e_edges,self.response.e_edges,self.response.rmf)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Initial Energy [keV]')
        plt.ylabel('Measured Energy [keV]')
        plt.colorbar(label='$\mathrm{cm^{2}\,sr^{-1}\,keV^{-1}}$')





def COSI_model_fit(theta, x, y, dx, y_err, rsp, e_lo, e_hi, function, eval=False):
    """
    Returns:
    Negative normal-distributed log-likelihood for COSI data fitted with 
    arbitraty spectral model, accounting for the spectral response (default).
    Or, with eval=True, the output model in data space (in units of cnts/s/keV)

    Parameters:
    :param theta:     Array of to-be-fitted parameters used in combination with defined 'function'
    :param x:         Energy array (dim=(n,); in keV)
    :param y:         Differential count rate array (dim=(n,); in cnts/s/keV)
    :param dx:        Energy bin sizes for x (dim=(n,); in keV)
    :param y_err:     Uncertainties on differential count rate (one sigma values) (dim=(n,); in cnts/s/keV)
    :param rsp:       Response matrix that fits to the dimensions of x (and dx, y, y_err; dim=(m,n);, units cm2 something)
    :param e_lo:      Lower interpolated energy edges for response calculation(dim=(m,); in keV)
    :param e_hi:      Upper interpolated energy edges for response calculation(dim=(m,); in keV)  
    :param function:  String of named function 
    """

    # Dec 10 2020; TS: here, I need to include the factor that I don't know yet: fudge factor normalisation
        
    # Integrate model with Simpson's rule over the interpolated energy bins
    # this works in general for smooth models, but not for cuts, steps or models 'within one bin', etc.
    if function != 'powerlaw_boxcar':
        integrated_model = (e_hi-e_lo)/6.0*(globals()[function](e_lo,theta)+
                                            4*globals()[function]((e_lo+e_hi)/2.0,theta)+
                                            globals()[function](e_hi,theta))
    else:
        integrated_model = integrate_powerlaw_boxcar(e_lo,e_hi,theta)

    # Apply response matrix
    folded_model = np.dot(integrated_model,rsp)
 
    # Return to differential model
    folded_differential_model = folded_model / dx

    # Evaluate either chi2
    if eval==False:
        return -0.5*np.nansum((y-folded_differential_model)**2/y_err**2)
    # or return the folded model itself at a certain set of parameters
    else:
        return folded_differential_model



def ln_prior_COSI_spectrum(theta,prior):
    """
    this function assumes prior dictionaries of the shape:
    prior = {0: ('normal_prior',1e-4,1e-5,'Amplitude'),
             1: ('uniform_prior',-1.8,-1.6,'Index')}
    with the key being numbered sequentially
    then there is a tag for the prior function (normal, uniform, whatever you define else)
    then the prior values as defined in the functions
    and a name to identifiy later
    NO CHECKS FOR CONSISTENCY ARE PERFORMED
    """
        
    # add prior for each parameter to be fitted according to definition in 'function'
        
    lnprior = 0.
    for p in range(len(theta)):
        lnprior += globals()[prior[p][0]](theta[p],prior[p][1:])
    return lnprior


def ln_posterior_COSI_spectrum(theta, x, y, dx, y_err, rsp, e_lo, e_hi, function, prior):
    """
    Combining likelihood (chi2) and prior information for full posterior asdf asdf asdf asdf 
    """
        
    lp = ln_prior_COSI_spectrum(theta,prior)
    if not np.isfinite(lp):
        return -np.inf
    return lp + COSI_model_fit(theta, x, y, dx, y_err, rsp, e_lo, e_hi, function)




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




def COSImodfit_bgonly(theta, data, background_model, eval=False):
    """
    Returns:
    cash-stat for model fit
    
    one-parameter background model used only

    Parameters:
    :param theta:     Array of to-be-fitted parameters 
    :param other: ...
    """

    total_model = theta[0]*background_model
    
    stat = -2*np.sum(total_model - data*(np.log(total_model)))
    if np.isnan(stat) | (theta[0] < 0):
        return -np.inf
    else:
        return stat
