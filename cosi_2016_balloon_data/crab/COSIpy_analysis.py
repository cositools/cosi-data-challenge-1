from COSIpy import *
import response
from fit import *
from COSIpy_tools import step_plot
import timeit
# EFF AREA
import mhealpy as hp
from mhealpy import HealpixMap, HealpixBase
from histpy import Histogram,Axes,Axis
from sparse import COO
import numpy as np
import ROOT as M
import math

M.gSystem.Load("$(MEGALIB)/lib/libMEGAlib.so")

G = M.MGlobal()
G.Initialize()
###

start = timeit.default_timer()

# define variables
data_dir = 'simulated_data' # directory of source data set (simulated data with all 4 sources)
filename = 'DC1_Combined.inc1.id1.extracted.extracted6.tra.gz' # data set file name (simulated data with all sources & bg)
response_filename = 'response_files/ContinuumResponse.npz' # response file name
background_filename = 'simulated_data/Scaled_Ling_BG.inc1.id1.extracted.npz' # simulated background file name
chris_background_filename = 'Scaled_Ling_BG.inc1.id1.extracted.extracted.tra.gz'
source_name = 'CygX1_sim3' # name of source (simulated data)
l,b = 71.334998, 3.06683 # longitude & latitude of source (Cyg X-1)
Delta_T = 7200 # time bin size
energy_bin_edges = np.array([150,  220,  325,  480,  520,  765, 1120, 1650, 2350, 3450, 5000]) # definition of energy bins
pixel_size = 6. # pixel size
background_mode = 'from file' # Background model from file
ul = 3 # SNR limit
prior = {0: ('truncated_normal_prior',1e-9,1e8,0,np.inf,r'$C_0$'), 1: ('normal_prior',-2.0,5.0,r'$\alpha$')} # prior
function = 'powerlaw' # function
bs = 60
num = 14

# read in & bin data
def read(data_dir,filename,Delta_T,energy_bin_edges,pixel_size):
	print('Reading in data...')
	analysis = COSIpy(data_dir,filename) # create analysis object
	analysis.read_COSI_DataSet() # read in data
	print('Binning data...')
	analysis.dataset.time_binning_tags(time_bin_size=Delta_T) # define time bins
	analysis.dataset.init_binning(energy_bin_edges=energy_bin_edges,pixel_size=pixel_size) # define energy and pixel binning
	analysis.dataset.get_binned_data() # bin data
	print('Number of time bins:',analysis.dataset.times.n_ph)
	return analysis

# plot spectrum & light curve
def plot_spec_light(analysis,source_name):
	print('Making plots...')
	analysis.dataset.plot_raw_spectrum() # plot spectrum of raw counts
	plt.xscale('log')
	plt.savefig('plots/spectrum_' + source_name + '.png',dpi=500)
	plt.clf()
	analysis.dataset.plot_lightcurve() # plot light curve
	plt.savefig('plots/lightcurve_' + source_name + '.png',dpi=500)
	plt.clf()

# define & plot pointings, plot elevation
def pointings(analysis,source_name,l,b):
	print('Defining pointings...')
	pointing = Pointing(dataset=analysis.dataset) # definition of pointings
	print('Making plots...')
	plt.plot(pointing.zpoins[:,0]+360,pointing.zpoins[:,1],'o') # create plot of pointings
	plt.plot(l,b,'*r',markersize=10)
	plt.xlabel('Longitude [deg]')
	plt.ylabel('Latitude [deg]')
	plt.savefig('plots/pointings_' + source_name + '.png',dpi=500)
	plt.clf()
	analysis.plot_elevation([l],[b],[source_name]) # create plot of elevation
	tot_time = 0
	tmp_elevation = analysis.horizon-angular_distance(l,b,analysis.dataset.l_pointing,analysis.dataset.b_pointing)
	for i in range(len(tmp_elevation)):
		if tmp_elevation[i] > 0:
			tot_time += Delta_T
	plt.savefig('plots/elevation_' + source_name + '.png',dpi=500)
	plt.clf()
	return pointing, tot_time

# read in background & response
def bkgd_resp(analysis,background_mode,response_filename,background_filename,pointing,l,b,pixel_size):
	print('Reading in background...')
	tracer = np.sum(analysis.dataset.binned_data,axis=(1,2,3)) # define tracer
	tracer = tracer/np.mean(tracer)
	if background_mode == 'from file':
		background1 = BG(dataset=analysis.dataset,mode=background_mode,filename=background_filename) # read in background (without tracer)
		background2 = BG(dataset=analysis.dataset,mode=background_mode,filename=background_filename,tracer=tracer) # read in background (with tracer)
	else:
	        background1 = BG(dataset=analysis.dataset,mode=background_mode) # read in background (without tracer)
	        background2 = BG(dataset=analysis.dataset,mode=background_mode,tracer=tracer) # read in background (with tracer)
	print('Reading in response...')
	rsp = response.SkyResponse(filename=response_filename,pixel_size=pixel_size) # read in response
	rsp.calculate_PS_response(analysis.dataset,pointing,l,b,1,background=background1,pixel_size=pixel_size,lookup=False)
	print('Making plots...')
	plt.plot(np.sum(analysis.dataset.binned_data[:,0,:,:],axis=(1,2))) # plot binned data light curve for first energy bin
	plt.savefig('plots/binnedlightcurve_' + source_name + '.png',dpi=500)
	plt.clf()   
	plt.plot(np.sum(analysis.dataset.binned_data[:,1,:,:],axis=(1,2))) # plot binned data light curve for second energy bin
	plt.plot(np.sum(background2.bg_model_reduced[1],axis=1)) # plot background response (first guess at background)
	plt.plot(np.sum(rsp.sky_response[1],axis=1)*1000) # plot sky response (first guess at source contribution)
	plt.savefig('plots/responses_' + source_name + '.png',dpi=500)
	plt.clf()
	return background1, background2, rsp

# fitting
def fitting(analysis,rsp,pointing,background,l,b,pixel_size):
	print('Creating fitting object...')
	result = fit(analysis.dataset,pointing,rsp,background,verbose=True) # create fitting object
	return result

# fit spectrum
def fit_spectrum(analysis,source_name,result,pointing,rsp,background,ul,function,prior):
	print('Fitting...')
	result.MAP_solution(scipy=False) # perform optimization of joint posterior distribution
	result.fit(iters=2000,use_emcee=True) # fit with pointing definition, background model, & point source response

# fit plots
def fit_plots(analysis,source_name,result,pointing,rsp,background,ul,function,prior):
	print('Making plots...')
	result.plot_MAP_spectrum() # plot MAP spectrum
	plt.savefig('plots/MAPspectrum_' + source_name + '.png',dpi=500)
	plt.clf()
	pc = result.plot_extracted_spectrum(ul=ul) # plot extracted spectrum
	plt.xlim(150,3450)
	plt.savefig('plots/extractedspectrum_' + source_name + '.png',dpi=500)
	plt.clf()
	result.fit_spectrum(function,prior,iters=10000,e_select=np.arange(0,9),with_systematics=False) # fit extended spectrum with energy redistribution matrix
	result.calculate_model_posteriors() # calculate model posterior
	result.plot_posterior_spectrum(ul=ul,with_systematics=False) # plot final result of spectral fit
	plt.savefig('plots/posteriorspectrum_' + source_name + '.png',dpi=500)
	plt.clf()
	result.plot_parameter_chains() # plot parameters as function of iteration
	plt.savefig('plots/parameterchains_' + source_name + '.png',dpi=500)
	plt.clf()
	result.corner_plot() # create corner plot of posteriors
	plt.savefig('plots/cornerplot_' + source_name + '.png',dpi=500)
	plt.clf()
	result.plot_posterior_chains() # plot posterior as function of iteration
	plt.savefig('plots/posteriorchains_' + source_name + '.png',dpi=500)
	plt.clf()
	return pc

# TS map
def tsmap(analysis,source_name,result,l,b,bs,num,background):
	l_gridg = np.linspace(l-bs,l+bs,num)
	b_gridg = np.linspace(b-bs,b+bs,num)
	l_grid = l_gridg[0:-1]+np.diff(l_gridg[0:-1])[0]/2
	b_grid = b_gridg[0:-1]+np.diff(b_gridg[0:-1])[0]/2
	l_grid[l_grid<-180] += 360
	L_GRIDg, B_GRIDg = np.meshgrid(l_gridg,b_gridg)
	L_GRID, B_GRID = np.meshgrid(l_grid,b_grid)
	print('Creating TS map...')
	result.TS_map([L_GRID,B_GRID],scipy=False,ebins=[-1,-1,2,-1,-1,-1,-1,-1,-1,-1],lookup=False) # create test statistics map
	print('Making plots...')
	result.plot_TS_map_results(ebins=[-1,-1,2,-1,-1,-1,-1,-1,-1,-1],l_src=l-360,b_src=b) # plot test statistics map
	plt.savefig('plots/tsmap_' + source_name + '.png',dpi=500)
	plt.clf()

def integral_func(E):
        return 10**(-3) * (E/100)**(-2.2)

def band_func(E, a, b, E0, A):
        return A * ((a-b)*E0/100)**(a-b) * math.exp(b-a) * (E/100)**b

# back of envelope calculation
def boe():
	print('Doing BoE calculation...')
	result.MAP_solution(scipy=True)
	dr = Histogram.open("simulated_data/COSI2016_6deg_10ebins.p1.binnedimaging.imagingresponse.area.nside16.h5")
	theta = 0.35
	phi = 0
	healpix_axis = HealpixBase(nside = hp.npix2nside(dr.axes['NuLambda'].nbins))
	skypix = healpix_axis.ang2pix(theta,phi)
	ax,plot = dr.slice[{'NuLambda': skypix}].project('Ei').plot(errorbars = False)
	plt.savefig('plots/effarea.png',dpi=500)
	plt.clf()
	analysis_bkgd = read(data_dir,chris_background_filename,Delta_T,energy_bin_edges,pixel_size)
	energy_bins = np.array([185,  272.5,  402.5,  500,  642.5,  942.5, 1385, 2000, 2625, 4225]) # definition of energy bins
	widths = np.array([70, 105, 155, 40, 245, 355, 530, 700, 550, 1550])        
	f_arr = np.array([])
	int_arr = np.array([])
	int_x_arr = np.array([])
	int_arr_chris = np.array([])
	tot_f = 0
	for i in range(150,5000):
        	int_x_arr = np.append(int_x_arr, i)
        	int_arr = np.append(int_arr, band_func(i,-1.99,-2.32,531,7.52e-4))
	for i in range(100,50000):
        	int_arr_chris = np.append(int_arr_chris, band_func(i,-1.99,-2.32,531,7.52e-4))
	for i in range(len(energy_bin_edges)-1):
        	N = np.sum(np.sum(analysis.dataset.binned_data[:,i,:,:],axis=(1,2)))
        	B = np.sum(np.sum(analysis_bkgd.dataset.binned_data[:,i,:,:],axis=(1,2)))
        	E = 0
        	for j in range(analysis.dataset.times.n_ph):
                	E += dr.slice[{'NuLambda': skypix}].project('Ei')[i] * Delta_T
        	F = (N - B)/E
        	f_arr = np.append(f_arr, F/widths[i])
        	print(i, N, B, E)
	result.plot_MAP_spectrum()
	plt.bar(energy_bins, f_arr, width=widths, color='white', edgecolor='green')
	plt.plot(int_x_arr, int_arr, color='orange')
	print(np.sum(int_arr))
	print(np.sum(int_arr_chris))
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel('Energy (keV)')
	plt.ylabel('ph/cm^2/s/keV')
	plt.savefig('plots/flux_spectrum6.png',dpi=500)

def compare_spectrum(data_dir,Delta_T,energy_bin_edges,pixel_size,source_name,background_mode,response_filename,background_filename,l,b):
	analysis1 = read(data_dir,'DC1_Combined.inc1.id1.extracted.extracted6.tra.gz',Delta_T,energy_bin_edges,pixel_size)
	analysis2 = read(data_dir,'GalacticScan.cygX1_and_BG.inc1.id1.extracted.tra.gz',Delta_T,energy_bin_edges,pixel_size)
	analysis3 = read(data_dir,'GalacticScan.cygX1_only.inc1.id1.extracted2.tra.gz',Delta_T,energy_bin_edges,pixel_size)
	pointing1,tot_time1 = pointings(analysis1,source_name,l,b)
	pointing2,tot_time2 =  pointings(analysis2,source_name,l,b)
	background1_1,background2_1,rsp_1 = bkgd_resp(analysis1,background_mode,response_filename,background_filename,pointing1,l,b,pixel_size)
	background1_2,background2_2,rsp_2 = bkgd_resp(analysis2,background_mode,response_filename,background_filename,pointing2,l,b,pixel_size)
	result1 = fitting(analysis1,rsp_1,pointing1,background1_1,l,b,pixel_size)
	result2 = fitting(analysis2,rsp_2,pointing2,background1_2,l,b,pixel_size)
	result1.fit(iters=2000,use_emcee=True)
	result2.fit(iters=2000,use_emcee=True)
	analysis3.dataset.plot_raw_spectrum()
	result1.plot_extracted_spectrum(ul=ul,col1='black',col2='black') # plot extracted spectrum
	result2.plot_extracted_spectrum(ul=ul,col1='blue',col2='blue') # plot extracted spectrum
	plt.xlim(150,5000)
	plt.ylabel('Counts [cnts/keV]')
	plt.savefig('plots/extractedspectrum_' + source_name + '.png',dpi=500)

def orig_ana(data_dir,filename,Delta_T,energy_bin_edges,pixel_size,source_name,background_mode,response_filename,l,b,ul,function,prior,bs,num):
	analysis = read(data_dir,filename,Delta_T,energy_bin_edges,pixel_size)
	plot_spec_light(analysis,source_name)
	pointing,tot_time = pointings(analysis,source_name,l,b)
	background1,background2,rsp = bkgd_resp(analysis,background_mode,response_filename,background_filename,pointing,l,b,pixel_size)
	result = fitting(analysis,rsp,pointing,background1,l,b,pixel_size)
	fit_spectrum(analysis,source_name,result,pointing,rsp,background1,ul,function,prior)
	pc = fit_plots(analysis,source_name,result,pointing,rsp,background1,ul,function,prior)
	#tsmap(analysis,source_name,result,l,b,bs,num,background1)

#orig_ana(data_dir,filename,Delta_T,energy_bin_edges,pixel_size,source_name,background_mode,response_filename,l,b,ul,function,prior,bs,num)
compare_spectrum(data_dir,Delta_T,energy_bin_edges,pixel_size,source_name,background_mode,response_filename,background_filename,l,b)

stop = timeit.default_timer()
print('Time: ' + '{:.2f}'.format((stop - start)/60) + ' mins')



