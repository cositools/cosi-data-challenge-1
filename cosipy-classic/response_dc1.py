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
from COSIpy_dc1 import FISBEL
from COSIpy_dc1 import dataset
from COSIpy_dc1 import GreatCircle
from COSIpy_dc1 import angular_distance
from COSIpy_dc1 import find_nearest

deg2rad = np.pi/180

class SkyResponse:

    def __init__(self,
                 filename,
                 pixel_size,
                 from_saved_file=True,
                 energy_bin_edges=np.array([506,516]),
                 keep_everything=False,
                 verbose=False): # energy bin edges not used right now (read in through MEGAlib file)
        
        self.filename = filename
        self.pixel_size = pixel_size
        # the energy bins that will be used for the analysis, NOT (necessarily) the ones defined in the rsp file
        self.energy_bin_edges = energy_bin_edges
        # initialise empty dataset for the sky response to use FISBEL and other binning
        self.rsp = dataset(name='SkyResponse')
        # init binning according the how the response is specified
        self.rsp.init_binning(pixel_size=self.pixel_size,energy_bin_edges=self.energy_bin_edges)

        if from_saved_file:
            # read everything in (very large in many cases)
            print('Reading complete continuum response.')
            self.LoadRegularBinnedMEGAlibResponse()
            print('Done.\n')

            # reduced energy matrix
            self.CDS_Summed_Response()

            # remove full matrix after it is read in and built for IRF and RMF (default: do that because it is huge)
            if not keep_everything:
                self.ReduceToImportantParts()


    def CDS_Summed_Response(self):

        self.rsp.ZA_energy_response = 0
        
        if self.n_e == 1:

            print('You probably used a line response, there is no RMF.')

        elif self.n_e > 1:
            
        #log.info('Creating general RMF matrices, stay tuned.')
            print('Creating general RMF matrices.')
            # sum over CDS (phi and psi/chi) axis 2 and 3 of matrix:
            # this allows quicker weightings for the RMF later

            for i in tqdm(range(self.rsp.phis.n_phi_bins),'Loop over phi bins:'):
                self.rsp.ZA_energy_response += np.sum(self.rsp.response_grid_normed[:,:,i,:,:,:],2)

        else:
            print('This should not happen.')
                
        print('Done.\n')

        
    def ReduceToImportantParts(self):

        # for IRF, sum over initial energies (axis 4 if not reduced data space)
        print('Creating general IRF.')
        
        if self.n_e == 1:
            self.rsp.response_grid_normed_efinal = np.copy(self.rsp.response_grid_normed)

        # multiple energy bins:
        elif self.n_e > 1:
            self.rsp.response_grid_normed_efinal = np.sum(self.rsp.response_grid_normed,axis=4)

        else:
            print('This should not happen, did you read in a response file?')

        print('Done.\n')

        # and delete the full matrix from the rsp object to save space
        # yes, bad style
        #@log.info('Deleting full matrix.')
        print('Deleting full matrix.')
        del self.rsp.response_grid_normed
        print('Done.')
            




    def LoadRegularBinnedMEGAlibResponse(self):

        try:

            with np.load(self.filename) as content:
                self.rsp.response_grid_normed = content['ResponseGrid']
                self.e_cen = content['e_cen']
                self.e_wid = content['e_wid']
                self.e_edges = content['e_edges']
                self.e_max = content['e_max']
                self.e_min = content['e_min']
                self.n_e = content['n_e']
                self.l_cen = content['l_cen']
                self.l_wid = content['l_wid']
                self.l_edges = content['l_edges']
                self.l_max = content['l_max']
                self.l_min = content['l_min']
                self.n_l = content['n_l']
                self.b_cen = content['b_cen']
                self.b_wid = content['b_wid']
                self.b_edges = content['b_edges']
                self.b_max = content['b_max']
                self.b_min = content['b_min']
                self.n_b = content['n_b']
                self.L_ARR = content['L_ARR']
                self.B_ARR = content['B_ARR']
                self.L_ARR_edges = content['L_ARR_edges']
                self.B_ARR_edges = content['B_ARR_edges']
                self.dL_ARR = content['dL_ARR']
                self.dB_ARR = content['dB_ARR']
                self.dL_ARR_edges = content['dL_ARR_edges']
                self.dB_ARR_edges = content['dB_ARR_edges']
                self.dOmega = content['dOmega']

        except FileNotFoundError:
            print('File '+str(self.filename)+' not found or is not a regular binned response array.')


    
    def calculate_PS_response(self,
                              dataset,
                              pointings,
                              l_src,b_src,flux_norm,
                              reduced=True,
                              pixel_size=5.,
                              cut=60.,
                              background=None,
                              lookup=True,
                              verbose=False):
        self.l_src = l_src
        self.b_src = b_src
        self.flux_norm = flux_norm
        self.verbose = verbose

        if self.verbose:
            print('Calculating zeniths and azimuths for all pointings ...')
            
        zens,azis = zenazi(pointings.xpoins[:,0],pointings.xpoins[:,1],
                           pointings.ypoins[:,0],pointings.ypoins[:,1],
                           pointings.zpoins[:,0],pointings.zpoins[:,1],
                           self.l_src,self.b_src)

        if self.verbose:
            print('Done.\n')

            print('Calculating CDS count expectations for all bins ...')
        # initialise sky response list to include all energies
        self.sky_response = []

        if reduced:
            
            for i in tqdm(range(dataset.energies.n_energy_bins),'Loop over energy bins: '):
                    
                # reshape background model to reduce CDS if possible
                # this combines the 3 CDS angles into a 1D array for all times at the chosen energy
                bg_tmp = background.bg_model[:,i,:,:].reshape(dataset.times.n_ph,self.rsp.n_phi_bins*self.rsp.n_fisbel_bins)
                
                # get indices of where no entries are there at all (will always be zero, so can be ignored)
                calc_this = np.where(np.sum(bg_tmp,axis=0) != 0)[0]

                # reshape response grid the same way and choose only non-zero indices
                # one energy (or same response for each chosen bin)
                if self.n_e == 1:
                    rsp_tmp = self.rsp.response_grid_normed_efinal.reshape(self.n_b,self.n_l,self.rsp.n_phi_bins*self.rsp.n_fisbel_bins)[:,:,calc_this]

                # multiple energy bins:
                # choose nearest neighbour (center to center) to get response (TS: interpolation maybe later? how to if only three bands and one is a strong line?)
                elif self.n_e > 1:
                    rsp_idx = find_nearest(self.e_cen,dataset.energies.energy_bin_cen[i])
                    # sum over initial energy axis already to get entry for measured energy with all initials possible

                    rsp_tmp = self.rsp.response_grid_normed_efinal.reshape(self.n_b,self.n_l,self.rsp.n_phi_bins*self.rsp.n_fisbel_bins,self.n_e)[:,:,calc_this,rsp_idx]#[:,:,:,rsp_idx]
                    
                # shouldnt happen
                else:
                    if self.verbose:
                        print('Something went wrong ...')
                    

                # sky response per pointing (weighted by (small) time interval in pointings to get counts
                # zenith pixel size at time:
                # conversion from fisbel to regular includes a factor sin(zenith): otherwise nrow of fisbel will get lost in the conversion (???)
                sky_response_pp = get_response_with_weights(rsp_tmp,zens,azis,cut=cut,binsize=pixel_size,lookup=lookup)*pointings.dtpoins[:,None]#/np.sin(zens*deg2rad)[:,None]#*zenith_pixel_size[:,None]
                sky_response_pp[np.isnan(sky_response_pp)] = 0.

                # pre-define response array per time bin to fill
                sky_response_hh = np.zeros((dataset.times.n_ph,len(calc_this)))
    
                # loop until all defined time bins of previous definition are included
                for c in range(dataset.times.n_ph):
                    cdx = np.where((pointings.cdtpoins > dataset.times.times_min[dataset.times.n_ph_dx[c]]) &
                                   (pointings.cdtpoins <= dataset.times.times_max[dataset.times.n_ph_dx[c]]))[0]
                        
                    sky_response_hh[c,:] = np.sum(sky_response_pp[cdx,:],axis=0)


                # normalise response to total effective area?
                sky_response_hh /= np.sum(sky_response_hh)

                # calculate sky model count expectaion
                self.sky_response.append(sky_response_hh*self.flux_norm)

            if self.verbose:
                print('Done.\n')
            
            if self.verbose:
                print('Calculating averaged RMF for object at (l,b) = (%.1f,%.1f)' % (self.l_src,self.b_src))

            if self.n_e == 1:
                print('No RMF still.')

            elif self.n_e > 1:
                
                # zenith indices of response
                zidx = np.floor(zens/dataset.pixel_size).astype(int)
                # azimuth indices of response
                aidx = np.floor(azis/dataset.pixel_size).astype(int)
            
                # remove out of bounds indices
                weights = np.ones(len(zidx))
                zidx[zidx < 0] = 0.
                weights[zidx < 0] = 0.
                aidx[aidx < 0] = 0.
                weights[aidx < 0] = 0.
            
                # energy normalisation matrix
                erg_mat = np.meshgrid(self.e_wid,self.e_wid)
            
                # weighting the response at each pointing
                self.rmf = 0

                for n in tqdm(range(len(pointings.dtpoins)),'Loop over pointings:'):
                    self.rmf += self.rsp.ZA_energy_response[zidx[n],aidx[n],:,:].T/erg_mat[0]*pointings.dtpoins[n]*weights[n]


        else:

            for i in range(dataset.energies.n_energy_bins):

                # reshape response grid the same way and choose only non-zero indices
                # one energy (or same response for each chosen bin)
                if self.n_e == 1:
                    rsp_tmp = self.rsp.response_grid_normed.reshape(self.n_b,self.n_l,self.rsp.n_phi_bins*self.rsp.n_fisbel_bins)[:,:,calc_this]


                # multiple energy bins:
                # choose nearest neighbour (center to center) to get response (TS: interpolation maybe later? how to if only three bands and one is a strong line?)
                elif self.n_e > 1:
                    rsp_idx = find_nearest(self.e_cen,dataset.energies.energy_bin_cen[i])
                    # sum over initial energy axis already to get entry for measured energy with all initials possible
                    rsp_tmp = self.rsp.response_grid_normed_efinal.reshape(self.n_b,self.n_l,self.rsp.n_phi_bins*self.rsp.n_fisbel_bins,self.n_e)[:,:,:,rsp_idx]

                # shouldnt happen
                else:
                    if self.verbose:
                        print('Something went wrong ...')


                # sky response per pointing (weighted by (small) time interval in pointings to get counts
                sky_response_pp = get_response_with_weights(rsp_tmp,zens,azis,cut=60,binsize=pixel_size,lookup=lookup)*pointings.dtpoins[:,None]

                # pre-define response array per time bin to fill
                sky_response_hh = np.zeros((dataset.times.n_ph,self.rsp.n_phi_bins*self.rsp.n_fisbel_bins))
                
                # loop until all defined time bins of previous definition are included
                for c in range(dataset.times.n_ph):
                    cdx = np.where((pointings.cdtpoins > dataset.times.times_min[dataset.times.n_ph_dx[c]]) &
                                   (pointings.cdtpoins <= dataset.times.times_max[dataset.times.n_ph_dx[c]]))[0]

                    sky_response_hh[c,:] = np.sum(sky_response_pp[cdx,:],axis=0)

                # calculate sky model count expectaion
                sky_response_hh = sky_response_hh.reshape(dataset.times.n_ph,self.rsp.n_phi_bins,self.rsp.n_fisbel_bins)
                
                # normalise response to total effective area?
                sky_response_hh /= np.sum(sky_response_hh)
                
                self.sky_response.append(sky_response_hh)

            sky_rsp_full_tmp = np.zeros((dataset.times.n_ph,dataset.energies.n_energy_bins,self.rsp.n_phi_bins,self.rsp.n_fisbel_bins))
            
            for i in range(dataset.energies.n_energy_bins):
                sky_rsp_full_tmp[:,i,:,:] = self.sky_response[i]

            self.sky_response = sky_rsp_full_tmp

            if self.verbose:
                print('Calculating averaged RMF for object at (l,b) = (%.1f,%.1f)' % (self.l_src,self.b_src))
            # zenith indices of response
            zidx = np.floor(zens/dataset.pixel_size).astype(int)
            # azimuth indices of response
            aidx = np.floor(azis/dataset.pixel_size).astype(int)
                
            # remove out of bounds indices
            weights = np.ones(len(zidx))
            zidx[zidx < 0] = 0.
            weights[zidx < 0] = 0.
            aidx[aidx < 0] = 0.
            weights[aidx < 0] = 0.

            # energy normalisation matrix
            erg_mat = np.meshgrid(self.e_wid,self.e_wid)

            # weighting the response at each pointing
            self.rmf = 0

            for n in tqdm(range(len(pointings.dtpoins)),'Loop over pointings:'):
                self.rmf += self.rsp.ZA_energy_response[zidx[n],aidx[n],:,:].T/erg_mat[0]*pointings.dtpoins[n]*weights[n]
            
                

    def plot_CDS_response(self,
                          phi=None,psi=None,chi=None,
                          zen=None,azi=None,erg=0):

        # dimension of response
        dim = len(self.rsp.response_grid_normed_efinal.shape)
        
        if (phi != None) & (psi != None) & (chi != None):
            print('Plotting Compton response from (phi/psi/chi) = (%.1f/%.1f/%.1f) to (Z/A):' % (phi,psi,chi))
        
            idx = find_CDS_indices(phi,
                                   psi-90.,
                                   chi,
                                   np.rad2deg(self.rsp.phis.phi_cen),
                                   np.rad2deg(self.rsp.fisbels.lat_cen)-90.,
                                   np.rad2deg(self.rsp.fisbels.lon_cen))
        
            #print(idx)
            if dim == 4:
                plt.pcolormesh(np.rad2deg(self.L_ARR),
                               np.rad2deg(self.B_ARR),
                               self.rsp.response_grid_normed_efinal[:,:,idx[0],idx[1]])
            elif dim == 5:
                plt.pcolormesh(np.rad2deg(self.L_ARR),
                               np.rad2deg(self.B_ARR),
                               self.rsp.response_grid_normed_efinal[:,:,idx[0],idx[1],erg])
            else:
                print('not a response?')
                
            plt.colorbar(label=r'$\mathrm{ph\,cm^{2}\,sr^{-1}}$')
            plt.xlabel('Azimuth [deg]')
            plt.ylabel('Zenith [deg]')
        
        elif (zen != None) & (azi != None) & (phi != None):
            print('Plotting Compton response from (Z/A) = (%.1f/%.1f) to (%.1f/psi/chi):' % (zen,azi,phi))
            print('This might take a little ...')
        
            idx = find_grid_indices(zen,
                                    azi,
                                    np.rad2deg(self.b_cen),
                                    np.rad2deg(self.l_cen))
            #print(idx)
        
            idx_phi = find_grid_indices(phi,
                                        0,
                                        np.rad2deg(self.rsp.phis.phi_cen),
                                        0)
        
            #print(idx_phi)
            if dim == 4:
                self.rsp.fisbels.plot_FISBEL_tessellation(values=self.rsp.response_grid_normed_efinal[idx[0],idx[1],idx_phi[0],:],
                                                          colorbar=True,
                                                          tiles=True,
                                                          deg=True)
            elif dim == 5:
                self.rsp.fisbels.plot_FISBEL_tessellation(values=self.rsp.response_grid_normed_efinal[idx[0],idx[1],idx_phi[0],:,erg],
                                                          colorbar=True,
                                                          tiles=True,
                                                          deg=True)
            else:
                print('not a response?')

            plt.xlabel('Chi local [deg]')
            plt.ylabel('180 deg - Psi local [deg]')
        
        else:
            print('Need to define (phi/psi/chi) or (zen/azi/phi) for plot.')

            

def get_response_with_weights(Response,zenith,azimuth,deg=True,binsize=5,cut=60.0,lookup=True):
    """
    Calculate response at given zenith/azimuth position of a source relative to COSI,
    using the angular distance to the 4 neighbouring pixels that overlap using a 
    certain binsize.
    Note that this also introduces some smoothing of the response as zeniths/azimuths
    on the edges or corners of pixels will not be weighted equally but still get 
    contributions from the remaining pixels.
    :param: Response      Response grid with regular sky pixel dimension (zenith x azimuth)
                          and unfolded 1D phi-psi-chi dimension.
                          Optional with energy redistribution matrix included.
    :param: zenith        Zenith positions of the source with respect to the instrument (in deg)
    :param: azimuth       Azimuth positions of the source with respect to the instrument (in deg)
    :option: deg          Default True, (right now not checked for any purpose)
    :option: binsize      Default 5 deg (matching the sky dimension of the response). If set
                          differently, make sure it matches the sky dimension as otherwise,
                          false results may be returned
    :option: cut          Threshold to cut the response calculation after a certain zenith angle.
                          Default 60 deg (~ COSI FoV)
    :option: lookup       Use only pixel that got hit to calculate response (default: True)
    Returns an array of length equal the response that is input.
    """

    if lookup == True:
        # look up, no weighting
        rsp_mean = get_response_from_pixelhit_vector(Response,zenith,azimuth,binsize=binsize,cut=cut)

    else:
    
        # calculate the weighting for neighbouring pixels using their angular distance
        # also returns the indices of which pixels to be used for response averaging
        widx = get_response_weights_vector(zenith,azimuth,binsize,cut=cut)
        # This is a vectorised function so that each entry gets its own weighting
        # at the correct positions of the input angles ([:, None] is the same as 
        # column-vector multiplcation of a lot of ones)

        # check for negative weights and indices and remove
        widx[1][widx[0][:,0,:] < 0] = 0.
        widx[1][widx[0][:,1,:] < 0] = 0.
        for i in range(4):
            widx[0][i,0,widx[0][i,0,:] < 0] = 0.
            widx[0][i,1,widx[0][i,1,:] < 0] = 0.
    
        # one energy bin
        #print(Response.shape,len(Response.shape))
        if len(Response.shape) < 4:
            rsp0 = Response[widx[0][0,1,:],widx[0][0,0,:],:]*widx[1][0,:][:, None]
            rsp1 = Response[widx[0][1,1,:],widx[0][1,0,:],:]*widx[1][1,:][:, None]
            rsp2 = Response[widx[0][2,1,:],widx[0][2,0,:],:]*widx[1][2,:][:, None]
            rsp3 = Response[widx[0][3,1,:],widx[0][3,0,:],:]*widx[1][3,:][:, None]
            # with energy matrix included
        elif len(Response.shape) >= 4:
            rsp0 = Response[widx[0][0,1,:],widx[0][0,0,:],:,:,:]*widx[1][0,:][:, None, None, None]
            rsp1 = Response[widx[0][1,1,:],widx[0][1,0,:],:,:,:]*widx[1][1,:][:, None, None, None]
            rsp2 = Response[widx[0][2,1,:],widx[0][2,0,:],:,:,:]*widx[1][2,:][:, None, None, None]
            rsp3 = Response[widx[0][3,1,:],widx[0][3,0,:],:,:,:]*widx[1][3,:][:, None, None, None]
        else:
            print('Something went wrong; this should not happen.')
        
        rsp_mean = rsp0 + rsp1 + rsp2 + rsp3

    # return response
    return rsp_mean

         


def get_response_weights_vector(zenith,azimuth,binsize=5,cut=57.4):
    """
    Get Compton response pixel weights (four nearest neighbours),
    weighted by angular distance to zenith/azimuth vector(!) input.
    Binsize determines regular(!!!) sky coordinate grid in degrees.

    For single zenith/azimuth pairs use get_response_weights()
    
    :param: zenith        Zenith positions of the source with respect to the instrument (in deg)
    :param: azimuth       Azimuth positions of the source with respect to the instrument (in deg)
    :option: binsize      Default 5 deg (matching the sky dimension of the response). If set
                          differently, make sure it matches the sky dimension as otherwise,
                          false results may be returned
    :option: cut          Threshold to cut the response calculation after a certain zenith angle.
                          Default 57.4 deg (0.1 deg before last pixel reaching beyon 60 deg)
    """

    # assuming useful input:
    # azimuthal angle is periodic in the range [0,360[
    # zenith ranges from [0,180[ 
    # checking azimuth range (can be exactly 360?)
    azimuth[azimuth == 360] -= 0.01
    
    # check which pixel (index) was hit on regular grid
    hit_pixel_zi = np.floor(zenith/binsize).astype(int)
    hit_pixel_ai = np.floor(azimuth/binsize).astype(int)

    # and which pixel centre
    hit_pixel_z = (hit_pixel_zi+0.5)*binsize
    hit_pixel_a = (hit_pixel_ai+0.5)*binsize

    # check which zeniths are beyond threshold
    bad_idx = np.where(hit_pixel_z > cut) 
    
    # calculate nearest neighbour pixels indices
    za_idx = np.array([[np.floor(azimuth/binsize+0.5),np.floor(zenith/binsize+0.5)],
                       [np.floor(azimuth/binsize+0.5),np.floor(zenith/binsize-0.5)],
                       [np.floor(azimuth/binsize-0.5),np.floor(zenith/binsize+0.5)],
                       [np.floor(azimuth/binsize-0.5),np.floor(zenith/binsize-0.5)]]).astype(int)

    # take care of bounds at zenith (azimuth is allowed to be -1!)
    (za_idx[:,1,:])[np.where(za_idx[:,1,:] < 0)] += 1
    (za_idx[:,1,:])[np.where(za_idx[:,1,:] >= 180/binsize)] = int(180/binsize-1)
    # but azimuth may not be larger than range [0,360/binsize[
    (za_idx[:,0,:])[np.where(za_idx[:,0,:] >= 360/binsize)] = 0
    
    # and pixel centres of neighbours
    azimuth_neighbours = (za_idx[:,0]+0.5)*binsize
    zenith_neighbours = (za_idx[:,1]+0.5)*binsize

    # calculate angular distances to neighbours
    dists = angular_distance(azimuth_neighbours,zenith_neighbours,azimuth,zenith)

    # inverse weighting to get impact of neighbouring pixels
    n_in = len(zenith)
    weights = (1/dists)/np.sum(1/dists,axis=0).repeat(4).reshape(n_in,4).T
    # if pixel is hit directly, set weight to 1.0
    weights[np.isnan(weights)] = 1
    # set beyond threshold weights to zero
    weights[:,bad_idx] = 0

    return za_idx,weights



                              
def zenazi(scx_l, scx_b, scy_l, scy_b, scz_l, scz_b, src_l, src_b):
    """
    # from spimodfit zenazi function (with rotated axes (optical axis for COSI = z)
    # calculate angular distance wrt optical axis in zenith (theta) and
    # azimuth (phi): (zenazi function)
    # input: spacecraft pointing directions sc(xyz)_l/b; source coordinates src_l/b
    # output: source coordinates in spacecraft system frame
    
    Calculate zenith and azimuth angle of a point (a source) given the orientations
    of an instrument (or similar) in a certain coordinate frame (e.g. galactic).
    Each point in galactic coordinates can be uniquely mapped into zenith/azimuth of
    an instrument/observer/..., by using three Great Circles in x/y/z and retrieving
    the correct angles
    
    :param: scx_l      longitude of x-direction/coordinate
    :param: scx_b      latitude of x-direction/coordinate
    :param: scy_l      longitude of y-direction/coordinate
    :param: scy_b      latitude of y-direction/coordinate
    :param: scz_l      longitude of z-direction/coordinate
    :param: scz_b      latitude of z-direction/coordinate
    :param: src_l      SOURCE longitude
    :param: src_b      SOURCE latitude
    
    Space craft coordinates can also be vectors for quick computation for arrays
    """
    # Zenith is the distance from the optical axis (here z)
    costheta = GreatCircle(scz_l,scz_b,src_l,src_b)                                                                        
    # Azimuth is the combination of the remaining two
    cosx = GreatCircle(scx_l,scx_b,src_l,src_b)
    cosy = GreatCircle(scy_l,scy_b,src_l,src_b)
    
    # theta = zenith
    theta = np.rad2deg(np.arccos(costheta))
    # phi = azimuth
    phi = np.rad2deg(np.arctan2(cosy,cosx))
    
    # make azimuth going from 0 to 360 deg
    if phi.size == 1:
        if (phi < 0):
            phi += 360
        phi[phi < 0] += 360
    
    return theta,phi   

                    

def find_CDS_indices(phi,psi,chi,phi_bins,fisbel_bins_lat,fisbel_bins_lon):

    fisbel_index = np.argmin(angular_distance(chi,
                                              psi,
                                              fisbel_bins_lon,
                                              fisbel_bins_lat))

    phi_index = np.argmin(np.abs(phi-phi_bins))

    return(phi_index,fisbel_index)


def find_grid_indices(zen,azi,zen_bins,azi_bins):

    zen_index = np.argmin(np.abs(zen-zen_bins))
    azi_index = np.argmin(np.abs(azi-azi_bins))

    return(zen_index,azi_index)
