import ROOT as M
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import matplotlib as mpl

from tqdm.autonotebook import tqdm
from IPython.display import HTML

import warnings
warnings.filterwarnings('ignore')

# Load MEGAlib into ROOT
M.gSystem.Load("$(MEGAlib)/lib/libMEGAlib.so")

# Initialize MEGAlib
G = M.MGlobal()
G.Initialize()


class COSIpy:
    
    def __init__(self, data_dir, filename):
        self.data_dir = data_dir
        self.filename = filename
        self.horizon = 60 # field of view
        self.dataset = dataset()


    def read_COSI_DataSet(self):
        """
        Reads in MEGAlib .tra (or .tra.gz) file (self.filename) in given directory (self.data_dir)
        Returns COSI data set as a dictionary of the form
        COSI_DataSet = {'Full filename':self.data_dir+'/'+self.filename,
                        'Energies':erg,
                        'TimeTags':tt,
                        'Xpointings':np.array([lonX,latX]).T,
                        'Ypointings':np.array([lonY,latY]).T,
                        'Zpointings':np.array([lonZ,latZ]).T,
                        'Phi':phi,
                        'Chi local':chi_loc,
                        'Psi local':psi_loc,
                        'Distance':dist,
                        'Chi galactic':chi_gal,
                        'Psi galactic':psi_gal}
        :param: data_dir    Path of directory in which data set is stored
        :param: filename    .tra or tra.gz file to read in
        """

        # read COSI data starts here

        # check if file exists
        Reader = M.MFileEventsTra()
        if Reader.Open(M.MString(self.data_dir+'/'+self.filename)) == False:
            print("Unable to open file " + self.data_dir+'/'+self.filename + ". Aborting!")

        else:

            # initialise empty lists to store to
            # total photon energy
            erg = []
            # Time tag in UNIX time
            tt = []
            # Event Type (0: CE; 4:PE; ...)
            et = []
            # latitude of X direction of spacecraft
            latX = []
            # lontitude of X direction of spacecraft
            lonX = []
            # latitude of Z direction of spacecraft
            latZ = []
            # longitude of Z direction of spacecraft
            lonZ = []
            # Compton scattering angle
            phi = []
            # measured data space angle chi (azimuth direction; 0..360 deg)
            chi_loc = []
            # measured data space angle psi (polar direction; 0..180 deg)
            psi_loc = []
            # First lever arm distance in cm
            dist = []
            # measured gal angle chi (lon direction)
            chi_gal = []
            # measured gal angle psi (lat direction)
            psi_gal = [] 

            # browse through .tra file, select events, and sort into corresponding list
            while True:
                # the Reader class from MEGAlib knows where an event starts and ends and
                # returns the Event object which includes all information of an event
                Event = Reader.GetNextEvent()
                if not Event:
                    break
                # here only select Compton events (will add Photo events later as optional)
                
                # all calculations and definitions taken from:
                # /MEGAlib/src/response/src/MResponseImagingBinnedMode.cxx
                    
                # Total Energy
                erg.append(Event.Ei())
                # Time tag in UNIX seconds
                tt.append(Event.GetTime().GetAsSeconds())
                # Event type (0 = Compton, 4 = Photo)
                et.append(Event.GetEventType())
                # x axis of space craft pointing at GAL latitude
                latX.append(Event.GetGalacticPointingXAxisLatitude())
                # x axis of space craft pointing at GAL longitude
                lonX.append(Event.GetGalacticPointingXAxisLongitude())
                # z axis of space craft pointing at GAL latitude
                latZ.append(Event.GetGalacticPointingZAxisLatitude())
                # z axis of space craft pointing at GAL longitude
                lonZ.append(Event.GetGalacticPointingZAxisLongitude()) 
                
                # note that the y axis can be calculated from the X and Z components
                # therefore it is not saved, and will be computed further down
                    
                if Event.GetEventType() == M.MPhysicalEvent.c_Compton:    
                    # Compton scattering angle
                    phi.append(Event.Phi()) 
                    # data space angle chi (azimuth)
                    chi_loc.append((-Event.Dg()).Phi())
                    # data space angle psi (polar)
                    psi_loc.append((-Event.Dg()).Theta())
                    # interaction length between first and second scatter in cm
                    dist.append(Event.FirstLeverArm())
                    # gal longitude angle corresponding to chi
                    chi_gal.append((Event.GetGalacticPointingRotationMatrix()*Event.Dg()).Phi())
                    # gal longitude angle corresponding to chi
                    psi_gal.append((Event.GetGalacticPointingRotationMatrix()*Event.Dg()).Theta())
                    
            # because everything is built upon numpy arrays later, we will initialise them here
            erg = np.array(erg)
            tt = np.array(tt)
            et = np.array(et)
                
            latX = np.array(latX)
            lonX = np.array(lonX)
            # change longitudes to from 0..360 deg to -180..180 deg
            lonX[lonX > np.pi] -= 2*np.pi
            
            latZ = np.array(latZ)
            lonZ = np.array(lonZ)
            # change longitudes to from 0..360 deg to -180..180 deg
            lonZ[lonZ > np.pi] -= 2*np.pi
            
            phi = np.array(phi)
            
            chi_loc = np.array(chi_loc)
            # change azimuth angle to 0..360 deg
            chi_loc[chi_loc < 0] += 2*np.pi
            
            psi_loc = np.array(psi_loc)
            
            dist = np.array(dist)
            
            chi_gal = np.array(chi_gal)
            psi_gal = np.array(psi_gal)
            
            # construct Y direction from X and Z direction
            lonlatY = construct_scy(np.rad2deg(lonX),np.rad2deg(latX),
                                    np.rad2deg(lonZ),np.rad2deg(latZ))
            lonY = np.deg2rad(lonlatY[0])
            latY = np.deg2rad(lonlatY[1])
            
            # avoid negative zeros
            chi_loc[np.where(chi_loc == 0.0)] = np.abs(chi_loc[np.where(chi_loc == 0.0)])
            
            # make observation dictionary
            COSI_DataSet = {'Full filename':self.data_dir+'/'+self.filename,
                            'Energies':erg,
                            'TimeTags':tt,
                            'Xpointings':np.array([lonX,latX]).T,
                            'Ypointings':np.array([lonY,latY]).T,
                            'Zpointings':np.array([lonZ,latZ]).T,
                            'Phi':phi,
                            'Chi local':chi_loc,
                            'Psi local':psi_loc,
                            'Distance':dist,
                            'Chi galactic':chi_gal,
                            'Psi galactic':psi_gal}
            
            self.dataset.data = COSI_DataSet

        

    def plot_raw_spectrum(self,bins=100):
        """
        Plot spectrum of the raw counts in the given data set.
        Returns energy bin information as well as counts per bin.
        
        :param: bins   number of bins to be used for plotting, or specific defintion of energy bin edges
                       Default = 100 equally spaced bins
        """
        
        try:
            
            # histogramming the energies in the data set
            self.dataset.energy_spec_data, self.dataset.energy_bin_edges = np.histogram(self.dataset.data['Energies'],bins=bins)
            # bin minima and maxima
            self.dataset.energy_bin_min = self.dataset.energy_bin_edges[0:-1]
            self.dataset.energy_bin_max = self.dataset.energy_bin_edges[1:]
            # center of bins
            self.dataset.energy_bin_cen = 0.5*(self.dataset.energy_bin_max+self.dataset.energy_bin_min)
            # half-width of bins
            self.dataset.energy_bin_wid = 0.5*(self.dataset.energy_bin_max-self.dataset.energy_bin_min)
            # number of energy bin
            self.dataset.n_energy_bins  = len(self.dataset.energy_bin_cen)
            
            # make a plot
            fig, ax = plt.subplots(1,1,figsize=(8,6))
            #ax.plot(self.dataset.energy_bin_cen,self.dataset.energy_spec_data,where='mid')
            for i in range(self.dataset.n_energy_bins):
                if i == 0 :
                    ax.plot([self.dataset.energy_bin_min[i],
                             self.dataset.energy_bin_min[i]],
                            [0,
                             self.dataset.energy_spec_data[i]/self.dataset.energy_bin_wid[i]/2],marker='',linestyle='-',color='black')
                elif i < self.dataset.n_energy_bins:
                    ax.plot([self.dataset.energy_bin_min[i],
                             self.dataset.energy_bin_min[i]],
                            [self.dataset.energy_spec_data[i-1]/self.dataset.energy_bin_wid[i-1]/2,
                             self.dataset.energy_spec_data[i]/self.dataset.energy_bin_wid[i]/2],marker='',linestyle='-',color='black')
    
                ax.plot([self.dataset.energy_bin_max[i],
                         self.dataset.energy_bin_max[i]],
                        [self.dataset.energy_spec_data[i]/self.dataset.energy_bin_wid[i]/2,
                         0],marker='',linestyle='-',color='black')
                ax.plot([self.dataset.energy_bin_min[i],
                         self.dataset.energy_bin_max[i]],
                        [self.dataset.energy_spec_data[i]/self.dataset.energy_bin_wid[i]/2,
                         self.dataset.energy_spec_data[i]/self.dataset.energy_bin_wid[i]/2],marker='',linestyle='-',color='black')
            ax.set_xlabel('Energy [keV]')
            ax.set_ylabel('Counts [cnts/keV]')


        except AttributeError:
            print('not working here ...')


        
    def plot_elevation(self,l_src,b_src,name_src):
        """
        Plot the elevation of a list of sources for the given data set.
        Returns the zenith pointings in longituda and latitude as well as the respective time tags

        :param: l_src      list of longitude coordinates for sources to calculate elevation
        :param: b_src      list of latitude coordinates for sources to calculate elevation
        :param: name_src   name of respective sources
        
        """

        try:
        
            # extracting the pointing direction information (zenith direction for each time tag)
            self.dataset.l_pointing, self.dataset.b_pointing = np.rad2deg(self.dataset.data['Zpointings'][:,0]), np.rad2deg(self.dataset.data['Zpointings'][:,1])
        
            # extracting time tags for corresponding pointings
            self.dataset.t_pointing = self.dataset.data['TimeTags']

            # make a plot for
            try:
            
                # set up figure and axis
                fig, ax = plt.subplots(1,1,figsize=(16,6))
            
                # loop over sources
                for l,b,name in zip(l_src,b_src,name_src):
                    # the elevation is just the edge of the field of view (= horizon) minus the angular distance to the source
                    tmp_elevation = self.horizon-angular_distance(l,b,self.dataset.l_pointing,self.dataset.b_pointing)
                    ax.plot(self.dataset.t_pointing,tmp_elevation,label=name,marker='.',color='k',linestyle='')
                
                ax.set_xlabel('Unix time [s]')
                ax.set_ylabel('Elevation above COSI horizon [deg]')
            
                # indicating what the maximum elevation is (TS: the y-axis might be confusing)
                ax.axhline(self.horizon,linestyle='--',color='black')
                ax.text(np.median(minmax(self.dataset.t_pointing)),self.horizon,
                        'Zenith',horizontalalignment='center',verticalalignment='bottom')
                ax.set_ylim(0,)
                ax.legend()
            
            except TypeError:
                print('Input longitudes, latitudes, and names should be lists.')

        except AttributeError:
            print('not working here ...')
        




    def plot_lightcurve(self):
        try:
            try:
                self.dataset.light_curve = [len(self.dataset.data_time_tagged[i]['Indices'])/self.dataset.data_time_tagged[i]['DeltaTime'] for i in range(self.dataset.n_time_bins)]
                self.dataset.times_bins  = [0] + [self.dataset.data_time_tagged[i]['DeltaTime'] for i in range(self.dataset.n_time_bins)]
                self.dataset.times_edges = np.cumsum(self.dataset.times_bins)
                self.dataset.times_min   = self.dataset.times_edges[0:-1]
                self.dataset.times_max   = self.dataset.times_edges[1:]
                self.dataset.times_cen   = 0.5*(self.dataset.times_max+self.dataset.times_min)
                self.dataset.times_wid   = 0.5*(self.dataset.times_max-self.dataset.times_min)
                
                fig, ax = plt.subplots(1,1,figsize=(8,6))
                ax.step(self.dataset.times_cen,self.dataset.light_curve,where='mid')
                ax.set_xlabel('Seconds since UNIX second '+str('%.2f' % self.dataset.data['TimeTags'][0])+' [s]')
                ax.set_ylabel('Count rate [cnts/s]')

            except:
                print('Time tags not binned, yet, binning for 1 hour intervals now ...')
                self.dataset.time_binning_tags()
                self.plot_lightcurve()

        except AttributeError:
            print('not working here ...')


  
    @property
    def available_attributes(self):
        """
        Print available attributes of current instance of COSIpy analysis object.
        """
        for key in self.__dict__.keys():
            print(key)

    @property
    def available_methods(self):
        """
        Print available methods (and attributes) of current instance.
        """
        print([thing for thing in dir(self) if not thing.startswith('__') and not np.any(np.isin(list(self.__dict__.keys()),thing))])

        

    def init_spi(self):
        self.spi = SPI()


class dataset(COSIpy):


    def __init__(self,name='name'):
        #super().__init__(*args, **kwargs)
        self.name = name
        self.data = None
        """self.data = {'Full filename':self.data_dir+'/'+self.filename,
                     'Energies':None,
                     'TimeTags':None,
                     'Xpointings':None,
                     'Ypointings':None,
                     'Zpointings':None,
                     'Phi':None,
                     'Chi local':None,
                     'Psi local':None,
                     'Distance':None,
                     'Chi galactic':None,
                     'Psi galactic':None}"""


    @property
    def plot_( self ):
        raise AttributeError( "'Bar' object has no attribute 'foo'" )
    
        
    @property
    def data_info(self):
        """
        Prints tabularised information of the individual COSI_DataSet dictionary entries.
        In particular mean, std, min, and max.
        TS: might be using pandas data frame in the future because it's doing that automatically.
        """

        # only if the data is already read in
        if isinstance(self.data,dict):
            # formatting the table
            row_format ='{:>15}' * 5

            # which file is used
            print('Full filename'+':',self.data['Full filename'])

            # first table row
            print(row_format.format('',*['mean','std','min','max']))

            # loop over remaining dictionary entries and calculate and print content
            for key in self.data.keys():
                # ignore strings (first entry)
                if isinstance(self.data[key],str):
                    pass
                else:
                    print(row_format.format(key+':',
                                            str('%.3f' % np.mean(self.data[key])),
                                            str('%.3f' % np.std(self.data[key])),
                                            str('%.3f' % np.min(self.data[key])),
                                            str('%.3f' % np.max(self.data[key]))))

        else:
            print('Data set not read in, use read_COSI_DataSet() first.')



    def time_binning_tags(self,time_bin_size=3600,time_binning=None):
        """
        Get COSI data reformatted to a data set per time bin.
        Output is a dictionary of the form:
        data = {'Full filename':COSI_Data['Full filename'],
                'Bin number':b,
                'Indices':tdx,
                'DeltaTime':time}
        :param: COSI_Data       dictionary from read_COSI_DataSet
        :param: time_bin_size   equi-sized time intervals to bin select photons in data set
                                Default: 3600 (in units of seconds)
        """

        self.init_time_bin_size = time_bin_size

        self.n_time_bins = int(np.ceil(np.diff(minmax(self.data['TimeTags']))/self.init_time_bin_size))

        s2b = 1./self.init_time_bin_size

        self.last_bin_size = np.diff(minmax(self.data['TimeTags']))[0]-(self.n_time_bins-1)/s2b

        self.data_time_tagged = []

        for b in range(self.n_time_bins):
            tdx = np.where( (minmin(self.data['TimeTags'])*s2b >= b) &
                            (minmin(self.data['TimeTags'])*s2b <  b+1) )[0]

            tmp_data = {'Full filename':self.data['Full filename'],
                        'Bin number':b,
                        'Indices':tdx,
                        'DeltaTime':self.init_time_bin_size if b < self.n_time_bins-1 else self.last_bin_size}

            self.data_time_tagged.append(tmp_data)

        self.times = TIME()
        self.times.times_bins  = [0] + [self.data_time_tagged[i]['DeltaTime'] for i in range(self.n_time_bins)]
        self.times.n_time_bins = self.n_time_bins
        self.times.times_edges = np.cumsum(self.times.times_bins)
        self.times.times_min   = self.times.times_edges[0:-1]
        self.times.times_max   = self.times.times_edges[1:]
        self.times.times_cen   = 0.5*(self.times.times_max+self.times.times_min)
        self.times.times_wid   = 0.5*(self.times.times_max-self.times.times_min)

        self.times.n_ph_t = np.array([len(self.data_time_tagged[i]['Indices']) for i in range(self.times.n_time_bins)])
        self.times.n_ph_dx = np.where(self.times.n_ph_t != 0)[0]
        self.times.n_ph   = len(self.times.n_ph_dx)
        self.times.total_time  = 2*np.sum(self.times.times_wid[self.times.n_ph_dx])




    ##########################################
    ## JB May 4, 2022: extend time_binning_tags function to take
    # arbitrary time bins rather than one time bin size
    def time_binning_tags_2(self,time_bin_size,time_binning=None):
        """
        Get COSI data reformatted to a data set per time bin.
        Output is a dictionary of the form:
        data = {'Full filename':COSI_Data['Full filename'],
                'Bin number':b,
                'Indices':tdx,
                'DeltaTime':time}
        :param: COSI_Data       dictionary from read_COSI_DataSet
        :param: time_bin_size   np.array of custom time bin edges
                                Default: none
        """
        # general time bin size
        # TS: January 25: can be set to one number or general shape
        #if not isinstance(time_bin_size, (list, tuple, np.ndarray)):
        self.init_time_bin_size = time_bin_size
        self.n_time_bins = len(self.init_time_bin_size)-1
        # time conversions seconds to bin-size
        s2b = 1./np.diff(self.init_time_bin_size)
        # calculate last time bin interval
        # ( how much is left in the data set inside the last bin )
        self.last_bin_size = np.diff(self.init_time_bin_size)[-1] #self.init_time_bin_size[-1]
        # fill data (as formatted per time tagged of individual time bins of size time_bin_size)
        self.data_time_tagged = []
        for b in range(self.n_time_bins):
            tdx = np.where( (minmin(self.data['TimeTags'])*s2b[b] >= b) &
                            (minmin(self.data['TimeTags'])*s2b[b] <  b+1) )[0]
            tmp_data = {'Full filename':self.data['Full filename'],
                        'Bin number':b,
                        'Indices':tdx,
                        'DeltaTime':np.diff(self.init_time_bin_size)[b] if b < self.n_time_bins else self.last_bin_size}
            self.data_time_tagged.append(tmp_data)
        self.times = TIME()
        self.times.times_bins  = [self.data_time_tagged[i]['DeltaTime'] for i in range(self.n_time_bins)]#-1)]
        self.times.n_time_bins = self.n_time_bins
        self.times.times_edges = time_bin_size #np.cumsum(self.times.times_bins)
        self.times.times_min   = self.times.times_edges[0:-1]
        self.times.times_max   = self.times.times_edges[1:]
        self.times.times_cen   = 0.5*(self.times.times_max+self.times.times_min)
        self.times.times_wid   = 0.5*(self.times.times_max-self.times.times_min)
        # need to take into account empty time bins
        self.times.n_ph_t = np.array([len(self.data_time_tagged[i]['Indices']) for i in range(self.times.n_time_bins)])
        #self.times.n_ph_t = np.array([len(self.data_time_tagged[i]['Indices']) for i in range(self.times.n_time_bins-1)])
        self.times.n_ph_dx = np.where(self.times.n_ph_t != 0)[0]
        self.times.n_ph   = len(self.times.n_ph_dx)
        self.times.total_time  = 2*np.sum(self.times.times_wid)
  
            

    def init_binning(self,energy_bin_edges=np.linspace(100,1000,10),pixel_size=5.):
        """
        Prepare data set to input energy, time, and angle information in a binned arrary
        """

        # energy
        self.energies = ENERGY(bin_edges=energy_bin_edges)

        # Comptel data space
        # pixel definitions
        self.pixel_size = pixel_size # (in deg)
        # pixel area in rad^2
        self.pixel_area = np.deg2rad(self.pixel_size)**2
        # number of fisbel bins
        self.n_fisbel_bins = int(np.floor(4*np.pi/self.pixel_area))
        # number of phi (Compton Scatter Angle) bins
        self.n_phi_bins = int(180./self.pixel_size)
        
        # psi, chi
        self.fisbels = FISBEL(n_bins=self.n_fisbel_bins)
        # phi (Compton scattering angle)
        self.phis = CSA(n_bins=self.n_phi_bins)




    def get_binned_data(self):
        """
        Bin data according to definitions in init_binning()
        """

        try:
    
            # init data array
            self.binned_data = np.zeros((self.times.n_ph,#self.times.n_time_bins,
                                         self.energies.n_energy_bins,
                                         self.phis.n_phi_bins,
                                         self.fisbels.n_fisbel_bins))

            ph_dx = 0

            for t in range(self.times.n_time_bins):

                if self.times.n_ph_t[t] != 0:
                
                    # use indexed photons for specific time interval
                    idx_tmp = self.data_time_tagged[t]['Indices']

                    # temporary CDS indexed for time interval
                    phi_tmp = self.data['Phi'][idx_tmp]
                    psi_tmp = self.data['Psi local'][idx_tmp]
                    chi_tmp = self.data['Chi local'][idx_tmp]
                    # and energies for indexed time interval
                    erg_tmp = self.data['Energies'][idx_tmp]
                    
                    # because FISBEL is not monotonic in both dimensions, loop over FISBEL bins
                    for f in range(self.n_fisbel_bins):

                        # select 2D pixel range where photons fall in
                        fisbel_idx_tmp = np.where((chi_tmp >= self.fisbels.lon_min[f]) &
                                                  (chi_tmp <  self.fisbels.lon_max[f]) &
                                                  (psi_tmp >= self.fisbels.lat_min[f]) &
                                                  (psi_tmp <  self.fisbels.lat_max[f]))[0]
                        
                        # multi-D histogramming of the events in energy and phi
                        tmp_hist = np.histogramdd(np.array([erg_tmp[fisbel_idx_tmp],
                                                            phi_tmp[fisbel_idx_tmp]]).T,
                                                  bins=np.array([self.energies.energy_bin_edges,
                                                                 self.phis.phi_edges]))#,

                        # fill into binned_data array
                        self.binned_data[ph_dx,:,:,f] = tmp_hist[0]
                    ph_dx += 1

                    
        except AttributeError:
        
            print('Something is not defined; run init_binning() first')




    def bin_for_angles(self,binned_array=None):
        """
        extract phi, psi, and chi information per time and energy
        """

        # temporary angle bin definition
        # chi
        ll = self.fisbels.lon_cen
        dll = self.fisbels.lon_wid

        # psi
        bb = self.fisbels.lat_cen
        dbb = self.fisbels.lat_wid

        # phi
        pp = self.phis.phi_cen
        dpp = self.phis.phi_wid
        n_pp = len(pp)
        
        # find indices for psi and chi
        uniq_bb = np.unique(bb)
        n_bb = len(uniq_bb)
        bb_idx = []
        for i in range(n_bb):
            bb_idx.append(np.where(bb == uniq_bb[i])[0])

        uniq_ll = np.unique(ll)
        n_ll = len(uniq_ll)
        ll_idx = []
        for i in range(n_ll):
            ll_idx.append(np.where(ll == uniq_ll[i])[0])


        if np.any(binned_array == None):
        
            self.phi_binned = np.zeros((self.times.n_time_bins,
                                        self.energies.n_energy_bins,
                                        n_pp))

            for t in range(self.times.n_time_bins):
                for e in range(self.energies.n_energy_bins):
                    for i in range(n_pp):
                        self.phi_binned[t,e,i] = np.sum(self.binned_data[t,e,i,:])
                        
                    
            self.psi_binned = np.zeros((self.times.n_time_bins,
                                        self.energies.n_energy_bins,
                                        n_bb))

            for t in range(self.times.n_time_bins):
                for e in range(self.energies.n_energy_bins):
                    for i in range(n_bb):
                        self.psi_binned[t,e,i] = np.sum(self.binned_data[t,e,:,bb_idx[i]])

                    
            self.chi_binned = np.zeros((self.times.n_time_bins,
                                        self.energies.n_energy_bins,
                                        n_ll))

            for t in range(self.times.n_time_bins):
                for e in range(self.energies.n_energy_bins):
                    for i in range(n_ll):
                        self.chi_binned[t,e,i] = np.sum(self.binned_data[t,e,:,ll_idx[i]])


        else:

            phi_binned = np.zeros((self.times.n_time_bins,
                                   self.energies.n_energy_bins,
                                   n_pp))

            for t in range(self.times.n_time_bins):
                for e in range(self.energies.n_energy_bins):
                    for i in range(n_pp):
                        phi_binned[t,e,i] = np.sum(binned_array[t,e,i,:])


            psi_binned = np.zeros((self.times.n_time_bins,
                                   self.energies.n_energy_bins,
                                   n_bb))

            for t in range(self.times.n_time_bins):
                for e in range(self.energies.n_energy_bins):
                    for i in range(n_bb):
                        psi_binned[t,e,i] = np.sum(binned_array[t,e,:,bb_idx[i]])


            chi_binned = np.zeros((self.times.n_time_bins,
                                   self.energies.n_energy_bins,
                                   n_ll))

            for t in range(self.times.n_time_bins):
                for e in range(self.energies.n_energy_bins):
                    for i in range(n_ll):
                        chi_binned[t,e,i] = np.sum(binned_array[t,e,:,ll_idx[i]])


            return phi_binned, psi_binned, chi_binned


                        
            
    def plot_raw_spectrum(self,mode='total'):
        """
        Plot spectrum of the raw counts in the given (binned) data set.
        """
        
        try:

            #if (self.binned_data.shape[0] != self.times.n_time_bins):
            if (self.binned_data.shape[0] != self.times.n_ph):
                print('The binning was updated since the last call! Run time_binning_tags() for your current binning.')
            else:
                
                if mode == 'total':
                    self.energies.energy_spec_data = np.sum(self.binned_data,axis=(0,2,3)).reshape(1,self.energies.n_energy_bins)
                    #self.times.total_time = np.sum(self.times.times_wid*2)
                    self.energies.energy_spec_data /= self.times.total_time
                elif mode == 'all':
                    self.energies.energy_spec_data = np.sum(self.binned_data,axis=(2,3))
                    self.energies.energy_spec_data /= (self.times.times_wid[:,None]*2)
                else:
                    print('Plotting spectrum for time bin number '+str(mode))
                    self.energies.energy_spec_data = np.sum(self.binned_data,axis=(2,3))[mode,:].reshape(1,self.energies.n_energy_bins)
                    self.energies.energy_spec_data /= (self.times.times_wid[mode]*2)

                self.energies.energy_spec_data /= self.energies.energy_bin_wid
                self.energies.energy_spec_data /= 2 # 2 for half bin widths
    
                # make a plot
                fig, ax = plt.subplots(1,1,figsize=(8,6))
                for t in range(self.energies.energy_spec_data.shape[0]):
                    for i in range(self.energies.n_energy_bins):
                        if i == 0 :
                            ax.plot([self.energies.energy_bin_min[i],
                                     self.energies.energy_bin_min[i]],
                                    [0,
                                     self.energies.energy_spec_data[t,i]],marker='',linestyle='-',color='black')
                        elif i < self.energies.n_energy_bins:
                            ax.plot([self.energies.energy_bin_min[i],
                                     self.energies.energy_bin_min[i]],
                                    [self.energies.energy_spec_data[t,i-1],
                                     self.energies.energy_spec_data[t,i]],marker='',linestyle='-',color='black')
    
                        ax.plot([self.energies.energy_bin_max[i],
                                 self.energies.energy_bin_max[i]],
                                [self.energies.energy_spec_data[t,i],
                                0],marker='',linestyle='-',color='black')
                        ax.plot([self.energies.energy_bin_min[i],
                                 self.energies.energy_bin_max[i]],
                                [self.energies.energy_spec_data[t,i],
                                 self.energies.energy_spec_data[t,i]],marker='',linestyle='-',color='black')
                ax.set_xlabel('Energy [keV]')
                ax.set_ylabel('Counts [cnts/s/keV]')
            
                if (mode == 'total') | (mode == 'all'):
                    text_time = self.times.total_time
                else:
                    text_time = self.times.times_wid[mode]*2
                
                ax.text(0.7,0.8,
                        'Time bin: '+str(mode)+
                        '\n ('+str('%.1f' % text_time)+' s)',
                        transform=ax.transAxes,fontsize=16)


        except AttributeError:
            print('Data not yet binned? Use get_binned_data() first.')




    def plot_lightcurve(self,mode='total'):
        """
        Plotting the light curve of a binned data set, either for the total spectrum,  all energy bins, or individually
        :option: energy   'total': plotting sum over angles, and energies (default)
                          'all':   plotting light curve of all energy bins on top of each other
                          int:     number of energy bin for which the light curve should be plotted
        """
        
        try:

            if (self.binned_data.shape[0] != self.times.n_ph):
                print('The binning was updated since the last call! Run time_binning_tags() for your current binning.')
            else:

                # sum over angles to get only time and energy array
                self.lc_all = np.sum(self.binned_data,axis=(2,3))
                #self.lc_all = self.lc_all/self.times.total_time
                self.lc_all = self.lc_all/self.init_time_bin_size # normalize by bin time
                
                if mode == 'total':
                    self.light_curve = np.sum(self.lc_all,axis=1)
                elif mode == 'all':
                    self.light_curve = np.copy(self.lc_all)
                else:
                    self.light_curve = self.lc_all[:,mode]
                
                fig, ax = plt.subplots(1,1,figsize=(8,6))
                ax.step(self.times.times_cen[self.times.n_ph_dx],self.light_curve,where='mid',color="black")
                ax.set_xlabel('Seconds since UNIX second '+str('%.2f' % self.data['TimeTags'][0])+' [s]')
                ax.set_ylabel('Count rate [cnts/s]')

                if (mode == 'total') | (mode == 'all'):
                    text_energy = minmax(self.energies.energy_bin_edges)
                else:
                    text_energy = np.array([self.energies.energy_bin_min[mode],
                                            self.energies.energy_bin_max[mode]])

                ax.text(0.6,0.8,
                        'Energy bin: '+str(mode)+
                        '\n ('+str('%.1f-%.1f keV' % (text_energy[0],text_energy[1]))+')',
                        transform=ax.transAxes,fontsize=16)
            
        except AttributeError:
            print('Data not yet binned? Use get_binned_data() first.')
            



class TIME(dataset):

    def __init__(self):
        pass



class ENERGY(dataset):

    def __init__(self,bin_edges=np.linspace(100,1000,10)):
        self.energy_bin_edges = bin_edges
        self.energy_bin_min = self.energy_bin_edges[0:-1]
        self.energy_bin_max = self.energy_bin_edges[1:]
        # center of bins
        self.energy_bin_cen = 0.5*(self.energy_bin_max+self.energy_bin_min)
        # half-width of bins
        self.energy_bin_wid = 0.5*(self.energy_bin_max-self.energy_bin_min)
        # number of energy bin
        self.n_energy_bins  = len(self.energy_bin_cen)
    

class CSA(dataset):

    def __init__(self,n_bins=36):
        self.n_phi_bins = n_bins
        
        self.phi_edges = np.linspace(0,np.pi,n_bins+1)

        self.phi_min = self.phi_edges[0:-1]
        self.phi_max = self.phi_edges[1:]
        
        self.phi_cen = 0.5*(self.phi_max+self.phi_min)

        self.phi_wid = 0.5*(self.phi_max-self.phi_min)


class FISBEL(dataset):

    # tolerance of binning (strange effects when binning in angle space)
    tol = 6 #1e-6
    
    def __init__(self,n_bins=1650):
        self.fisbelbins = self.get_FISBEL_binning(n_bins=n_bins)
        self.n_fisbel_bins = n_bins
    
        self.lat_cen = (np.array(self.fisbelbins[0]))[:,0]#-np.pi/2.
        self.lon_cen = (np.array(self.fisbelbins[0]))[:,1]#-np.pi

        self.lat_wid = (np.array(self.fisbelbins[1]))[:,0].reshape(n_bins,)
        self.lon_wid = (np.array(self.fisbelbins[1]))[:,1].reshape(n_bins,)

        self.lat_min = np.around(self.lat_cen-self.lat_wid/2,self.tol)
        self.lon_min = np.around(self.lon_cen-self.lon_wid/2,self.tol)

        self.lat_max = np.around(self.lat_cen+self.lat_wid/2,self.tol)
        self.lon_max = np.around(self.lon_cen+self.lon_wid/2,self.tol)

        self.lat_edges = np.concatenate([np.array([0.]),self.lat_max])
        self.lon_edges = np.concatenate([np.array([0.]),self.lon_max])
        
    
    @staticmethod
    def get_FISBEL_binning(n_bins=1650,lon_shift=0,verbose=False):
        """
        MEGAlib's FISBEL spherical axis binning
        /MEGAlib/src/global/misc/src/MBinnerFISBEL.cxx
        Used to make images with more information (pixels) in the centre, and less at higher latitudes
        // CASEBALL - Constant area, squares at equator, borders along latitude & longitude
        
        // Rules:
        // (1) constant area
        // (2) squarish at equator
        // (3) All borders along latitude & longitude lines
        // (4) latitude distance close to equal
        
        Don't know where "lon_shift" is actually required, keeping it for the moment
        :param: n_bins      number of bins to populate the sky (typically 1650 (~5 deg) or 4583 (~3 deg))
        :param: lon_shift   ???
        :option: verbose    True = more verbose output
        
        Returns CoordinatePairs and respective Binsizes
        """
        
        # if only one bin, full sky in one bin
        if (n_bins == 1):
            LatitudeBinEdges = [0,np.pi]
            LongitudeBins = [1]
            n_collars = 1
        else:
            FixBinArea = 4*np.pi/n_bins
            SquareLength = np.sqrt(FixBinArea)
            
            n_collars = np.int((np.pi/SquareLength-1)+0.5) + 2
            # -1 for half one top AND Bottom, 0.5 to round to next integer
            # +2 for the half on top and bottom
            
            verb(verbose,'Number of bins: %4i' % n_bins)
            verb(verbose,'Fix bin area: %6.3f' % FixBinArea)
            verb(verbose,'Square length: %6.3f' % SquareLength)
            verb(verbose,'Number of collars: %4i' % n_collars)
            
            LongitudeBins = np.zeros(n_collars)
            LatitudeBinEdges = np.zeros(n_collars+1)
            
            # Top and bottom first
            LatitudeBinEdges[0] = 0
            LatitudeBinEdges[n_collars] = np.pi
            
            # Start on top and bottom with a circular region:
            LongitudeBins[0] = 1
            LongitudeBins[n_collars - 1] = LongitudeBins[0]
            
            LatitudeBinEdges[1] = np.arccos(1 - 2.0/n_bins)
            LatitudeBinEdges[n_collars - 1] = np.pi - LatitudeBinEdges[1]
            
            # now iterate over remaining bins
            for collar in range(1,np.int(np.ceil(n_collars/2))):
                UnusedLatitude = LatitudeBinEdges[n_collars-collar] - LatitudeBinEdges[collar]
                UnusedCollars = n_collars - 2*collar
                
                NextEdgeEstimate = LatitudeBinEdges[collar] + UnusedLatitude/UnusedCollars
                NextBinsEstimate = 2*np.pi * (np.cos(LatitudeBinEdges[collar]) - np.cos(NextEdgeEstimate)) / FixBinArea
                
                # roundgind
                NextBins = np.int(NextBinsEstimate+0.5)
                NextEdge = np.arccos(np.cos(LatitudeBinEdges[collar]) - NextBins*FixBinArea/2/np.pi)
            
                # insert at correct position
                LongitudeBins[collar] = NextBins
                if (collar != n_collars/2):
                    LatitudeBinEdges[collar+1] = NextEdge
                LongitudeBins[n_collars-collar-1] = NextBins
                if (collar != n_collars/2):
                    LatitudeBinEdges[n_collars-collar-1] = np.pi - NextEdge


        LongitudeBinEdges = []
        for nl in LongitudeBins:
            if (nl == 1):
                LongitudeBinEdges.append(np.array([0,2*np.pi]))
            else:
                n_lon_edges = int(nl+1)
                LongitudeBinEdges.append(np.linspace(0,2*np.pi,n_lon_edges))
                
        CoordinatePairs = []
        Binsizes = []
        for c in range(n_collars):
            for l in range(np.int(LongitudeBins[c])):
                CoordinatePairs.append([np.mean(LatitudeBinEdges[c:c+2]),np.mean(LongitudeBinEdges[c][l:l+2])])
                Binsizes.append([np.diff(LatitudeBinEdges[c:c+2]),np.diff(LongitudeBinEdges[c][l:l+2])])
     
        return CoordinatePairs,Binsizes



        
class Pointing():

    def __init__(self,
                 dataset=None,
                 angle_threshold=5.):
        self.angle_threshold = angle_threshold
        self.construct_pointings(dataset)


    def construct_pointings(self,dataset):
        """
        Construct a list of 'stable' pointings during a steadily moving COSI observation:
        Depending on the angular threshold, a number of triggers is collected into one time bin,
        making an average observation, Xpointing, Ypointing, ZPointing, during this time.
        This defines the exposure time (not dead-time corrected) for a particular pointing.
        TS: need a cleaner pointing definition from ori files (or auxiliary information)
        :param:  COSI_Data.dataset  COSI data set dictionary from read_COSI_DataSet()
        :option: angle_threshold    Default = 5 (deg): if COSI changes by more than angle_threshold, the accumulation
                                    of triggers starts over, defining a new pointing.
        
        Output: xpoins,ypoins,zpoins,dtpoins
        Pointings in x-, y-, z-direction, and observation time (interval) for each set of triggers.
        TS Warning: there could be something wrong, but i don't see where
        """
    
        self.xpoins  = []
        self.ypoins  = []
        self.zpoins  = []
        self.dtpoins = []

        self._n_ph          = []
        self._n_pointings   = []
        self._pointing_cuts = []
    
        for t in range(dataset.n_time_bins):

            # indices for this time bin
            idx_tmp = dataset.data_time_tagged[t]['Indices']
        
            # get number of photons as number of indexed triggers
            n_ph = len(idx_tmp)
        
            # temporary pointing coordinates of indexed events in deg
            zpoins_tmp = np.rad2deg(dataset.data['Zpointings'][idx_tmp])
            ypoins_tmp = np.rad2deg(dataset.data['Ypointings'][idx_tmp])
            xpoins_tmp = np.rad2deg(dataset.data['Xpointings'][idx_tmp])
        
            # time tages associated with indices
            times_tmp  = dataset.data['TimeTags'][idx_tmp]
        
            # sorting of times since, especially for simulations,
            # the time stamps may be arbitrary and not in sequence
            sort_idx = np.argsort(times_tmp)
            # sort times and x/y/z coordinates according to that
            times_tmp  = times_tmp[sort_idx]
            zpoins_tmp = zpoins_tmp[sort_idx,:]
            ypoins_tmp = ypoins_tmp[sort_idx,:]
            xpoins_tmp = xpoins_tmp[sort_idx,:]
        
            # cross-check number of photons
            n_ph2 = len(zpoins_tmp[:,0])
        
            if n_ph != n_ph2:
                print('something went wrong ...; this should never happen here!')
            else:
            
                # begin iterative combination of triggers up to angle_threshold
                # first photon

                # if no photon in time bin:
                if n_ph == 0:

                    save_zpoins = [[np.nan,np.nan]]
                    save_ypoins = [[np.nan,np.nan]]
                    save_xpoins = [[np.nan,np.nan]]
                    save_times = []

                else:
                
                    i = 0
                    save_zpoins = [zpoins_tmp[i,:]]
                    save_ypoins = [ypoins_tmp[i,:]]
                    save_xpoins = [xpoins_tmp[i,:]]
                    save_times = []
                    save_f = [0]
        
                    # next photon wrt to first trigger in block
                    for f in range(n_ph):
                
                        # calculate angular distance in all 3 directions
                        # (very important because non-circular response)
                        dz = angular_distance(zpoins_tmp[f,0],zpoins_tmp[f,1],
                                              zpoins_tmp[i,0],zpoins_tmp[i,1])
                        dy = angular_distance(ypoins_tmp[f,0],ypoins_tmp[f,1],
                                              ypoins_tmp[i,0],ypoins_tmp[i,1])
                        dx = angular_distance(xpoins_tmp[f,0],xpoins_tmp[f,1],
                                              xpoins_tmp[i,0],xpoins_tmp[i,1])
                
                        # time bin (see warning above, and work-around below)
                        dt = times_tmp[f]-times_tmp[i]
                        # if any of the angles is greater than threshold,
                        # save the first trigger as average pointing
                        # this is arbitrary and could be any point inside the block;
                        # it only defines a point in which the stability is
                        # angle_threshold degrees.
                        if ((dz >= self.angle_threshold) | (dy >= self.angle_threshold) | (dx >= self.angle_threshold)) | (f == n_ph-1):
                            save_zpoins.append(zpoins_tmp[f,:])
                            save_ypoins.append(ypoins_tmp[f,:])
                            save_xpoins.append(xpoins_tmp[f,:])
                            save_times.append(dt)
                            save_f.append(f)
                            # make last to first and repeat until all photons are inside
                            i = f

                    # append last time segment as delta between time bin edge and 
                    # that was already accounted for
                save_times.append(dataset.data_time_tagged[t]['DeltaTime']-np.sum(save_times)) 

                self.xpoins.extend(save_xpoins)
                self.ypoins.extend(save_ypoins)
                self.zpoins.extend(save_zpoins)
                self.dtpoins.extend(save_times)

                self._n_ph.append(n_ph)
                #self._n_pointings.append(len(save_f))
                #self._pointing_cuts.extend(save_f)

        self.zpoins  = np.array(self.zpoins)
        self.ypoins  = np.array(self.ypoins)
        self.xpoins  = np.array(self.xpoins)
        self.dtpoins = np.array(self.dtpoins)
        # cumulative pointing times for resrponse weighting later
        self.cdtpoins = np.cumsum(self.dtpoins)
    
    
    

class BG():

    def __init__(self,
                 mode='default',
                 dataset=None, # COSIpy.dataset object
                 filename='',
                 tracer=None):

        self.bg_mode = mode
        self.filename = filename
        self.n_time_bins = dataset.times.n_ph#dataset.times.n_time_bins
        self.times_wid = dataset.times.times_wid*2 # half-widths * 2

            
        # elif self.bg_mode == 'sim 6deg despina':
        if self.bg_mode == 'default 6deg':
            print('Reading in flight-average background response for 6 deg CDS binning ...')
            self.default_bg_response_file = '../../data_products/flight_bg_all_v1_fine_6deg.npz'
        
        if self.bg_mode == 'sim 6deg despina':
            print('Reading in simulated Ling-model (1973) background response for 6 deg CDS binning from despina only...')
            self.default_bg_response_file = 'LingModel_bg_all_v1_fine_6deg_despina_only.npz'
            
        else:
            print('Using background mode: '+self.bg_mode)

        # construct / read response
        if self.bg_mode == 'from data':
            self.bg_response = self.construct_bg_response_from_data(dataset)

        elif self.bg_mode == 'from file':
            try:
                with np.load(self.filename) as content:
                    self.bg_response_all = content['bg_response']
                    self.n_bg_ph_per_bin = content['n_bg_ph_per_bin']
                    self.energy_bin_edges_bg = content['energy_bin_edges']

                self.bg_response = self.construct_general_bg_response(dataset)
                    
            except FileNotFoundError:
                print('File '+str(self.filename)+' not found.')
        else:
            try:
                with np.load(self.default_bg_response_file) as content:
                    self.bg_response_all = content['bg_response']
                    self.n_bg_ph_per_bin = content['n_bg_ph_per_bin']
                    self.energy_bin_edges_bg = content['energy_bin_edges']
            
                self.bg_response = self.construct_general_bg_response(dataset)
            
            except FileNotFoundError:
                print('File '+str(self.default_bg_response_file)+' not found.')

                
        # initialise background cuts (default: constant model)
        self.make_bg_cuts(cuts=[])


        # create constant tracer if none is given
        self.tracer = tracer

        if isinstance(self.tracer,(type(None))):
            self.tracer = np.ones(self.n_time_bins)
        else:
            if isinstance(self.tracer,(int, float, complex)):
                self.tracer = np.repeat(self.tracer,self.n_time_bins)
            else:
                if len(self.tracer) != self.n_time_bins:
                    print('Input tracer has '+str(len(self.tracer))+
                          ' entries, but should be '+str(self.n_time_bins)+
                          ', check your input.\nUsing constant tracer instead.')
                    self.tracer = np.ones(self.n_time_bins)
                else:
                    pass
        # repeat (constant) response for all time bins
        self.bg_model = np.repeat(np.expand_dims(self.bg_response,axis=0),
                                  self.n_time_bins,axis=0)

        # normalise to time bins
        #self.bg_model *= (self.times_wid)[:,None,None,None]/np.mean(self.times_wid)
        # TS: not needed?

        # normalise to mean rate in time bin
        self.bg_model *= np.mean(np.sum(dataset.binned_data,axis=(2,3)),axis=0)[None,:,None,None]

        # apply tracer (TS: should be energy dependent in the future)
        self.bg_model *= self.tracer[:,None,None,None]

        # reduce CDS to filter all always-zero elements and speed up fits
        # init list for all energy bins since number of elements will change
        self.bg_model_reduced = []
        # save indices for later use
        self.calc_this = []
        # loop over energy bins
        for i in range(dataset.energies.n_energy_bins):
            # reshape background model to reduce CDS if possible
            # this combines the 3 CDS angles into a 1D array for all times at the chosen energy
            bg_tmp = self.bg_model[:,i,:,:].reshape(self.n_time_bins,
                                                    self.bg_response.shape[1]*self.bg_response.shape[2])

            # get indices of where no entries are there at all (will always be zero, so can be ignored)
            self.calc_this.append(np.where(np.sum(bg_tmp,axis=0) != 0)[0])

            # choose only non-zero entries for all times
            bg_tmp = bg_tmp[:,self.calc_this[i]]

            # append to list for energy bins
            self.bg_model_reduced.append(bg_tmp)
            
    
    def construct_general_bg_response(self,dataset):
        """
        Construct a standard background response from flight data
        :param: dataset  COSIpy.dataset object with binned dataset attribute
        """

        #try:
        
        # array to be filled with response per energy bin
        bg_response_per_bin = np.zeros(dataset.binned_data.shape[1:])

        # making a copy of energy bins in data set to dead with exceptions (oob)
        energy_bin_edges_tmp = np.copy(dataset.energies.energy_bin_edges)

        # loop over energy bins in data set
        for e in range(dataset.energies.n_energy_bins):

            # exception handling for oob
            if (energy_bin_edges_tmp[e] < 143):
                print('Energy bin '+str(e)+' below minimum threshold (143 keV), setting to 143 keV for BG response.')
                energy_bin_edges_tmp[e] = 143.
            if (energy_bin_edges_tmp[e+1] > 5000):
                print('Energy bin '+str(e+1)+' above maximum threshold (2500 keV), setting to 2500 keV for BG response.')        
                energy_bin_edges_tmp[e+1] = 5000.

            # finding indices in finely-binned background response
            idx = np.where((self.energy_bin_edges_bg >= energy_bin_edges_tmp[e]) &
                           (self.energy_bin_edges_bg <= energy_bin_edges_tmp[e+1]))[0]

            # in-between case
            if (len(idx) == 0):
                # only use the last index as representative
                idx = np.where(self.energy_bin_edges_bg <= energy_bin_edges_tmp[e])[0][-1]
                # gets weight one
                weights = 1.
                # is then also response index
                idx_rsp = idx

                # exception when last bin would be included
                if (idx_rsp == 174):
                    idx_rsp = 174-1

            # case where boundaries would fall into neighbouring bins
            elif (len(idx) == 1):

                # gives at most 2 response entries
                weights = np.ones(2)

                # check if data set edges fall onto background edges
                if not (self.energy_bin_edges_bg[idx[0]] == energy_bin_edges_tmp[e+1]):
                    # calculate how much overlap is there
                    delta_e_bg = np.diff(self.energy_bin_edges_bg[idx[0]:idx[0]+2])[0]
                    delta_e_data = np.abs(self.energy_bin_edges_bg[idx[0]]-energy_bin_edges_tmp[e+1])
                    # fractional weighting of last bin
                    ratio_last_bin = delta_e_data/delta_e_bg
                    weights[-1] = ratio_last_bin
                else:
                    # if exactly at boundary, neighbouring bins are weighted zero
                    weights[-1] = 0.

                # same as above just for the first bin
                if not (self.energy_bin_edges_bg[idx[0]] == energy_bin_edges_tmp[e]):
                    delta_e_bg = np.diff(self.energy_bin_edges_bg[idx[0]-1:idx[0]+1])[0]
                    delta_e_data = np.abs(self.energy_bin_edges_bg[idx[0]]-energy_bin_edges_tmp[e])
                    ratio_first_bin = delta_e_data/delta_e_bg
                    weights[0] = ratio_first_bin
                else:
                    weights[0] = 0.
                
                # response index is the same as energy boundary index -1
                idx_rsp = np.concatenate([np.array([idx[0]-1]),idx])

            # data set energy bin matches more than 2 backgorund response bins
            else:

                # then gives one more weight
                weights = np.ones(len(idx)+1)

                # same as above with slightly different boundary handling
                if not (self.energy_bin_edges_bg[idx[0]] == energy_bin_edges_tmp[e]):
                    delta_e_bg = np.diff(self.energy_bin_edges_bg[idx[0]-1:idx[0]+1])[0]
                    delta_e_data = np.abs(self.energy_bin_edges_bg[idx[0]]-energy_bin_edges_tmp[e])
                    ratio_first_bin = delta_e_data/delta_e_bg
                    weights[0] = ratio_first_bin
                else:
                    weights[0] = 0.

                if not (self.energy_bin_edges_bg[idx[-1]] == energy_bin_edges_tmp[e+1]):
                    delta_e_bg = np.diff(self.energy_bin_edges_bg[idx[-1]-1:idx[-1]+1])[0]
                    delta_e_data = np.abs(self.energy_bin_edges_bg[idx[-1]]-energy_bin_edges_tmp[e+1])
                    ratio_last_bin = delta_e_data/delta_e_bg
                    weights[-1] = ratio_last_bin
                else:
                    weights[-1] = 0.

                idx_rsp = np.concatenate([np.array([idx[0]-1]),idx])

                if (idx_rsp[-1] == 174):
                    idx_rsp = np.delete(idx_rsp,-1)
                    weights = np.delete(weights,-1)

            # finally, weighted sum over response entries
            # number of photons is important since we dont have a flat spectrum
            if not isinstance(idx_rsp,np.int64):
                for i in range(len(idx_rsp)):
                    bg_response_per_bin[e,:,:] += self.bg_response_all[idx_rsp[i],:,:]*weights[i]*self.n_bg_ph_per_bin[idx_rsp[i]]
            else:
                bg_response_per_bin[e,:,:] = self.bg_response_all[idx_rsp,:,:]*weights*self.n_bg_ph_per_bin[idx_rsp]
            # normalise to one again
            bg_response_per_bin[e,:,:] /= np.sum(weights*self.n_bg_ph_per_bin[idx_rsp])

        # return to global BG object
        return(bg_response_per_bin)

        print('######################################')
          
        
    def make_bg_cuts(self,
                     cuts):    # list of change points to renormalise the background at

        # add minimum and maximum to allow for empty list
        cuts = [1] + cuts + [1e99]
        # make sure input is somewhat reasonable
        cuts = list(np.unique(cuts))
        self.bg_cuts = np.zeros(self.n_time_bins)
        cidx = 0
        for i in range(1,self.n_time_bins+1):
            if (cuts[cidx] <= i < cuts[cidx+1]):
                self.bg_cuts[i-1] = cuts[cidx]
            else:
                cidx += 1
                self.bg_cuts[i-1] = cuts[cidx]

        self.Ncuts = len(np.unique(self.bg_cuts))
        self.idx_arr = np.ones(self.n_time_bins)
        for i in range(self.Ncuts):
            self.idx_arr[np.where(self.bg_cuts == cuts[i])[0]] = i+1

        self.bg_cuts = self.bg_cuts.astype(int)
        self.idx_arr = self.idx_arr.astype(int)



# # misc. functions required for the stuff to work
def minmax(x):
    """
    Return minimum and maximum of array
    :param: x   Input array
    """
    return np.array([x.min(),x.max()])


def minmin(x):
    """
    Return array shifted to zero
    :param: x   Input array
    """
    return x - x.min()


def polar2cart(ra,dec):
    """
    Coordinate transformation of ra/dec (lon/lat) [phi/theta] polar/spherical coordinates
    into cartesian coordinates
    :param: ra   angle in deg
    :param: dec  angle in deg
    """
    x = np.cos(np.deg2rad(ra)) * np.cos(np.deg2rad(dec))
    y = np.sin(np.deg2rad(ra)) * np.cos(np.deg2rad(dec))
    z = np.sin(np.deg2rad(dec))
    
    return np.array([x,y,z])


def cart2polar(vector):
    """
    Coordinate transformation of cartesian x/y/z values into spherical (deg)
    :param: vector   vector of x/y/z values
    """
    ra = np.arctan2(vector[1],vector[0]) 
    dec = np.arcsin(vector[2])
    
    return np.rad2deg(ra), np.rad2deg(dec)


def construct_scy(scx_l, scx_b, scz_l, scz_b):
    """
    Construct y-coordinate of spacecraft/balloon given x and z directions
    Note that here, z is the optical axis
    :param: scx_l   longitude of x direction
    :param: scx_b   latitude of x direction
    :param: scz_l   longitude of z direction
    :param: scz_b   latitude of z direction
    """
    x = polar2cart(scx_l, scx_b)
    z = polar2cart(scz_l, scz_b)
    
    return cart2polar(np.cross(z,x,axis=0))


def GreatCircle(l1,b1,l2,b2,deg=True):
    """
    Calculate the Great Circle length on a sphere from longitude/latitude pairs to others
    in units of rad on a unit sphere
    :param: l1    longitude of point 1 (or several)
    :param: b1    latitude of point 1 (or several)
    :param: l2    longitude of point 2
    :param: b2    latitude of point 2
    :option: deg  Default True to convert degree input to radians for trigonometric function use
                  If False, radian input is assumed
    """
    if deg == True:
        l1,b1,l2,b2 = np.deg2rad(l1),np.deg2rad(b1),np.deg2rad(l2),np.deg2rad(b2)

    return np.sin(b1)*np.sin(b2) + np.cos(b1)*np.cos(b2)*np.cos(l1-l2)


def angular_distance(l1,b1,l2,b2,deg=True):
    """
    Calculate angular distance on a sphere from longitude/latitude pairs to other using Great circles
    in units of deg
    :param: l1    longitude of point 1 (or several)
    :param: b1    latitude of point 1 (or several)
    :param: l2    longitude of point 2
    :param: b2    latitude of point 2
    :option: deg  option carried over to GreatCricle routine
    """
    # calculate the Great Circle between the two points
    # this is a geodesic on a sphere and describes the shortest distance
    gc = GreatCircle(l1,b1,l2,b2,deg=deg)
    
    # check for exceptions
    if gc.size == 1:
        if gc > 1:
            gc = 1.
    else:
        gc[np.where(gc > 1)] = 1.

    return np.rad2deg(np.arccos(gc))


def verb(q,text):
    """
    Print text of q is True
    :param:   q (boolean)
    :param:   text
    """
    if q:
        print(text)

def circle_on_the_sky(ls,bs,th,n_points=100):
    """
    Returns (galactic) coordinates of a circle with with radius th
    with its centre at a position on a sphere (ls/bs).
    Default are n_points=100 points.
    All angles are to be given in degree!
    
    :param: ls          longitude centre of the circle (deg)
    :param: bs          latitude centre of the circle (deg)
    :param: th          radius of circle on the sky (deg)
    :option: n_points   Default 100, number of points to calculate
    """
    from scipy.spatial.transform import Rotation as R

    thr = np.deg2rad(th)
    
    # start from the circle centre point at galactic coordiantes 0/0 on that sphere
    vec = np.array([np.cos(thr),0,0])
    # rotate that point to the wanted position
    r = R.from_euler('yz',[bs+180,ls+180],degrees=True)
    rot_vec = r.apply(vec)
    # initial and rotated point are NOT UNIT VECTORS, thus normalise when required

    # get points of that circle (radius sin(th), AT position cos(th))
    alpha = np.linspace(-np.pi,np.pi,n_points)
    circle_vec = np.array([np.ones(len(alpha))*np.cos(thr),
                           np.sin(thr)*np.cos(alpha),
                           np.sin(thr)*np.sin(alpha)])
    # rotate these points in the same way
    rot_circle_vec = []
    for i in range(len(alpha)):
        rot_circle_vec.append(r.apply(circle_vec[:,i]))
    rot_circle_vec = np.array(rot_circle_vec).T
    # should not happen, but let's make sure
    rot_circle_vec[2,rot_circle_vec[2,:] < -1] = -1
    
    # calculate l and b coordiantes from (cartesian to spherical on unit sphere)
    b_calc = np.rad2deg(np.arcsin(rot_circle_vec[2,:]/
                                  vector_length(rot_circle_vec[0,:],
                                                rot_circle_vec[1,:],
                                                rot_circle_vec[2,:])))
    l_calc = np.rad2deg(np.arctan2(rot_circle_vec[1,:],rot_circle_vec[0,:]))
    
    return l_calc,b_calc



def vector_length(x,y,z):
    """
    Return Euclidean (L2) norm of a 3D vector (or series of vectors in x/y/z coordiantes):
    
    :param: x    x value(s)
    :param: y    y value(s)
    :param: z    z value(s)
    """
    return np.sqrt(x**2+y**2+z**2)


def find_nearest(array, value):
    """
    Find nearest index for value in array (where for non-existent values)
    :param: array   Input array
    :param: value   value to search for
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

