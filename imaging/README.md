# Welcome to image deconvolution with cosipy-classic

In this example, we'll use a Richardson-Lucy deconvolution algorithm to image the full simulated sky as observed during the balloon flight. There are three notebooks in this directory:
1. [RL-DataChallenge-Point_Sources-10XFlux-Ling.ipynb](RL-DataChallenge-Point_Sources-10XFlux-Ling.ipynb) which images 4 simulated point sources.
2. [RL-DataChallenge-511keV_10xFlux-Ling.ipynb](RL-DataChallenge-511keV_10xFlux-Ling.ipynb) performs the diffuse imaging of the 511 keV emission.
3. [RL-DataChallenge-Al26_10xFlux-Ling.ipynb](RL-DataChallenge-Al26_10xFlux-Ling.ipynb) performs the diffuse imaging of the Al-26 1.8 MeV decay line.
The 3 notebooks are almost identical in their execution, but different data sets will be uploaded, different response matrixes will be used, and slightly different imaging parameters will be required. This README should act as a guide offering additional information for each step of the analysis. Please refer to the [data_products](../data_products) README for more details on the scientific background for these sources, and the simulated source models. The following sections align with the different steps within the notebooks.

## Intial Setup

Import relavent packages, define file names, and read in the data files.

## Bin the data

### Define the bins for the Compton data space
**Time bins:** The COSI-balloon instrument freely floated on the balloon platform. This means that, unlike a space or ground-based telescope with well-defined pointings and slewing schedule, its orientation was largely dependent on the unconstrained path of the balloon. It was a zenith-pointing instrument, meaning that its vertical orientation pointed straight above the hanging instrument, towards the balloon above it.

The exception to this freedom is that during the day time, COSI's azimuthal orientation was fixed such that its solar panels remained oriented facing the Sun. At nighttime, though, the instrument freely rotated about its azimuth. 

This is all to say that COSI's orientation (e.g. roll/pitch/yaw) changed rapidly during flight. As such, we might prefer to bin the data into very small (~seconds) time bins to preserve an accurate orientation of the instrument tied to the data. However, this would require massive computational resources. Also, time bins which are too small may contain too few photons for meaningful analysis. 

Through extensive testing, **1800 second** (30 minute) time bins were found to strike a practical balance between a sufficiently precise treatment of instrument orientation and computational means. 

You can feel free to play with the time binning. You may choose to decrease the size to 900 s, or increase it to 3600 s, for example, to see if there's an effect on the image (e.g. does it look less/more blurred, respectively, as you lump more data into more/fewer time bins?)

**Energy bins:** We need to define the energy bins exactly as they are defined in the response.

For point source imaging, we use a continuum response simulation which spans several energy bins across COSI's 0.2-5 MeV bandpass: [150, 220, 325, 480, 520, 765, 1120, 1650, 2350, 3450, 5000] keV.

For positron-electron annihilation keV, we use a response simulation with only one energy bin around the 511 keV signature: **501-521 keV**.

For Al-26, we use a response simulation with only one energy bin around the 1809 keV photopeak signature: **1803-1817 keV**.

**Sky pixel size:** As with the energy binning, the pixel size here must match that of the response. The responses that have been simulated for COSI-balloon assume $6^{\circ} \times 6^{\circ}$ resolution.

### Binning
Calling `.get_binned_data()` will loop through all of the events in the MEGAlib photon list contained in the .tra.gz file to fill the bins of the Compton data sapce.

### Examining the shape

The binned data are contained in "analysis1.dataset.binned_data." This is a 4-dimensional object representing the 5 dimensions of the Compton data space: (time, energy, $\phi$, FISBEL).

The number of bins in each dimension are shown by calling "shape."

Per the binning definitions above, there are 2240 time bins, 10 energy bins for the continuum analysis or 1 bin for line analysis (as governed by those in the response), 30 $\phi$ bins ($\phi$ is the Compton scattering angle; 30 bins of $6^{\circ}$ spanning the full $0-180^{\circ}$ range of possible Compton scattering angles), and 1145 FISBEL bins. 

FISBEL is a unique index which specifies the $\chi$ and $\psi$ dimensions of the Compton Data Space (CDS), whose third dimension is $\phi$. The $\chi$ and $\psi$ dimensions specify the direction of the scattered photon in a Compton interaction. 

How do we end up with 1145 FISBEL bins? Consider a sphere which is $4 \pi( 180^{\circ}/ \pi)^2 = 41252.96 \textrm{ deg}^{2}$ \
Given our $6 \textrm{ deg}^{2}$ binning, we have $41252.96 \textrm{ deg}^{2}$ / $6 \textrm{ deg}^{2}$ $\sim$ 1145 bins.

The notebook will show you how to get the shape of the data set, and how to extract the bin sizes.

### Inspecting the data

Since we have the simulated data read and binned into the Compton data space, we can now make raw spectra, light curves, and other projections of the data. Two examples of this are shown in the notebook. The first is the total raw spectrum from the simulated data set. This is for the duration of the balloon flight (total time = 4031996 seconds = 46.6 days), and the majority of photons in this spectrum are from the background simulation. The 511 keV line, which fills the narrow 4th bin of the full continuum spectrum, is clearly visible and has contributions from the Ling background and from the Galactic center source simulation. The light curve is ...

## Pointing Class

The pointing class handles the aspect information of the balloon platform. 

In the point source notebook, 
Here we see the source move in and out of the field of view. The black dots show when the Crab falls above the horizon, i.e. within COSI's field of view. The Crab is more visible in the latter part of the flight. We notice, too, that the Crab is always somewhat off-axis; it is never directly overhead the instrument at zenith.

COSI's field of view extends ~60 deg from zenith. Hence, the "horizon" lies at the maximum extent of what COSI can see beyond zenith. This means that COSI's zenith lies 60 deg above the horizon.

This is also evident in the above plot: the black dots tracing COSI's zenith pointing don't overlap with the Crab's position, marked by a red star.

## Background Model

We model the background using extensive simulations of Earth's atmospheric $\gamma$-ray background. The simulations assume the Ling model of atmospheric $\gamma$-ray emission ([Ling 1975](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/JA080i022p03241), which is often adopted for this purpose in MeV $\gamma$-ray astrophysics experiments.

The simulations use an accurate mass model of the COSI-balloon instrument during flight and follow the true orientation of the instrument as it traveled along its flight path.

Notably, this background model excludes the significant background from instrumental activation. Instrumental activation refers to the excitation of instrument materials by bombarding high-energy particles, e.g. cosmic rays. The instrument materials subsequently de-excite via the emission of $\gamma$-rays which fall in COSI's energy bandpass. Often, the de-excitation lines exactly overlap with astrophysical lines of interest, including 511 keV and 1809 keV. 

A complete treatment of the background, therefore, would include both atmospheric and instrumental activation simulations. For simplicity, however, in this imaging tutorial we model only atmospheric background. Future data challenges will include instrumental activation.

## Read in response matrix

The response is created through large simulations in MEGAlib's ResponseCreator program. After we read in the response, which is a 5D numpy array (.npz format) ...

The shape of the response spans (Galactic latitude $b$, Galactic longitude $\ell$, Compton scattering angle $\phi$,  FISBEL, energy). The size of each dimension depends on the chosen pixel size. Here, we've chosen $6^{\circ}$ pixels. 

Galactic latitude $b \in [-90^{\circ}, 90^{\circ}] \rightarrow$ 30 bins.\
Galactic longitude $\ell \in [-180^{\circ}, 180^{\circ}] \rightarrow$ 60 bins.\
Compton scattering angle $\phi \in [0^{\circ}, 180^{\circ}] \rightarrow$ 30 bins ("analysis1.dataset.phis.n_phi_bins").\
See above for explanation of 1145 FISBEL bins ("rsp.rsp.n_fisbel_bins").\
The continuum response has 10 energy bins.

It's also helpful at this time to remind ourselves of the shapes of the data and background model.

The shape of the data and background objects span (time, energy, Compton scattering angle, FISBEL).

Given the time bin size "Delta_T" which we defined at the beginning of the notebook, there are 2240 time bins.

## RL Imaging

Up until this point, the analysis steps have been completely generic and paraellel what is required for standard COSI analysis. The spectral fitting notebook is almost identical for the above steps. Everything that comes next is specific for the image deconvolution algorithims, and as you can see the RL algorithim and supporting functions are hard-coded in this notebook. These tools will be ....

Many of the comments in the code within the notebook can help shed light on the various steps, but...

### Setup for imaging

Defining a grid on the sky to make images, linear 6x6 grid, and convert the grid to zenith and azimuth pairs for each pointing.

Define a function to get the respones of an image for arbitraty time binning. 

We're selecting only one energy bin here for the full analysis 320-480 keV for the continuum since it has the most counts.
And then reduce the response dimensions to only use non-zero bins. The shape after this shows the reduction and re-shape to latitude x longitude x CDS bins.

Define the function to get the image response for arbitrary time binning. The convolves the pointings and the response, to fill the CDS with the expected
Calculate the general image response for the current data set - ...

### Exposure map

Here, we plot the explosure map, which we're representing with the response weighted by the time.

### RL Algorithim

The steps follow the algorithm as outlined in [Knödlseder et al. 1999](https://ui.adsabs.harvard.edu/abs/1999A%26A...345..813K/abstract). Refer to that paper for a mathematical description of the algorithm.

The total memory used during these iterations is about 94 GB for the continuum response, 75 GB for the 511 keV, and 105 GB for the Al-26!! You might be limited on your personal computer, and at the very least, you might not be able to do much else with your machine while this is running. 

#### Adjustable parameters
There are three parameters at the beginning of this RL cell which we encourage you to adjust. In fact, it is often necessary to adjust these parameters depending on the data being studied.

- map_init\
This is the flux value of the initial, isotropic map. Typically, a value of 0.01 works well. For stronger sources, you can try increasing it to 0.1 or 1.0. As an example, running the algorithm on a source-only (no BG) simulation of the Crab, Cen A, Cygnus X-1, and Vela works well with map_init = 0.01. However, when imaging these sources each simulated with 10X their true flux values, the algorithm fails at 0.01 and work when map_init = 1.0.

- iterations\
This is the number of RL iterations. You can set this to a small value, say 50, as you get used to using the algorithm. In our testing, though, for fully converged images we usually let the algorithm run for 150 iterations. ***This can take anywhere from several hours (usually simulations without background) to overnight (simulations with background) to run.***

- afl_scl\
This is a scaling factor for the delta map which we call the "acceleration parameter." This allows the delta map to be afl_scl times stronger than the original RL algorithm suggests (c.f. Knoedlseder+1997).\
The default value here is 2000, though 1000 also works well. If you find that the algorithm returns "Fit failed" messages after running for awhile, for example, lowering this acceleration parameter to 1000 can help.

Other parameters you can adjust:
- mu_Abg, sigma_Abg\
There is a prior in the background fit defined by mu_Abg +/- sigma_Abg. By default, mu_Abg and sigma_Abg are set to fitted_bg and most testing has been done with this setting. You can try constraining the fit by decreasing sigma_Abg, for example, to sigma_Abg/2., sigma_Abg/10., which would enable to fit to vary by 50%, 10% of the initial guess.

- delta_map_tot_old\
You can change the exponent of the denominator. By default, it is set to 0.25 to help avoid exposure edge effects. All testing has been done with this fourth root. However, you can try setting it to 0, 0.5, etc. to see what happens. You can also try smoothing delta_map_tot_old with a Gaussian filter.


### Truncating the response at $90^{\circ}$
As discussed above, COSI's field of view extends $60^{\circ}$ beyond its zenith. The data/simulations themselves, however, only have hard cut manually applied at $90^{\circ}$: this is the Earth Horizon Cut (EHC). The EHC is applied to the data in order to remove background $\gamma$-ray emanating from the Earth's atmosphere below the instrument. 

To preserve any photons which may scatter in just beyond COSI's nominal $60^{\circ}$ field of view but not beyond the $90^{\circ}$ EHC, we define a "cut" which trucates, i.e. zeroes out, the response at $90^{\circ}$. 

Note that because the EHC has removed all photons beyond $90^{\circ}$, setting the cut to a value greater than $90^{\circ}$ will behave identically to cut = $90^{\circ}$.


