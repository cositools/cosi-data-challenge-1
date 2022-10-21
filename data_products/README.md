# Data Challenge Simulations

For the first Data Challenge, we wanted to give the users a basic look at COSI data analysis, so we chose 3 straightforward examples based on the COSI science goals:
- Extracting the spectra from the Crab, Cen A, Cyg X-1, and Vela
- Imaging bright point sources, such as the Crab and Cyg X-1
- Imaging diffuse emission from 511 keV and the 26Al 1.8 MeV gamma-ray line

For each of these examples, we have provided a detailed description of the simulated sources and data products here in the data_products directory. Each of the sources was simulated at 10x the astrophysical flux. Having a strong signal simplifies the analysis and allows us to focus on the workflow of the procedures. The COSI SMEX mission is expected to be 50x more sensitive than the balloon-borne mission

The simulations were all performed in MEGAlib, with an accurate mass model of the COSI Balloon instrument. The [COSIBalloon.9Detector.geo.setup](https://github.com/cositools/massmodel-cosi-balloon/blob/main/COSIBalloon.9Detector.geo.setup)  model, which accounts for the failure of three  GeD detectors at different times during flight. Each of the continuum simulations was performed for 100 keV – 10 MeV, and an energy range selection of <5 MeV was used in MEGAlib’s mimrec event selection tool. 

## Data Products:

We have included many combinations of the source simulations and background to allow for flexibility and further testing with these files. All files are either the .npz zipped numpy array format for the CDS-binned response matrix or background model, or the MEGAlib photon list .tra.gz format. MEGAlib was used to perform all of these simulations, and the source models are described in detail below.

Within this directory, there are 2 background models and 3 response matrices. 
- **Scaled_Ling_BG_1x.npz**: background model generated from C. Karwin's scaled 1x Ling background simulation 
- **Scaled_Ling_BG_3x.npz**: background model generated from C. Karwin's scaled 3x Ling background simulation for more statistics
- **Continuum_Response.npz**: 6º response used for spectral analysis and imaging continuum sources
- **511keV_imaging_response.npz**: 6º imaging response required for RL imaging of Galactic positron annihilation
- **1809keV_imaging_response.npz**: 6º imaging response required for RL imaging of Galactic Al-26

There is the full-sky simulation with background and all of the sources combined:
- **DC1_combined_10x.tra.gz**: 4 point sources with 10x flux (Crab, Cyg X1, Cen A, Vela), 511 kev & Al-26 lines, 1x Ling background

There is each of the sources individually with the background:
- **Point_sources_10x_BG.tra.gz**: 4 point sources with 10x flux (Crab, Cyg X1, Cen A, Vela) and 1x Ling background
- **Crab_BG_10x.tra.gz**: Crab with 10x flux and 1x Ling background
- **CenA_BG_10x.tra.gz**: Cen A with 10x flux and 1x Ling background
- **GC511_10xFlux_and_Ling.inc1.id1.extracted.tra.gz**: 511 emission with 10x flux with Ling background
- **DC1_Al26_10xFlux_and_Ling.inc1.id1.extracted.tra.gz**: Al26 emission with 10x flux with Ling background

Each of the sources is also included without background:
- **Point_sources_10x.tra.gz**: 4 point sources with 10x flux (Crab, Cyg X1, Cen A, Vela)
- **Crab_only_10x.tra.gz**: Crab with 10x flux
- **CygX1_only_10x.tra.gz**: Cyg X1 with 10x flux
- **CenA_only_10x.tra.gz**: Cen A with 10x flux
- **Vela_only_10x.tra.gz**: Vela with 10x flux
- **GC_511_10xFlux_only.inc1.id1.extracted.tra.gz**: 511 emission with 10x flux
- **Al26_10xFlux_Only.inc1.id1.extracted.tra.gz**: Al26 emission with 10x flux

And finally, there is 1 background simulation (this is not required for any analysis, but is included here for posterity):
- **Scaled_Ling_BG_1x.tra.gz**: C. Karwin's scaled 1x Ling background simulation

All files have the same start and stop time, in unix time:
start time: 1463443400.0 s
stop time: 1467475400.0 s
total time: 4032000.0 s = 46.67 days

This corresponds to May 17 2016 00:03:20 GMT to Jul 02 2016 16:03:20 GMT, covering the full COSI Balloon flight in 2016.

As described in the main cosi-data-challenge-1 README, these files are stored on Git’s Large File Server, and one needs to install git-lfs to access them.

## Source Models

### Point Sources 

There are four bright point sources that are included in these simulations: the Crab nebula, Cygnus X-1, Centaurus A, and Vela. The spectra and flux for each of these sources was determined from the literature. As explained above, we have used 10x the flux values for these sources in order to simplify the analysis. Please keep that in mind when you see the detections in this Data Challenge.

Crab:   
&emsp; (l,b) = (184.56, -5.78)  
&emsp; Spectral shape: Band function from [Jourdain et al. 2020](https://ui.adsabs.harvard.edu/abs/2020ApJ...899..131J/abstract)  
&emsp; Flux = 0.48977 ph/cm2/s between 100 keV and 10 MeV  
Cen A:   
&emsp; (l,b) = (309.52, 19.42)   
&emsp; Spectral shape: SED from [HESS+LAT collaboration 2018](https://ui.adsabs.harvard.edu/abs/2018A%26A...619A..71H/abstract)  
&emsp; Flux: 0.03609 ph/cm2/s between 100 keV and 10 MeV  
Cyg X1:  
&emsp; (l,b) = (71.33, 3.07)  
&emsp; Spectral shape: SED from [Kantzas+21](https://ui.adsabs.harvard.edu/abs/2021MNRAS.500.2112K/abstract)  
&emsp; Flux = 0.40644 ph/cm2/s between 100 keV - 10 MeV  
Vela:  
&emsp; (l,b) = (263.55, -2.79)  
&emsp; Spectral shape: Power law extrapolation of the Fermi-LAT Vela pulsar (4FGL J0835.3-4510), see [Abdollahi+20](https://iopscience.iop.org/article/10.3847/1538-4365/ab6bcb) </br>
&emsp; Flux = 0.00120 ph/cm2/s between 100 keV - 10 MeV  

The input spectra for these point sources is shown below. Note that the plot shows the true spectra, whereas the flux values reported above are for the 10x simulations. 

<img width="760" alt="Screen Shot 2022-10-17 at 2 17 37 AM" src="https://user-images.githubusercontent.com/33991471/196102631-06aca78d-31da-4363-9f43-b59a6c515786.png">

The simulations are run in MEGAlib’s cosima tool, and then a list-mode image is created in mimrec to confirm the correct point source locations:

![Screen Shot 2022-10-17 at 2 18 30 AM](https://user-images.githubusercontent.com/33991471/196102769-1ac5ffb5-b3e5-45e0-afe1-952c068bf693.png)

### Positron Annihilation at 511 keV

The morphology of the 511 keV emission from positron annihilation is not well constrained. For this first Data Challenge, we have used the model defined in [Knödlseder et al. 2005](https://ui.adsabs.harvard.edu/abs/2005A%26A...441..513K/abstract), where the emission was fit with a 2-D asymmetric Gaussian spatial model with the following parameters (the reported flux is the true value and not the 10x value):

<img width="350" alt="Screen Shot 2022-10-17 at 2 21 06 AM" src="https://user-images.githubusercontent.com/33991471/196103180-3b25375b-3153-44bb-8c1a-a0c68e56585d.png">

The emission was simulated as a 511 keV mono-energetic source, and the image from mimrec confirms the extended emission in the Galactic Center.

![Screen Shot 2022-10-17 at 12 49 54 AM](https://user-images.githubusercontent.com/33991471/196103404-b4f7d089-88bf-4560-be8b-94f0cebb23e8.png)

### Aluminum-26 Decay at 1.8 MeV

The Diffuse Infrared Background Experiment (DIRBE) 240 um map has been shown to be a good tracer for the Al-26 emission, as measured by COMPTEL and INTEGRAL/SPI. We use that distribution as the spatial model for the Al-26 emission, and the inner Galaxy flux was normalized to 3.3×10−4 ph/cm2/s ([Diehl et al. 2006](https://ui.adsabs.harvard.edu/abs/2006Natur.439...45D/abstract)). This model is described further in [Beechert et al. 2022](https://ui.adsabs.harvard.edu/abs/2022ApJ...928..119B/abstract).

![Screen Shot 2022-10-17 at 2 23 44 AM](https://user-images.githubusercontent.com/33991471/196103586-bccf5846-c66f-4a40-bd1c-942e8cfbe672.png)

The emission was simulated with a 1.8 MeV mono-energetic source, and the image from mimrec confirms the extended emission along the Galactic disk.

![Screen Shot 2022-10-17 at 12 49 46 AM](https://user-images.githubusercontent.com/33991471/196103634-818c032e-cae0-40f8-9bc0-b4c997799a9b.png)

Note that the list-mode imaging method employed in MEGAlib is not optimized for diffuse sources, and the extra structure out of the Galactic plane is an imaging artifact.

### Background Radiation

The background radiation model used for the Data Challenge is based on the semi-empirical model from [Ling 1975](https://ui.adsabs.harvard.edu/abs/1975JGR....80.3241L/abstract), which has been used for all COSI balloon analyses. The model describes the angular and atmospheric depth dependance of gamma-rays within 0.3 – 10 MeV and includes 3 primary gamma-ray sources: continuum emission (Bremsstrahlung from cosmic ray interactions in the atmosphere), 511 keV line component (from electron-positron annihilation in the atmosphere, and cosmic diffuse component (down-scattered gamma-rays from the extragalactic background and Galactic plane). The Ling model requires a description of the atmosphere’s density, interaction depth, and mass absorption coefficients, which are obtained through the NRLMSISE-00 model. An altitude of 33.5 km was assumed for these models - altitude drops and changes in longitude and latitude were not taken into account for this simulation.

The amplitude of the Ling background was scaled so the total integrated background spectrum from simulation matched closely to what was measured during flight, as can been seen in the figure below. The flight data (“All Data” label) shows a time-variable background count rate that was influenced by the geomagnetic cutoff, and balloon altitude drops. The cyan line shows the scaled Ling background model. As a reference, the count rate from the 1x flux point sources are shown, and the signal is at most a few percent of the background count rate. 

![Screen Shot 2022-10-17 at 2 04 27 AM](https://user-images.githubusercontent.com/33991471/196103711-13981256-8577-4b1f-9ce8-701874ba10c2.png)

An image of the background simulation traces the exposure map, since the orientation of the COSI Balloon in Galactic coordinates was included in the simulation:

![Screen Shot 2022-10-17 at 2 04 44 AM](https://user-images.githubusercontent.com/33991471/196103767-dcd05934-5b73-48a2-9cbf-74c57520b982.png)

## Flight path and transmission

The source simulations include the real flight aspect information so that the balloon path and source exposure time is accurate. This can be seen in the below plot showing the Galactic longitude and latitude of the zenith direction of the COSI balloon as a function of time.

![Screen Shot 2022-10-17 at 1 55 33 AM](https://user-images.githubusercontent.com/33991471/196103809-2eef45e3-f889-4049-81ec-2a83db5865fb.png)

Furthermore, the transmission probability of the source photons in the atmosphere are calculated for each instance of the simulation. The probability of transmission is taken at a constant altitude of 33 km, and is shown as a function of zenith angle and energy in the below figure.

![Screen Shot 2022-10-17 at 1 56 10 AM](https://user-images.githubusercontent.com/33991471/196103855-2e805235-c4a1-4d82-a568-edf7bca6727e.png)


 

























