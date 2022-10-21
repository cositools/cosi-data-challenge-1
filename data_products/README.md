# Data Challenge Simulations

For the first Data Challenge, we wanted to give the users a basic look at COSI data analysis, so we chose 3 straightforward examples based on the COSI science goals:
- Extracting the spectra from the Crab, Cen A, Cyg X-1, and Vela
- Imaging bright point sources, such as the Crab and Cyg X-1
- Imaging diffuse emission from 511 keV and the 26Al 1.8 MeV gamma-ray line

For each of these examples, we have provided a detailed description of the simulated sources and data products here in the data_products directory. Each of the sources was simulated at 10x the astrophysical flux. Having a strong signal simplifies the analysis and allows us to focus on the workflow of the procedures. The COSI SMEX mission is expected to be 50x more sensitive than the balloon-borne mission

The simulations were all performed in MEGAlib, with an accurate mass model of the COSI Balloon instrument. The [COSIBalloon.9Detector.geo.setup](https://github.com/cositools/massmodel-cosi-balloon/blob/main/COSIBalloon.9Detector.geo.setup)  model, which accounts for the failure of three  GeD detectors at different times during flight. Each of the continuum simulations was performed for 100 keV – 10 MeV, and an energy range selection of <5 MeV was used in MEGAlib’s mimrec event selection tool. 

## Data Products:

We have included many combinations of the source simulations and background to allow for flexibility and further testing with these files. All files are either the .npz zipped numpy array format for the CDS-binned response matrix or background model, or the MEGAlib photon list .tra.gz format. MEGAlib was used to perform all of these simulations, and the source models are described in detail below.

Within this directory, there is 1 background model and 3 response matrices. 
- **Scaled_Ling_BG_1x.npz**: background model generated from C. Karwin's scaled 1x Ling background simulation 
- **Continuum_Response.npz**: 6º response used for spectral analysis and imaging continuum sources
- **511keV_imaging_response.npz**: 6º imaging response required for RL imaging of Galactic positron annihilation
- **1809keV_imaging_response.npz**: 6º imaging response required for RL imaging of Galactic Al-26

There is the full-sky simulation with background and all of the sources combined:
- **DC1_combined_10x.tra.gz**: 4 point sources with 10x flux (Crab, Cyg X1, Cen A, Vela), 511 kev & Al-26 lines, 1x Ling background

There is each of the sources individually with the background:
- **Point_sources_10x_BG.tra.gz**: 4 point sources with 10x flux (Crab, Cyg X1, Cen A, Vela) and 1x Ling background
- **Crab_BG_10x.tra.gz**: Crab with 10x flux and 1x Ling background
- **CygX1_BG_10x.tra.gz**: Cyg X1 with 10x flux and 1x Ling background
- **CenA_BG_10x.tra.gz**: Cen A with 10x flux and 1x Ling background
- **Vela_BG_10x.tra.gz**: Vela with 10x flux and 1x Ling background
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

All files have the same start and stop time, in unix time: \
start time: 1463443400.0 s \
stop time: 1467475400.0 s \
total time: 4032000.0 s = 46.67 days \

This corresponds to May 17 2016 00:03:20 GMT to Jul 02 2016 16:03:20 GMT, covering the full COSI Balloon flight in 2016.

As described in the main cosi-data-challenge-1 README, these files are stored on Git’s Large File Server, and one needs to install git-lfs to access them.

## Science Background

### Point Sources

There are four bright point sources that are included in these simulations: the Crab nebula, Cygnus X-1, Centaurus A, and Vela.

The Crab nebula is considered both a pulsar wind nebula (PWN) and supernova remnant. The PWN surrounds the Crab Pulsar, a rapidly rotating and magnetized neutron star in the Milky Way constellation Taurus. The supernova remnant was produced by SN 1054. The Crab entered COSI-balloon's field of view for only ~12 of the 46 days of the 2016 flight ([Sleator 2019](https://www.proquest.com/docview/2313733159?pq-origsite=gscholar&fromopenview=true)); the balloon remained largely in Earth's Southern Hemisphere and the Crab is more easily viewed from the Northern Hemisphere. Nevertheless, as the brightest persistent $\gamma$-ray source in the sky, the Crab is detectable in the balloon data and is imaged in this notebook with 10x its true flux of 0.049 ph cm $^{-2}$ s $^{-1}$ (100 keV-50 MeV).

Cygnus X-1 is a bright hard X-ray source in the Cygnus constellation of the Milky Way. It is believed to be a black hole in an X-ray binary system. Cygnus X-1 emits in COSI's bandpass as well, and like the Crab is simulated here at 10x its true flux of 0.041 ph cm $^{-2}$ s $^{-1}$ (100 keV-50 MeV). This data challenge thus helps establish expectations for COSI-balloon observations of Cygnus X-1 during the 2016 flight.

Centaurus A is a galaxy in the constellation of Centaurus. Also called NGC 5128, Centaurus A has a supermassive black hole which emits X-rays and radio waves. This notebook attempts to image Centaurus A as seen with 10x its true flux of 0.0036 ph cm $^{-2}$ s $^{-1}$ (100 keV-50 MeV) during the 2016 balloon flight.

The Vela supernova remnant (Type II supernova in the constellation Vela) is the final point source inlcuded at 10x its true flux of 0.00014 ph cm $^{-2}$ s $^{-1}$ (100 keV-50 MeV). It is fainter than the other three sources but is included as a source which could be observed by an instrument like COSI-balloon with more observation time.

### Positron Annihilation at 511 keV

keV emission from the center of the Milky Way Galaxy. As MeV gamma-ray instruments, COSI-balloon and the COSI satellite are uniquely equipped to study this signal, which traces one of the biggest unsolved mysteries in gamma-ray astrophysics. 

For context, a strong 511 keV signal emanating from the direction of the Galactic Center $(\ell = 0^{\circ}, b = 0^{\circ})$ was first discovered in the 1970s on a series of balloon missions [(Johnson III et al. 1972,](https://adsabs.harvard.edu/full/1972ApJ...172L...1J?TB_iframe=true&width=370.8&height=658.8) [ Leventhal et al. 1978,](https://adsabs.harvard.edu/full/record/seri/ApJ../0225/1978ApJ...225L..11L.html) [ Johnson III & Haymes 1973,](https://adsabs.harvard.edu/pdf/1973ApJ...184..103J) [ Haymes et al. 1975,](https://adsabs.harvard.edu/full/record/seri/ApJ../0201/1975ApJ...201..593H.html) [ Ling et al. 1977,](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/JA082i010p01463) [ Albernhe et al. 1981)](https://adsabs.harvard.edu/pdf/1981A%26A....94..214A). Despite decades of observation since the initial measurements, several questions remain regarding the nature of this abundant positron-electron annihilation. 

1. What is the underlying nature of the emission morphology?

The 511 keV image from the INTEGRAL SPI satellite ([Bouchet et al. 2010](https://iopscience.iop.org/article/10.1088/0004-637X/720/2/1772/meta), shown below as a significance map) reveals two seemingly distinct components: an extended "disk" which traces the Galactic Plane and a central, bright "bulge" around the Galactic Center. The notable difference between these structures is not understood. The extended disk emission may suggest that positrons are propogating away from and annihilating at a distance from their production sites, thereby smearing out the emission into a diffuse presentation. It is also possible, however, that there could be a collection of point-like sources emitting positrons which annihilate close to their progenitors; it is the collection of these sources together which may form a total diffuse structure. However, no point source of positrons has been detected. It is likely that positrons propogate away from their production sites and slow down to low enough energies to form positronium.

Imaging with fine angular resolution can help disentangle this dichotomy. The detection of an individual point source, for example, could revolutionize our understanding of the morphology. High-resolution spectroscopy is also critical to understanding the transport and eventual annihilation sites of these positrons: a significant ortho-Positronum (o-Ps; $\leq 511$ keV) component of the spectrum would suggest annihilation in colder regions of the interstellar medium (ISM). A stronger signature at the line energy of 511 keV (para-Positronium, p-Ps) indicates annhilation in warmer regions of the ISM. Data from SPI currently favor the latter scenario, though additional measurements of this o-Ps to p-Ps fraction are necessary.

<img width="619" alt="Bouchet_2010_SPI_511keV" src="https://user-images.githubusercontent.com/33991471/196822990-8698d761-bce7-4b5e-8ae7-3aca6163380f.png">


2. Where are all of these positrons coming from?

The bulge emission exhibits a positron annihilation rate of $10^{43} e^+ s^{-1}$. Positrons from the $\beta^+$ decay of nucleosynthesis products may account for the $\sim 10^{42} e^+ s^{-1}$ in the disk and some of the bulge, but the origin of the remaining postirons in the bulge is unknown. Important to note is that there is no lack of potential positron sources; rather, there are too many possible sources to explain the emission! Positrons are readily created in a wide variety of astrophysical objects and processes. Stellar flares, massive stars, supernovae (core-collapse and Type Ia), classical novae, and neutron star mergers all synthesis radioactive isotopes which can $\beta^+$ decay. Secondary interactions with cosmic rays produce positrons, positrons are created through pair creation channels in strong photon and magnetic fields, evaporating black holes may potentially produce positrons, dark matter density profiles may trace positron annihilation,...and more. 

Ultimately, we can boil the "positron puzzle" down to uncertainty around the *source, transport, and sink* (the annihilation itself) of these particles. 

The COSI-balloon flight in 2016 detected the 511 keV signature of the positron puzzle with $7.2\sigma$ significance [(Kierans et al. 2020)](https://iopscience.iop.org/article/10.3847/1538-4357/ab89a9/meta). It also clearly imaged the bright bulge emission near the Galactic Center [(Siegert et al. 2020)](https://iopscience.iop.org/article/10.3847/1538-4357/ab9607/meta). Both measurements are consistent with those from SPI and indicate additional extended emission. 

In this notebook, you will image the Galactic 511 keV emission, simulated at 10X its true flux (10X flux = $1.1 \times 10^{-2}$ ph cm $^{-2}$ s $^{-1}$) for robust statistics, as seen during the COSI-balloon flight in 2016. You should expect to see the bright bulge emission near the Galactic Center (but not the disk emission, which is too weak for COSI-balloon to see during its 46-day flight; SPI was able to detect the disk's $\sim 1$ ph/week with over a decade of observation time).


### Aluminum-26 Decay at 1.8 MeV

As MeV gamma-ray instruments, COSI-balloon and the COSI satellite are uniquely equipped to study the 1.809 MeV signature emission from this radioisotope, which traces stellar nucleosynthesis over millions of years. 

In the 1980s, NASA's High Energy Astrophysical Observatory (HEAO-3) satellite mission detected 1.809 MeV emission emanating from the direction of the Galactic Center [Mahoney et al. 1984](https://adsabs.harvard.edu/pdf/1984ApJ...286..578M). This marked the discovery of Galactic Al-26. The spectrum is shown below.

<img width="350" alt="Mahoney_1984_HEAO-3_1809keV" src="https://user-images.githubusercontent.com/33991471/196823369-33d78bd1-b84e-49a1-af73-fc0e2168a5ab.png">

Subsequent observations by the Compton telescope (COMPTEL) on-board NASA's Compton Gamma Ray Observatory (CGRO) yielded the first image of Al-26 emission ([Oberlack et al. 1996](https://ui.adsabs.harvard.edu/abs/1996A%26AS..120C.311O/abstract), [Oberlack 1997](https://ui.adsabs.harvard.edu/abs/1997PhDT.........8O/abstract), [Pluschke et al. 2001](https://ui.adsabs.harvard.edu/abs/2001ESASP.459...55P/abstract)). Emission is concentrated in the Inner Galaxy $(|\ell| \leq 30^{\circ}, |b| \leq 10^{\circ}$) with enhanced emission in regions of massive star activity, including Cygnus, Carina, and Vela. 

The SPectrometer on INTEGRAL (SPI) largely corroborated the features seen in the COMPTEL image with over a decade of observation time from ESA's INTEGRAL satellite ([Bouchet et al. 2015](https://iopscience.iop.org/article/10.1088/0004-637X/801/2/142/meta)). Emission is concentrated in the Inner Galaxy with a reported flux of $\sim 3.3 \times 10^{-4}$ ph cm $^{-2}$ s $^{-1}$. As in the COMPTEL image, there is enhanced emission in regions of massive star activity, including Perseus/Taurus, Cygnus/Cepheus, Carina, Vela, and Scorpius-Centaurus. 

<img width="400" alt="COMPTEL_1 8MeV_image" src="https://user-images.githubusercontent.com/33991471/196823447-6373d7da-1ddd-4a72-bfa6-465bebf1013c.png"> | <img width="400" alt="SPI_1 8MeV_image" src="https://user-images.githubusercontent.com/33991471/196823463-61772883-35ed-4151-8cab-c69872f75a7c.png"> \
COMPTEL 1.8 MeV image (Pluschke et al. 2001) | SPI 1.8 MeV image (Bouchet et al. 2015)

The COSI-balloon flight in 2016 measure the 1.809 MeV signature of Al-26 with $3.7\sigma$ significance, corresponding to about 106 Al-26 photons [(Beechert et al. 2022)](https://iopscience.iop.org/article/10.3847/1538-4357/ac56dc/meta). The reported Inner Galaxy flux of $(8.6 \pm 2.5) \times 10^{-4}$ ph cm $^{-2}$ s $^{-1}$ and line centroid of $1811.2 \pm 1.8$ keV are consistent with results from SPI and COMPTEL within 2 $\sigma$ uncertainties. Future observations with the COSI satellite (significantly increased observation time, greater effective area at 1.8 MeV, better constraints on high-latitude emission, and finer angular resolution) will comprise an important comparison to the balloon measurement, which was the first measurement of Al-26 on a compact Compton telescope.

Furthermore, the COSI satellite's full-sky observations with fine angular resolution have potential to more closely study individual regions of massive star activity; in particular, resolving individual sites of emission within Cygnus is a promising goal of the mission. Detailed imaging and spectrosocpic studies of the region may inform better understanding of the dynamics of Al-26 after it is produced and ejected from massive stars. 

In this notebook, you will image the Galactic Al-26 emission (traced by the DIRBE 240 $\mu m$ image) as seen during the COSI-balloon flight in 2016. The flux is simulated at 10X the observed Al-26 flux (10X Inner Galaxy flux = $3.3 \times 10^{-3}$ ph cm $^{-2}$ s $^{-1}$; 10X total map flux = $1.2 \times 10^{-2}$ ph cm $^{-2}$ s $^{-1}$) for robust statistics. You should expect to see extended emission along the Galactic Plane, similar to that revealed by COMPTEL and SPI's 1.8 MeV images. The massive star regions of Cygnus, Carina, and Vela will not be as easily identifiable in this simulation of the balloon flight; for those, be sure to participate in the data challenge (and real data analysis) of the COSI satellite!



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


 

























