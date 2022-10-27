# COSI 2016 Balloon Data

Here we will perform a spectral fit and imaging of the Crab using the 2016 COSI balloon data. The analyses will mostly follow the notebooks used for the simulated data. For more details on steps of the spectral analysis see the main spectral-fit README [(here)](../../spectral-fit/README.md), and likewise for the imaging analysis README [(here)](../../imaging/README.md). This README summarizes the key differences that are needed for analyzing the balloon data. 

## The Balloon Flight
COSI launched from Wanaka, New Zealand on NASA's Super Pressure Balloon on May 16th, 2016. After a 46 day flight, COSI made a landing in Southern Peru. The trajectory of the mission is shown in the image below. COSI experienced 3 high-voltage detector failures, the last of which is indicated with the start of the red curve. For more details of the COSI Balloon instrument see ([Kierans et al. 2017](https://ui.adsabs.harvard.edu/abs/2017arXiv170105558K/abstract)).

<img width="500" alt="Screen Shot 2022-10-26 at 7 06 49 PM" src="https://user-images.githubusercontent.com/54562666/198155231-b5653624-41aa-45c6-9e44-d7a1c818a974.png"><img width="337" alt="Screen Shot 2022-10-26 at 7 09 08 PM" src="https://user-images.githubusercontent.com/54562666/198155463-1068d09e-cf94-4677-9a71-35c80a037aa7.png">


## Data Selection
The data has been extracted using a $60^\circ$ pointing cut with respect to the position of the Crab. In practice, this means that we only use times for which the location of the Crab was within $60^\circ$ from COSI's zenith. This can be seen in the images below. The plot on the left shows the pointing in Galactic coordinates as a function of time for the simulated data (which doesn't use a pointing cut), and the plot on the right is for the balloon data, with the $60^\circ$ pointing cut. The green star shows the position of the Crab.

<img width="450" alt="Screen Shot 2022-10-26 at 6 40 01 PM" src="https://user-images.githubusercontent.com/54562666/198152120-ab33efe4-8c6e-4515-bae9-e7f781a6fe5a.png"><img width="450" alt="Screen Shot 2022-10-26 at 6 29 17 PM" src="https://user-images.githubusercontent.com/54562666/198151093-1e4ca7bd-c5b4-4c2e-8731-78b454eb7d21.png">

Additionally, we cut the first 21 days of the mission, when the instrument background rate was relatively high and variable. Note that the Crab was mostly outside the field of view during these first 3 weeks. Additionally, we also cut the last day of the mission because of high background rates. This can be seen in the plots below, which show the elevation of the Crab above COSI's horizon. The left plot is for the simulated data, using the entire flight time, and the right plot is for the balloon data, with the time cuts applied. As can be seen, the Crab was in COSI's field of view mainly during the latter part of the mission. 

<img width="500" alt="Screen Shot 2022-10-26 at 6 42 54 PM" src="https://user-images.githubusercontent.com/54562666/198152487-a08bf99c-9b05-460a-8561-c57f4d5cfcc6.png"><img width="500" alt="Screen Shot 2022-10-26 at 6 44 29 PM" src="https://user-images.githubusercontent.com/54562666/198152633-b61e5bd8-d7e6-4c75-b3df-5649ac868c22.png">

Otherwise, the data selection is identical to that of the simulated data. For more details on the 2016 COSI balloon background rates see [Siegert et al. 2020](https://iopscience.iop.org/article/10.3847/1538-4357/ab9607).


## Background
For the data challenge the simulated background is constant in time. With the actual balloon data, however, the background rate varies significantly, mostly driven by changes in the latiude (i.e. geomagnetic cutoff) and the altitude (i.e. atmospheric depth). Additionally, the application of the pointing cut produces drastic variation in the light curve as the Crab moves in and out of COSI's field of view. To accomodate these variations, in the fit we apply a so-called "tracer". The tracer is a factor that is applied to the background response, and it effectively normalizes the background rate to the actual data in each time bin. More specifically, for each time bin, the tracer is calculated by summing the counts over the energy and Compton data space dimensions, and then normalizing by the mean of this sum over all time bins. The difference between the background rates and models without and with a tracer can be seen below, for the simuted data analysis (left) and balloon data analysis, respectively. 

<img width="500" alt="" src="https://user-images.githubusercontent.com/54562666/198153505-cfdb30f7-0503-46f1-9e3a-2426a4470d18.png"><img width="500" alt="Screen Shot 2022-10-26 at 6 53 00 PM" src="https://user-images.githubusercontent.com/54562666/198153646-396fb0d5-ea2f-4cf4-890a-755982b10d26.png">

## Results of the Spectral Fit

Results for the Crab fit are shown below. Upper limits are shown for bins with SNR < 3. The fit will output a .dat file, where you can examine the numerical results, including the SNR of each bin. We have verified that bins 4, 5, and 6 are consistent with the physical values (based on our simulated Crab spectrum). The other bins (and most notably, bin 3) are within roughly a factor of ~2 of expectations. It is important to note, however, that the extracted spectrum has a significant dependence on the background model being used. Here we are using a highly simplified constant Ling model (with tracer). Additionally, we haven't accounted for activation backgrounds, nor other astrophysical sources (i.e. 511, $^{26}$ Al, Cen A, etc.). Thus, the precision of this result should not be over-interpreted at this point. The COSI team is currently working on developing improved response handling, background models, and analysis tools, which will be part of future Data Challenges!  

<img width="500" alt="Screen Shot 2022-10-26 at 8 21 03 PM" src="https://user-images.githubusercontent.com/54562666/198162612-3b02b0f7-058b-4c76-8d5a-2a6ca84f8150.png">

## Results of the Imaging
The primary outputs of the imaging, after 100 iterations, are shown below. Note that there is an artifact in the southern pole, where bright emission can be seen. This feature should be surpressed with more iterations. The cyan star in the map shows the position of the Crab, and as can be seen, it corresponds to the brightest pixel in the image (and the dimmer surrounding pixels). 

<img width="450" alt="Screen Shot 2022-10-26 at 8 44 09 PM" src="https://user-images.githubusercontent.com/54562666/198165282-fd87f5c6-bd25-49f5-b559-739dbcebd5c3.png"><img width="450" alt="Screen Shot 2022-10-26 at 8 44 24 PM" src="https://user-images.githubusercontent.com/54562666/198165010-21b027f3-a6c1-472a-93d4-4966f74db219.png">


<img width="700" alt="Screen Shot 2022-10-26 at 8 42 28 PM" src="https://user-images.githubusercontent.com/54562666/198164730-01d29ed7-eb36-4526-9617-0b7ceefca17e.png">



