# COSI 2016 Balloon Data

Here we will perform a spectral fit and imaging of the Crab using the 2016 COSI balloon data. The analyses will mostly follow the notebooks used for the simulated data. For more details on steps of the spectral analysis see the main spectral-fit README [(here)](../../spectral-fit/README.md), and likewise for the imaging analysis README [(here)](../../imaging/README.md). This README summarizes the key differences that are needed for analyzing the balloon data. 

## The Balloon Flight

## Data Selection
The data has been extracted using a $60^\circ$ pointing cut with respect to the position of the Crab. In practice, this means that we only use times for which the location of the Crab was within $60^\circ$ from COSI's zenith. This can be seen in the images below. The plot on the left shows the pointing in Galactic coordinates as a function of time for the simulated data (which doesn't use a pointing cut), and the plot on the right is for the balloon data, with the $60^\circ$ pointing cut. The green star shows the position of the Crab.

<img width="450" alt="Screen Shot 2022-10-26 at 6 40 01 PM" src="https://user-images.githubusercontent.com/54562666/198152120-ab33efe4-8c6e-4515-bae9-e7f781a6fe5a.png"><img width="450" alt="Screen Shot 2022-10-26 at 6 29 17 PM" src="https://user-images.githubusercontent.com/54562666/198151093-1e4ca7bd-c5b4-4c2e-8731-78b454eb7d21.png">

Additionally, we cut the first 21 days of the mission, when the instrument background rate was relatively high and variable. Note that the Crab was mostly outside the field of view during these first 3 weeks. Additionally, we also cut the last day of the mission because of high background rates. This can be seen in the plots below, which show the elevation of the Crab above COSI's horizon. The left plot is for the simulated data, using the entire flight time, and the right plot is for the balloon data, with the time cuts applied. As can be seen, the Crab was in COSI's field of view mainly during the latter part of the mission. 

<img width="500" alt="Screen Shot 2022-10-26 at 6 42 54 PM" src="https://user-images.githubusercontent.com/54562666/198152487-a08bf99c-9b05-460a-8561-c57f4d5cfcc6.png"><img width="500" alt="Screen Shot 2022-10-26 at 6 44 29 PM" src="https://user-images.githubusercontent.com/54562666/198152633-b61e5bd8-d7e6-4c75-b3df-5649ac868c22.png">

Otherwise, the data selection is identical to that of the simulated data. For more details on the 2016 COSI balloon background rates see [Siegert et al. 2020](https://iopscience.iop.org/article/10.3847/1538-4357/ab9607).


## Background
For the data challenge the simulated background is constant in time. With the actual balloon data, however, the background rate varies significantly, mostly driven by changes in the balloon's altitude and geomagnetic cutoff. To accomodate these variations in the fit we apply a so-called "tracer". The tracer is a factor that is applied to the background response, and it effectively normalizes the background rate to the actual data in each time bin. More specifically, for each time bin, the tracer is calculated by summing the counts over the energy and Compton data space dimensions, and then normalizing by the mean of this sum over all time bins. 
