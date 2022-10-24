# Spectral fit using the COSI 2016 balloon data. 

For more details on steps of the analysis see the main spectral-fit README [(here)](../../spectral-fit/README.md).

## Data Selection:
The data has been extracted using a $60^\circ$ pointing cut with respect to the position of the Crab. Additionally, we cut the first 21 days of the mission, when the instrument background rate was relatively high and variable. Note that the Crab was mostly outside the field of view during these first 3 weeks. Additionally, we also cut the last day of the mission because of a high background rate. Otherwise, the data selection is identical to that of the simulated data. For more details on the 2016 COSI balloon background rates see [Siegert et al. 2020](https://iopscience.iop.org/article/10.3847/1538-4357/ab9607).

## Background:
For the data challenge the simulated background is constant in time. With the actual balloon data, however, the background rate varies significantly, mostly driven by changes in the balloon's altitude. To accomodate these variations in the fit we apply a so-called "tracer". The tracer is a factor that is applied to the background response, and it effectively normalizes the background rate to the actual data in each time bin. More specifically, for each time bin, the tracer is calculated by summing the counts over the energy and Compton data space dimensions, and then normalizing by the mean of this sum over all time bins. 
