import numpy as np
import pandas as pd
from scipy.integrate import trapz

def powerlaw(energy, p):
    """
    Returns:
    Differential flux shaped as power-law function "in units of ph/cm2/s/keV"
    
    Parameters:
    :param energy: 1D-array of energies where power-law is evaluated (in keV)
    :param p:      p[0] = Normalisation (in ph/cm2/s/keV)
    :param p:      p[1] = Power-law index (unitless)

    Misc:
    E0 = 100 keV  Pivotal energy, i.e. where the power-law is normalised at (in keV)
    """
    
    E0 = 100.
    
    return p[0]*np.power(energy/E0,p[1])


def powerlaw_boxcar(energy, p):
    """
    Returns:
    Differential flux shaped as power-law function "in units of ph/cm2/s/keV"
    
    Parameters:
    :param energy: 1D-array of energies where power-law is evaluated (in keV)
    :param p:      p[0] = Normalisation (in ph/cm2/s/keV)
    :param p:      p[1] = Power-law index (unitless)
    :param p:      p[2] = Lower energy bound
    :param p:      p[3] = Higher energy bound

    Misc:
    E0             Pivotal energy, i.e. where the power-law is normalised at (in keV)
    """
    
    E0 = 100.

    return p[0]*np.power(energy/E0,p[1])*np.heaviside(energy-p[2],0)*np.heaviside(p[3]-energy,0)


def cutoff_powerlaw(energy,p):
    """
    Returns:
    Differential flux shaped as power-law function with high-energy
    cut-off "in units of ph/cm2/s/keV"
    
    Parameters:
    :param energy: 1D-array of energies where power-law is evaluated (in keV)
    :param p:      p[0] = Normalisation (in ph/cm2/s/keV)
    :param p:      p[1] = Power-law index (unitless)
    :param p:      p[2] = Cut-off energy (in keV)

    Misc:
    E0 = 100 keV   Pivotal energy, i.e. where the power-law is normalised at (in keV)
    """
    
    E0 = 100.
    
    return p[0]*np.power(energy/E0,p[1])*np.exp(-energy/p[2])


def hilo_cutoff_powerlaw(energy,p):
    """
    Returns:
    Differential flux shaped as power-law function with high- and
    low-energy cut-off "in units of ph/cm2/s/keV"
    
    Parameters:
    :param energy: 1D-array of energies where power-law is evaluated (in keV)
    :param p:      p[0] = Normalisation (in ph/cm2/s/keV)
    :param p:      p[1] = Power-law index (unitless)
    :param p:      p[2] = Low-energy cut-off energy (in keV)
    :param p:      p[3] = High-energy cut-off energy (in keV)

    Misc:
    E0 = 100 keV   Pivotal energy, i.e. where the power-law is normalised at (in keV)
    """
    
    E0 = 100.    
    
    return p[0]*np.exp(-p[2]/energy)*np.power(energy/E0,p[1])*np.exp(-energy/p[3])


def hilo_supercutoff_powerlaw(energy,p):
    """
    Returns:
    Differential flux shaped as power-law function with super-exponential
    high- and low-energy cut-off "in units of ph/cm2/s/keV"
    
    Parameters:
    :param energy: 1D-array of energies where power-law is evaluated (in keV)
    :param p:      p[0] = Normalisation (in ph/cm2/s/keV)
    :param p:      p[1] = Power-law index (unitless)
    :param p:      p[2] = Low-energy cut-off energy (in keV)
    :param p:      p[3] = Exponential low-energy cut-off index (unitless)
    :param p:      p[4] = High-energy cut-off energy (in keV)
    :param p:      p[5] = Exponential high-energy cut-off index (unitless)
    
    Misc:
    E0 = 100 keV   Pivotal energy, i.e. where the power-law is normalised at (in keV)
    """    
        
    E0 = 100.    
        
    return p[0]*np.exp(-(p[2]/energy)**p[3])*np.power(energy/E0,p[1])*np.exp(-(energy/p[4])**p[5])


from scipy.special import erf

def Gaussian(energy,p):
    """
    Returns:
    Differential flux shaped as symmetric Gaussian "in units of ph/cm2/s/keV"
    
    Parameters:
    :param energy: 1D-array of energies where power-law is evaluated (in keV)
    :param p:      p[0] = Normalisation (area under Gaussian) (in ph/cm2/s)
    :param p:      p[1] = Line centroid (in keV)
    :param p:      p[2] = Line width (1-sigma value) (in keV)
    """    

    return p[0]/(np.sqrt(2*np.pi)*p[2])*np.exp(-(energy-p[1])**2/(2*p[2]**2))


def trunc_Gaussian(x,p):
    N0 = p[0]
    mu = p[1]
    sigma = p[2]
    lo = p[3]
    hi = p[4]
    val = np.zeros(len(x))
    idx = np.where((x >= lo) & (x <= hi))[0]
    phi = Gaussian(x[idx],[1,mu,sigma])
    theta_hi = 0.5*(1+erf((hi-mu)/(np.sqrt(2)*sigma)))
    theta_lo = 0.5*(1+erf((lo-mu)/(np.sqrt(2)*sigma)))
    val[idx] = N0*phi/((theta_hi-theta_lo))
    return(val)



def fit_powerlaw_Gauss(energy,p):
    plaw = powerlaw(energy,p[0:2])
    line = Gaussian(energy,[p[2],1275,8.5])
    return plaw+line


def ortho_posi_int(p0,E0):

    """
    Returns:
    Integrated flux of an ortho-positronium continuum "in units of ph/cm2/s"

    Parameters:
    :param energy: 1D-array of energies where power-law is evaluated (in keV)
    :param p:      p[0] = Amplitude (in ph/cm2/s/keV)
    :param p:      p[1] = red-shift (not working properly, only in non-relativistic case!)
    """
    return p0*E0*(np.pi**2 - 9)

def ortho_posi(energy,p):
    """
    Returns:
    Differential flux shaped as ortho-positronium continuum "in units of ph/cm2/s/keV"
    after Ore & Powell (1949)

    Parameters:
    :param energy: 1D-array of energies where power-law is evaluated (in keV)
    :param p:      p[0] = Normalisation (area under continuum) (in ph/cm2/s)
    :param p:      p[1] = red-shift (not working properly, only in non-relativistic case!)
    """

    m = p[1]

    t1 = (energy*(m-energy)) / (2*m-energy)**2
    t2 = (2*m*(m-energy)**2)/(2*m-energy)**3
    t3 = (m-energy)/m
    t4 = (2*m-energy)/energy
    t5 = 2*m*(m-energy)/energy**2
    t6 = (m-energy)/m

    val = p[0] / ortho_posi_int(1,m) * 2 * (t1 - t2*np.log(t3) + t4 + t5*np.log(t6))
    val[~np.isfinite(val)] = 0

    return val


def posi_fraction(I2,I3):
    """
    Returns:
    Positronium fraction given the 511 keV line flux (direct+pPs) and the oPs flux
    after Prantzos et al. (2011) and references therein

    Parameters:
    :param I2:     I2 = 511 keV line flux in (in ph/cm2/s)
    :param I3:     I3 = oPs flux (in ph/cm2/s)
    """
    return (8 * I3 / I2) / (9 + 6 * I3 / I2)



def posi_fraction_from_ratio(R32):
    """
    Returns:
    Positronium fraction given the 511 keV line flux (direct+pPs) and the oPs flux ratio
    after Prantzos et al. (2011) and references therein

    Parameters:
    :param R32:     R32 = 511 keV line flux to oPs flux ratio
    """
    return (8 * R32) / (9 + 6 * R32)

def ratio_from_posi_fraction(fps):
    """
    Returns:
    511 keV line flux (direct+pPs) and the oPs flux ratio for a given positronium fraction
    in in ph/cm2/s/keV
    similar to Jean et al. (2006)

    Parameters:
    :param fps:     fps = positronium fraction
    """
    return (9 * fps) / (8 - 6 * fps)




def AiF_thermal_pair_annihilation(energy,p):
    """
    Returns:
    Differential spectrumm of thermal pair annihilation with approximation formula of Svensson (1996)
    No Doppler-shifts included.

    Parameters:
    :param p:     p[0] = total(!) integrated flux (in ph/cm2/s)
    :param p:     p[1] = MJ temperature (in keV)
    """
    # e+- mass in units of keV
    mec2 = 510.998910 
    
    # reduced temperature
    theta = p[1] / (mec2)
    
    # reduced energu
    x = energy / (mec2)
    
    # temperature * energy
    y = x*theta
    
    # approximation coefficients (Svensson 1996)
    c1 = 6.8515487
    c2 = 1.4251694
    c3 = 0.017790149
    d1 = 4.63115589
    d2 = 1.5253007
    d3 = 0.04522338
    
    # rational approximation formula
    C = (1 + c1*y + c2*y**2 + c3*y**3)/(1 + d1*y + d2*y**2 + d3*y**3)

    K2 = 4 * theta**4 * np.exp(-2 / theta) * (1 + 2.0049 / theta + 1.4774 / theta**2 + np.pi / (2 * theta)**3 )

    # unnormalised flux
    tmp = (np.pi / theta)**(0.5) * x**(3 / 2) * np.exp(- (x + x**(-1)) / theta) * C / K2

    # flux normalisation
    norm = trapz(tmp,x)
    
    # normalised flux
    val = p[0]*tmp/norm
    
    return(val)



def internal_bremsstrahlung_spectrum(x,p):
    alpha = 1./137
    s = 4*mdm
    sprime = 4*mdm*(mdm-x)
    me = 511.
    return alpha/np.pi * 1/x * (np.log(sprime/me**2) - 1) * (1 + (sprime/s)**2)


def fit_positronium_at511(energy,p):
    """
    Returns:
    Differential spectrumm of 511 keV line plus the oPs continuum for given positronium fraction
    after Prantzos et al. (2011) and references therein
    Note that here, no Dopper-shifts or -broadenings are allowed and lab (and astro) values are
    used instead

    Parameters:
    :param p:     p[0] = absolute normalisation (chosen as line flux)
    :param p:     p[1] = positronium fraction
    """
    line = Gaussian(energy,[1.,511.,1.4])
    oPs = ortho_posi(energy,[1.,511.])
    ratio = ratio_from_posi_fraction(p[1])
    
    return p[0]*(line + ratio*oPs)


def fit_IGspec511low(energy,p):
    """
    Returns:
    Differential spectrumm of 511 keV line plus the oPs continuum for given positronium fraction
    plus a general power-law for diffuse emission and a cut-off powerlaw describing unresolved
    point sources
    Note that here, no Dopper-shifts or -broadenings are allowed and lab (and astro) values are
    used instead
     
    Parameters:
    :param p:     p[0] = absolute normalisation (chosen as line flux)
    :param p:     p[1] = positronium fraction
    :param p:     p[2] = normalisation of overall power-law
    :param p:     p[3] = power-law index of overall power-law
    :param p:     p[4] = normalisation of cutoff power-law
    :param p:     p[5] = power-law index of cutoff power-law
    :param p:     p[6] = cutoff energy (average of all unresolved point sources)
    """

    five11 = fit_positronium_at511(energy,[p[0],p[1]])
    plaw = powerlaw(energy,[p[2],p[3]])
    cutoffpl = cutoff_powerlaw(energy,[p[4],p[5],p[6]])

    return cutoffpl + plaw + five11


def fit_IGspec511_highEplaw_26Al(energy,p):
    """
    Returns:
    Differential spectrumm of 511 keV line plus the oPs continuum for given positronium fraction
    plus a general power-law for diffuse emission and a cut-off powerlaw describing unresolved
    point sources
    Note that here, no Dopper-shifts or -broadenings are allowed and lab (and astro) values are
    used instead
     
    Parameters:
    :param p:     p[0] = absolute normalisation (chosen as line flux)
    :param p:     p[1] = positronium fraction
    :param p:     p[2] = normalisation of overall power-law
    :param p:     p[3] = power-law index of overall power-law
    :param p:     p[4] = normalisation of cutoff power-law
    :param p:     p[5] = power-law index of cutoff power-law
    :param p:     p[6] = cutoff energy (average of all unresolved point sources)
    :param p:     p[7] = flux of 26Al line
    """

    five11 = fit_positronium_at511(energy,[p[0],p[1]])
    plaw = powerlaw(energy,[p[2],p[3]])
    cutoffpl = cutoff_powerlaw(energy,[p[4],p[5],p[6]])
    al26 = Gaussian(energy,[p[7],1808.74,1.36])
    
    return cutoffpl + plaw + five11 + al26


def fit_IGspec511_highEplaw_26Al_60Fe(energy,p):
    """
    Returns:
    Differential spectrumm of 511 keV line plus the oPs continuum for given positronium fraction
    plus a general power-law for diffuse emission and a cut-off powerlaw describing unresolved
    point sources
    Note that here, no Dopper-shifts or -broadenings are allowed and lab (and astro) values are
    used instead
     
    Parameters:
    :param p:     p[0] = absolute normalisation (chosen as line flux)
    :param p:     p[1] = positronium fraction
    :param p:     p[2] = normalisation of overall power-law
    :param p:     p[3] = power-law index of overall power-law
    :param p:     p[4] = normalisation of cutoff power-law
    :param p:     p[5] = power-law index of cutoff power-law
    :param p:     p[6] = cutoff energy (average of all unresolved point sources)
    :param p:     p[7] = flux of 26Al line
    :param p:     p[8] = ratio of 60Fe/26Al (to set normalisation of 60Fe lines)
    """

    five11 = fit_positronium_at511(energy,[p[0],p[1]])
    plaw = powerlaw(energy,[p[2],p[3]])
    cutoffpl = cutoff_powerlaw(energy,[p[4],p[5],p[6]])
    al26 = Gaussian(energy,[p[7],1808.74,1.36])
    fe60 = p[8]*(Gaussian(energy,[p[7],1172.5,1.15]) + Gaussian(energy,[p[7],1332.5,1.19]))
    
    return cutoffpl + plaw + five11 + al26 + fe60


def fit_IGspec511_highEplaw_26Al_60Fe_AiF(energy,p):
    """
    Returns:
    Differential spectrumm of 511 keV line plus the oPs continuum for given positronium fraction
    plus a general power-law for diffuse emission and a cut-off powerlaw describing unresolved
    point sources
    Note that here, no Dopper-shifts or -broadenings are allowed and lab (and astro) values are
    used instead
     
    Parameters:
    :param p:     p[0] = absolute normalisation (chosen as line flux)
    :param p:     p[1] = positronium fraction
    :param p:     p[2] = normalisation of overall power-law
    :param p:     p[3] = power-law index of overall power-law
    :param p:     p[4] = normalisation of cutoff power-law
    :param p:     p[5] = power-law index of cutoff power-law
    :param p:     p[6] = cutoff energy (average of all unresolved point sources)
    :param p:     p[7] = flux of 26Al line
    :param p:     p[8] = ratio of 60Fe/26Al (to set normalisation of 60Fe lines)
    :param p:     p[9] = total(!) flux of AiF component (integrated from 511/2 keV to infinity)
    :param p:    p[10] = lg(injection energy) (temperature) of positron population in lg(MeV)
    """

    five11 = fit_positronium_at511(energy,[p[0],p[1]])
    plaw = powerlaw(energy,[p[2],p[3]])
    cutoffpl = cutoff_powerlaw(energy,[p[4],p[5],p[6]])
    al26 = Gaussian(energy,[p[7],1808.74,1.36])
    fe60 = p[8]*(Gaussian(energy,[p[7],1172.5,1.15]) + Gaussian(energy,[p[7],1332.5,1.19]))
    aif = thermal_cosmic_ray_annihilation_spectrum(energy/511,[p[9],10**p[10]])
    #aif = thermal_cosmic_ray_annihilation_spectrum(energy/511,[p[9],p[10]])
    
    return cutoffpl + plaw + five11 + al26 + fe60 + aif


def fit_IGspec511_highEplaw_26Al_60Fe_AiFBY06(energy,p):
    """
    Returns:
    Differential spectrumm of 511 keV line plus the oPs continuum for given positronium fraction
    plus a general power-law for diffuse emission and a cut-off powerlaw describing unresolved
    point sources
    Note that here, no Dopper-shifts or -broadenings are allowed and lab (and astro) values are
    used instead
     
    Parameters:
    :param p:     p[0] = absolute normalisation (chosen as line flux)
    :param p:     p[1] = positronium fraction
    :param p:     p[2] = normalisation of overall power-law
    :param p:     p[3] = power-law index of overall power-law
    :param p:     p[4] = normalisation of cutoff power-law
    :param p:     p[5] = power-law index of cutoff power-law
    :param p:     p[6] = cutoff energy (average of all unresolved point sources)
    :param p:     p[7] = flux of 26Al line
    :param p:     p[8] = ratio of 60Fe/26Al (to set normalisation of 60Fe lines)
    :param p:     p[9] = annihilation probabiliy from Eq. (5) in BY06
    :param p:    p[10] = lg(injection energy) (temperature) of positron population in lg(MeV)
    """

    five11 = fit_positronium_at511(energy,[p[0],p[1]])
    plaw = powerlaw(energy,[p[2],p[3]])
    cutoffpl = cutoff_powerlaw(energy,[p[4],p[5],p[6]])
    al26 = Gaussian(energy,[p[7],1808.74,1.36])
    fe60 = p[8]*(Gaussian(energy,[p[7],1172.5,1.15]) + Gaussian(energy,[p[7],1332.5,1.19]))
    
    # BY06 aif flux normalisation with Ps flux and annihilation probability:
    aif_norm = (1./(1-3.*p[1]/4.))*(1-p[9])/p[9]*p[0]
    aif = thermal_cosmic_ray_annihilation_spectrum(energy/511,[aif_norm,10**p[10]])
    #aif = thermal_cosmic_ray_annihilation_spectrum(energy/511,[p[9],p[10]])
    
    return cutoffpl + plaw + five11 + al26 + fe60 + aif


def fit_IGspec511_highEplaw_26Al_60Fe_CR(energy,p):
    """
    Returns:
    Differential spectrumm of 511 keV line plus the oPs continuum for given positronium fraction
    plus a general power-law for diffuse emission and a cut-off powerlaw describing unresolved
    point sources
    Note that here, no Dopper-shifts or -broadenings are allowed and lab (and astro) values are
    used instead
     
    Parameters:
    :param p:     p[0] = absolute normalisation (chosen as line flux)
    :param p:     p[1] = positronium fraction
    :param p:     p[2] = normalisation of overall power-law
    :param p:     p[3] = power-law index of overall power-law
    :param p:     p[4] = normalisation of cutoff power-law
    :param p:     p[5] = power-law index of cutoff power-law
    :param p:     p[6] = cutoff energy (average of all unresolved point sources)
    :param p:     p[7] = flux of 26Al line
    :param p:     p[8] = ratio of 60Fe/26Al (to set normalisation of 60Fe lines)
    :param p:     p[9] = normalisation of Benhabiles-Mezhoud+2013 CR excitation spectrum
    """

    five11 = fit_positronium_at511(energy,[p[0],p[1]])
    plaw = powerlaw(energy,[p[2],p[3]])
    cutoffpl = cutoff_powerlaw(energy,[p[4],p[5],p[6]])
    al26 = Gaussian(energy,[p[7],1808.74,1.36])
    fe60 = p[8]*(Gaussian(energy,[p[7],1172.5,1.15]) + Gaussian(energy,[p[7],1332.5,1.19]))
    cr = lecr_spectrum(energy,[p[9]])

    return cutoffpl + plaw + five11 + al26 + fe60 + cr


def fit_IGspec511_highEplaw_26Al_60Fe_TPA(energy,p):
    """
    Returns:
    Differential spectrumm of 511 keV line plus the oPs continuum for given positronium fraction
    plus a general power-law for diffuse emission and a cut-off powerlaw describing unresolved
    point sources
    Note that here, no Dopper-shifts or -broadenings are allowed and lab (and astro) values are
    used instead
     
    Parameters:
    :param p:     p[0] = absolute normalisation (chosen as line flux)
    :param p:     p[1] = positronium fraction
    :param p:     p[2] = normalisation of overall power-law
    :param p:     p[3] = power-law index of overall power-law
    :param p:     p[4] = normalisation of cutoff power-law
    :param p:     p[5] = power-law index of cutoff power-law
    :param p:     p[6] = cutoff energy (average of all unresolved point sources)
    :param p:     p[7] = flux of 26Al line
    :param p:     p[8] = ratio of 60Fe/26Al (to set normalisation of 60Fe lines)
    :param p:     p[9] = normalisation of TPA spectrum
    :param p:    p[10] = temperature (in units of 511 keV) of TPA spectrum
    """

    five11 = fit_positronium_at511(energy,[p[0],p[1]])
    plaw = powerlaw(energy,[p[2],p[3]])
    cutoffpl = cutoff_powerlaw(energy,[p[4],p[5],p[6]])
    al26 = Gaussian(energy,[p[7],1808.74,1.36])
    fe60 = p[8]*(Gaussian(energy,[p[7],1172.5,1.15]) + Gaussian(energy,[p[7],1332.5,1.19]))
    tpa = AiF_thermal_pair_annihilation(energy,[p[9],p[10]])

    return cutoffpl + plaw + five11 + al26 + fe60 + tpa




def fit_IGspec511high(energy,p):
    """
    Returns:
    Differential spectrumm of 511 keV line plus the oPs continuum for given positronium fraction
    plus a general power-law for diffuse emission and a cut-off powerlaw describing unresolved
    point sources
    Note that here, no Dopper-shifts or -broadenings are allowed and lab (and astro) values are
    used instead
     
    Parameters:
    :param p:     p[0] = absolute normalisation (chosen as line flux)
    :param p:     p[1] = positronium fraction
    :param p:     p[2] = normalisation of overall power-law
    :param p:     p[3] = power-law index of overall power-law
    :param p:     p[4] = normalisation of cutoff power-law
    :param p:     p[5] = power-law index of cutoff power-law
    :param p:     p[6] = cutoff energy (average of all unresolved point sources)
    :param p:     p[7] = total (!) flux of AiF spectrum
    :param p:     p[8] = MJ-distribution temperature
    """

    five11 = fit_positronium_at511(energy,[p[0],p[1]])
    plaw = powerlaw(energy,[p[2],p[3]])
    cutoffpl = cutoff_powerlaw(energy,[p[4],p[5],p[6]])
    aif = thermal_cosmic_ray_annihilation_spectrum(energy/511,[p[7],p[8]])

    return cutoffpl + plaw + five11 + aif



def fit_IGspec511_nucsyn(energy,p):
    """
    Returns:
    Differential spectrumm of 511 keV line plus the oPs continuum for given positronium fraction
    plus a general power-law for diffuse emission and a cut-off powerlaw describing unresolved
    point sources
    Note that here, no Dopper-shifts or -broadenings are allowed and lab (and astro) values are
    used instead
     
    Parameters:
    :param p:     p[0] = absolute normalisation (chosen as line flux)
    :param p:     p[1] = positronium fraction
    :param p:     p[2] = normalisation of overall power-law
    :param p:     p[3] = power-law index of overall power-law
    :param p:     p[4] = normalisation of cutoff power-law
    :param p:     p[5] = power-law index of cutoff power-law
    :param p:     p[6] = cutoff energy (average of all unresolved point sources)
    :param p:     p[7] = flux of 26Al line
    :param p:     p[8] = ratio of 60Fe/26Al (to set normalisation of 60Fe lines)
    :param p:     p[9] = flux of broadened 22Na line
    """

    five11 = fit_positronium_at511(energy,[p[0],p[1]])
    plaw = powerlaw(energy,[p[2],p[3]])
    cutoffpl = cutoff_powerlaw(energy,[p[4],p[5],p[6]])
    #aif = thermal_cosmic_ray_annihilation_spectrum(energy/511,[p[7],p[8]])
    al26 = Gaussian(energy,[p[7],1808.74,1.36])
    fe60 = p[8]*(Gaussian(energy,[p[7],1172.5,1.15]) + Gaussian(energy,[p[7],1332.5,1.19]))
    na22 = Gaussian(energy,[p[9],1275.0,8.41])

    return cutoffpl + plaw + five11 + al26 + fe60 + na22



def fit_IGspec511high_nucsyn(energy,p):
    """
    Returns:
    Differential spectrumm of 511 keV line plus the oPs continuum for given positronium fraction
    plus a general power-law for diffuse emission and a cut-off powerlaw describing unresolved
    point sources
    Note that here, no Dopper-shifts or -broadenings are allowed and lab (and astro) values are
    used instead
     
    Parameters:
    :param p:     p[0] = absolute normalisation (chosen as line flux)
    :param p:     p[1] = positronium fraction
    :param p:     p[2] = normalisation of overall power-law
    :param p:     p[3] = power-law index of overall power-law
    :param p:     p[4] = normalisation of cutoff power-law
    :param p:     p[5] = power-law index of cutoff power-law
    :param p:     p[6] = cutoff energy (average of all unresolved point sources)
    :param p:     p[7] = total (!) flux of AiF spectrum
    :param p:     p[8] = MJ-distribution temperature
    :param p:     p[9] = flux of 26Al line
    :param p:     p[10] = ratio of 60Fe/26Al (to set normalisation of 60Fe lines)
    :param p:     p[11] = flux of broadened 22Na line
    """

    five11 = fit_positronium_at511(energy,[p[0],p[1]])
    plaw = powerlaw(energy,[p[2],p[3]])
    cutoffpl = cutoff_powerlaw(energy,[p[4],p[5],p[6]])
    aif = thermal_cosmic_ray_annihilation_spectrum(energy/511,[10**p[7],10**p[8]])
    al26 = Gaussian(energy,[p[9],1808.74,1.36])
    fe60 = p[10]*(Gaussian(energy,[p[9],1172.5,1.15]) + Gaussian(energy,[p[9],1332.5,1.19]))
    na22 = Gaussian(energy,[p[11],1275.0,8.41])

    return cutoffpl + plaw + five11 + aif + al26 + fe60 + na22



def fit_IGspec511high_nucsyn_CR(energy,p):
    """
    Returns:
    Differential spectrumm of 511 keV line plus the oPs continuum for given positronium fraction
    plus a general power-law for diffuse emission and a cut-off powerlaw describing unresolved
    point sources
    Note that here, no Dopper-shifts or -broadenings are allowed and lab (and astro) values are
    used instead
     
    Parameters:
    :param p:     p[0] = absolute normalisation (chosen as line flux)
    :param p:     p[1] = positronium fraction
    :param p:     p[2] = normalisation of overall power-law
    :param p:     p[3] = power-law index of overall power-law
    :param p:     p[4] = normalisation of cutoff power-law
    :param p:     p[5] = power-law index of cutoff power-law
    :param p:     p[6] = cutoff energy (average of all unresolved point sources)
    :param p:     p[7] = total (!) flux of AiF spectrum
    :param p:     p[8] = MJ-distribution temperature
    :param p:     p[9] = flux of 26Al line
    :param p:     p[10] = ratio of 60Fe/26Al (to set normalisation of 60Fe lines)
    :param p:     p[11] = flux of broadened 22Na line
    :param p:     p[12] = normalisation of Benhabiles-Mezhoud+2013 CR excitation spectrum
    """

    five11 = fit_positronium_at511(energy,[p[0],p[1]])
    plaw = powerlaw(energy,[p[2],p[3]])
    cutoffpl = cutoff_powerlaw(energy,[p[4],p[5],p[6]])
    aif = thermal_cosmic_ray_annihilation_spectrum(energy/511,[10**p[7],10**p[8]])
    al26 = Gaussian(energy,[p[9],1808.74,1.36])
    fe60 = p[10]*(Gaussian(energy,[p[9],1172.5,1.15]) + Gaussian(energy,[p[9],1332.5,1.19]))
    na22 = Gaussian(energy,[p[11],1275.0,8.41])
    cr = lecr_spectrum(energy,[p[12]])

    return cutoffpl + plaw + five11 + aif + al26 + fe60 + na22 + cr



# TS: made a vectorised version out of that
def diff_spec_2rel_v2(k, energy_pos, ne):
    k = k/511.
    gam_plus = energy_pos/511.
    mom_plus = np.sqrt((gam_plus**2)-1)

    emin = 0.5*(gam_plus+1-mom_plus)
    emax = 0.5*(gam_plus+1+mom_plus)

    a = np.pi*(r0**2)*(3E10)*ne/(gam_plus*mom_plus)
    b = k/(2*gam_plus-k)
    d = (2*gam_plus-k)/k
    f = 2*((1/k)+((1)/(2*gam_plus-k)))
    g = ((1/k)+((1)/(2*gam_plus-k)))
    h = a*((b+d)+f-g)
    
    h[np.where((emin > k) | (emax < k))[0]] = 0
    
    return h



from scipy.special import kn as Bessel2
from scipy.special import expi as E1
from scipy.integrate import quad as quad

def thermal_cosmic_ray_flux(p):
    """
    Cosmic-ray annihilation spectrum flux of positrons (electrons) with (relativistic) temperature theta
    (in units of 511 keV) onto electrons (positrons) at rest. Following Svensson 1982, Eq. (65)
    Calculates area under curve as normalisation to the spectrum
#    :param: k     Energy array in units of 511 keV
    :param: p[0]  Total flux (normalisation of spectrum, area under the curve)
    :param: p[1]  Temperature of Maxwell-Juettner-distribution of particles *not* at rest (theta)
    """
    
    return (quad(lambda x: thermal_cosmic_ray_annihilation_spectrum_all(x,p),0.5,np.inf))[0]*511



def thermal_cosmic_ray_annihilation_spectrum_all(k,p):
    """
    Cosmic-ray annihilation spectrum of positrons (electrons) with (relativistic) temperature theta
    (in units of 511 keV) onto electrons (positrons) at rest. Following Svensson 1982, Eq. (65)
    :param: k     Energy array in units of 511 keV
    :param: p[0]  Total flux (normalisation of spectrum, area under the curve)
    :param: p[1]  Temperature of Maxwell-Juettner-distribution of particles *not* at rest (theta)
    """
    
    # temperature
    theta = p[1]
    
    # classical electron radius
    re = 2.82e-13
    # speed of light
    c  = 3e10
    # classical annihilation cross section
    ann0 = np.pi * re**2 * c

    # often occurring energy term
    k21 = 2*k-1
    # from MJ-dist
    K2 = Bessel2(2,1/theta)
    # amplitude of Eq.(65) in Svensson 1982
    amp = ann0/(k*K2*theta**2)
    # Eq. (36)
    gamma_B = (2*k**2-k21)/k21
    
    # individual terms of the final formula
    exp_term = np.exp(-gamma_B/theta)
    t1 = theta**3
    t2 = (k/k21 + 2 - 1/k)*theta**2
    t3 = k21*theta
    omega = k/(theta*k21)
    #print(omega)
    #print(np.exp(omega))
    #print(E1(omega))
    #t4 = np.exp(omega)*E1(omega)
    # TS: so the problem is that exp(omega) and E1(omega) can be very large and very small, but multiplied
    #     together, like 1e1000*1e-999 equals 1e-1, so a reasonable number. The fact that there is no easy
    #     approximation for ln(E1(omega)) makes it difficult to evaluate for arbitrary things, so that we 
    #     need to introduce (would have done that anyway) a scaling parameter, as I think the low temperature
    #     values get to many annihilations than the high-energy ones
    # for real arguments, E1 behaves like exp(x)+log(abs(x))
    # so that t4 = np.exp(omega)*E1(omega) = np.exp(omega)*(np.exp(omega)+np.log(np.abs(omega))), thus
    # ln(t4) = omega + np.log(E1(omega)) w
    #t4 = np.exp(omega+np.log(E1(omega)))
    # TS: well doesn't work either; best approximation:
    # np.exp(omega)*E1(omega) ~ 1 ...: overestimates for small theta? no shoulder at high energy?
    t4 = 2.
    t5 = (k**2 + k21 - 1)*theta + k
    
    # everything together
    val = amp*exp_term*(t1+t2-t3+t4*t5)
    
    # cutoff at 511/2 keV
    #val[np.where(k <= 0.5)[0]] = 0
    
    return val



def thermal_cosmic_ray_annihilation_spectrum(k,p):
    """
    Cosmic-ray annihilation spectrum of positrons (electrons) with (relativistic) temperature theta
    (in units of 511 keV) onto electrons (positrons) at rest. Following Svensson 1982, Eq. (65)
    :param: k     Energy array in units of 511 keV
    :param: p[0]  Total flux (normalisation of spectrum, area under the curve)
    :param: p[1]  Temperature of Maxwell-Juettner-distribution of particles *not* at rest (theta)
    """
    
    spec = thermal_cosmic_ray_annihilation_spectrum_all(k,p)
    
    norm = thermal_cosmic_ray_flux(p)
    
    # everything together
    val = p[0]*spec/norm
    
    # cutoff at 511/2 keV
    val[np.where(k <= 0.5)[0]] = 0
    
    return val




def MB_dist(gamma,theta):
    beta = np.sqrt(1 - 1 / gamma**2)
    val = (theta*Bessel2(2,1/theta))**(-1)*beta*gamma**2*np.exp(-gamma/theta)
    return val




def gamma2beta(gamma):
    return np.sqrt(1 - 1/gamma**2)



def ann_spec(x, p):
    
    n_x = len(x)
    k  = x
    
    gp = p[0]
    gm = p[1]
    
    bp = gamma2beta(gp)
    bm = gamma2beta(gm)
    
    re = 2.82e-13
    c  = 3e10
    ann0 = np.pi * re**2 * c

    n = bp * gp**2 * bm * gm**2

    ku = 0.5 * (gp*(1 + bp) + gm*(1 + bm))
    kp = 0.5 * (gp*(1 + bp) + gm*(1 - bm))
    km = 0.5 * (gp*(1 - bp) + gm*(1 + bm))
    kl = 0.5 * (gp*(1 - bp) + gm*(1 - bm))
    kt = 0.5 * (gp + gm)
  
    gcm_star = np.sqrt(k*(gp+gm-k))

    gcm_max = np.sqrt(0.5 * (1 + gp*gm + gp*gm*bp*bm))
    gcm_min = np.sqrt(0.5 * (1 + gp*gm - gp*gm*bp*bm))

    gcmu = np.minimum(np.repeat(gcm_max,n_x),gcm_star)

    gcml = np.repeat(gcm_min,n_x)

    val = np.zeros(n_x)

    if ((gp > 1) & (gm > 1)):
        for i in range(n_x):
            pprimeu = np.array(list(p)+[gcmu[i]])
            pprimel = np.array(list(p)+[gcml[i]])
            val[i] = ann0 / n * (sven55(x[i],pprimeu) - sven55(x[i],pprimel))
    else:
        gr = max([gp,gm])
        br = gamma2beta(gr)
        
        val = ann0 / (br*gr**2) * ((-(3 + gr)/(1 + gr) + (3 + gr)/k - 1 / k**2)/(1 - k/(1 + gr))**2 - 2)

    out = np.where((k < kl) | (k > ku))

    val[out[0]] = 0.
    
#    norm = integrate(val,x=x,dx=np.diff(x),)
    

    return val


def sven55(x, p):
    
    k  = x
    gp = p[0]
    gm = p[1]
    gcm = p[2]
    bp = gamma2beta(gp)
    bm = gamma2beta(gm)

    re = 2.82e-13
    c  = 3e10
    ann0 = np.pi * re**2 * c
    
    n = bp * gp**2 * bm * gm**2

    ku = 0.5 * (gp*(1 + bp) + gm*(1 + bm))
    kp = 0.5 * (gp*(1 + bp) + gm*(1 - bm))
    km = 0.5 * (gp*(1 - bp) + gm*(1 + bm))
    kl = 0.5 * (gp*(1 - bp) + gm*(1 - bm))
    kt = 0.5 * (gp + gm)

    gcm_star = np.sqrt(k*(gp+gm-k))

    gcm_max = np.sqrt(0.5 * (1 + gp*gm + gp*gm*bp*bm))
    gcm_min = np.sqrt(0.5 * (1 + gp*gm - gp*gm*bp*bm))

    dp = gm*(gp + gm) + k*(gp - gm)
    dm = gp*(gp + gm) - k*(gp - gm)

    cp = (gm - k)**2 - 1
    cm = (gp - k)**2 - 1

    up = np.sqrt(cp*gcm**2 + gcm_star**2)
    um = np.sqrt(cm*gcm**2 + gcm_star**2)

#  ;sven57
  
    if (cp > 0):
        Ip = 1 / np.sqrt(cp) * np.log(gcm * np.sqrt(cp) + up)
    else:
        Ip = 1 / np.sqrt(-cp) * np.arcsin(gcm/gcm_star*np.sqrt(-cp))

    if (cm > 0):
        Im = 1 / np.sqrt(cm) * np.log(gcm * np.sqrt(cm) + um)
    else:
        Im = 1 / np.sqrt(-cm) * np.arcsin(gcm/gcm_star*np.sqrt(-cm))

    if (cp == 0):
        Hp = (2 / 3 * gcm**3 + 3 * gcm + 1 / gcm)/gcm_star + 0.5 * (2 / 3 * gcm**3 - dp*gcm)/gcm_star**3
    else:
        Hp = (2 + (1 - gcm_star**2)/cp)*Ip + (1/gcm - gcm/cp + gcm/(2 * gcm_star**2)*(2 * cp - dp))/up + gcm/cp*up

    if (cm == 0):
        Hm = (2 / 3 * gcm**3 + 3 * gcm + 1 / gcm)/gcm_star + 0.5 * (2 / 3 * gcm**3 - dm*gcm)/gcm_star**3
    else:
        Hm = (2 + (1 - gcm_star**2)/cm)*Im + (1/gcm - gcm/cm + gcm/(2 * gcm_star**2)*(2 * cm - dm))/um + gcm/cm*um

    val = np.sqrt((gp+gm)**2 - 4 * gcm**2) + Hp + Hm
    
    return val



def ann_spec_flux_norm(gm,gp,p):
    GM, GP = np.meshgrid(gm,gp)
    return p[0]*GM**p[1]*np.exp(-p[2]/GM)*GP**p[1]*np.exp(-p[2]/GP)*(GM+GP)**p[3]*np.exp(-p[4]/(GM+GP))



def ann_spec_with_norm(x, p):
    
    # parameters from fitted result.x
    q = np.array([20.11990836e-12, -0.93042488,  0.1984951 ,  0.08272161,  1.12874467])
    
    n_x = len(x)
    k  = x
    
    gp = p[0]
    gm = p[1]
    
    bp = gamma2beta(gp)
    bm = gamma2beta(gm)
    
    re = 2.82e-13
    c  = 3e10
    ann0 = np.pi * re**2 * c

    n = bp * gp**2 * bm * gm**2

    ku = 0.5 * (gp*(1 + bp) + gm*(1 + bm))
    kp = 0.5 * (gp*(1 + bp) + gm*(1 - bm))
    km = 0.5 * (gp*(1 - bp) + gm*(1 + bm))
    kl = 0.5 * (gp*(1 - bp) + gm*(1 - bm))
    kt = 0.5 * (gp + gm)
  
    gcm_star = np.sqrt(k*(gp+gm-k))

    gcm_max = np.sqrt(0.5 * (1 + gp*gm + gp*gm*bp*bm))
    gcm_min = np.sqrt(0.5 * (1 + gp*gm - gp*gm*bp*bm))

    gcmu = np.minimum(np.repeat(gcm_max,n_x),gcm_star)

    gcml = np.repeat(gcm_min,n_x)

    val = np.zeros(n_x)

    if ((gp > 1) & (gm > 1)):
        for i in range(n_x):
            pprimeu = np.array(list(p[0:2])+[gcmu[i]])
            pprimel = np.array(list(p[0:2])+[gcml[i]])
            val[i] = ann0 / n * (sven55(x[i],pprimeu) - sven55(x[i],pprimel))
    else:
        gr = max([gp,gm])
        br = gamma2beta(gr)
        
        val = ann0 / (br*gr**2) * ((-(3 + gr)/(1 + gr) + (3 + gr)/k - 1 / k**2)/(1 - k/(1 + gr))**2 - 2)

    out = np.where((k < kl) | (k > ku))

    val[out[0]] = 0.
    
    norm = ann_spec_flux_norm(gp,gm,q)
    
    val *= p[2]/norm[0]

    return val



lecr = pd.read_csv('gam_lecr_pdcompatible.spc',skiprows=3,sep=',')
e_cr = lecr['Eg(MeV)']
f_cr = lecr['Fg(10^-5 ph/cm^2/s/sr/MeV)']
sr_norm = (np.sin(np.deg2rad(25))-np.sin(np.deg2rad(-25)))*np.deg2rad(50)
from scipy.interpolate import interp1d as interpol
lecr_interpol = interpol(e_cr*1000,f_cr/1000*sr_norm*1e-5,bounds_error=False,fill_value=0)

def lecr_spectrum(energy,p):
    return p*lecr_interpol(energy)



"""
Created on Wed Sep  6 15:41:22 2017

@author: mpleinti

Collection of different secondary functions 

adapted by tsiegert from mpleinti from tsiegert
"""

from scipy.special import erfc
from scipy.optimize import newton
from scipy.optimize import fmin

###############################################################################
### TOOLS FOR DEG_GAUSSFIT ####################################################
###############################################################################

def ln_erfc(x):
    """ logarithmic approximation of the complementary errorfunction"""
    x = np.array(x)
    
    a = [-1.26551223, 1.00002368, 0.37409196, 0.09678418, -0.18628806, 
         0.27886807, -1.13520398, 1.48851587, -0.82215223, 0.17087277]
      
    t = 1.0 / (1.0 + 0.5*np.abs(x))
    
    tau = t * np.exp( -x*x + (a[0] + t*(a[1] + t*(a[2] + t*(a[3] + t*(a[4] +\
                    t*(a[5] + t*(a[6] + t*(a[7] + t*(a[8] + t*a[9]))))))))))
    
    y = np.log(t) + ( -x*x + (a[0] + t*(a[1] + t*(a[2] + t*(a[3] + t*\
            (a[4] + t*(a[5] + t*(a[6] + t*(a[7] + t*(a[8] + t*a[9]))))))))))
    
    lt0 = np.where(x < 0)[0]
    if len(lt0) > 0:
        y[lt0] = y[lt0] + np.log(2./tau[lt0] - 1.)
    
    return(y)

def conv_line_shape(x, amp, mu, sig, tau, ln=True):
    """
    Calculate differential spectrum of a degraded (asymmetric) Gaussian line
    shape
    """
    x = np.array(x)
    
    #if tau < 1e-3:                                  # if tau too small
    #    val = amp*np.exp(-(x - mu)**2./(2.*sig**2)) # use standard Gaussian
    
    #else:
        
    if ln:                                      # logarithmic approximation
        q = np.log(np.sqrt(np.pi/2.)*amp*sig/tau)
        w = (2.*tau*(x - mu) + sig**2)/(2.*tau**2)  # logarithmic term
        e = ln_erfc((tau*(x - mu) + sig**2)/(np.sqrt(2.)*sig*tau))
        edx = np.where(np.isfinite(e) == 0)[0]     # set inf values to first
                                                # non-inf value
        e[edx] = 0#e[np.where(np.isfinite(e) == 1)[0][0]]

        val = np.exp(q + w + e)


    else:                                       # standard degraded Gauss
        val = np.sqrt(np.pi/2.)*amp*sig*np.exp((2.*tau*(x - mu) +\
              sig**2)/(2.*tau**2))*(erfc((tau*(x - mu) +\
              sig**2)/(np.sqrt(2.)*sig**2*tau)))/tau
        
    return(val)

def conv_line_shape_nosig(x, amp, mu, tau, ln=True):
    """
    Calculate differential spectrum of a degraded (asymmetric) Gaussian line
    shape without sigma
    """
    x = np.array(x)
    
    if tau < 1e-3:                                  # if tau too small
        val = amp*np.exp(-(x - mu)**2./(2.*1.34**2)) # use standard Gaussian
    
    else:
        
        if ln:                                      # logarithmic approximation
            q = np.log(np.sqrt(np.pi/2.)*amp*1.34/tau)
            w = (2.*tau*(x - mu) + 1.34**2)/(2.*tau**2)  # logarithmic term
            e = ln_erfc((tau*(x - mu) + 1.34**2)/(np.sqrt(2.)*1.34*tau))
            edx = np.where(np.isfinite(e) == 0)     # set inf values to first
                                                    # non-inf value
            e[edx] = e[np.where(np.isfinite(e) == 1)[0][0]]

            val = np.exp(q + w + e)
            
    
        else:                                       # standard degraded Gauss
            val = np.sqrt(np.pi/2.)*amp*1.34*np.exp((2.*tau*(x - mu) +\
                  1.34**2)/(2.*tau**2))*(erfc((tau*(x - mu) +\
                  1.34**2)/(np.sqrt(2.)*1.34**2*tau)))/tau
        
    return(val)

def conv_line_shape_cont(x, inc, slp, amp, mu, sig, tau, ln=True):
    """
    Calculate differential spectrum of a degraded (asymmetric) Gaussian line
    shape
    """
    x = np.array(x)
    
    if tau < 1e-3:                                  # if tau too small
        val = inc + slp*x + amp*np.exp(-(x - mu)**2./(2.*sig**2)) # standard
        
    else:
        
        if ln:                                      # logarithmic approximation
            q = np.log(np.sqrt(np.pi/2.)*amp*sig/tau)
            w = (2.*tau*(x - mu) + sig**2)/(2.*tau**2)  # logarithmic term
            e = ln_erfc((tau*(x - mu) + sig**2)/(np.sqrt(2.)*sig*tau))
            edx = np.where(np.isfinite(e) == 0)     # set inf values to first
                                                    # non-inf value
            e[edx] = e[np.where(np.isfinite(e) == 1)[0][0]]
 
            val = inc + slp*x + np.exp(q + w + e)
            
    
        else:                                       # standard degraded Gauss
            val = inc + slp*x +\
                  np.sqrt(np.pi/2.)*amp*sig*np.exp((2.*tau*(x - mu) +\
                  sig**2)/(2.*tau**2))*(erfc((tau*(x - mu) +\
                  sig**2)/(np.sqrt(2.)*sig**2*tau)))/tau
        
    return(val)

def conv_line_shape_cont_nosig(x, inc, slp, amp, mu, tau, ln=True):
    """
    Calculate differential spectrum of a degraded (asymmetric) Gaussian line
    shape
    """
    x = np.array(x)
    
    if tau < 1e-3:                                  # if tau too small
        val = inc + slp*x + amp*np.exp(-(x - mu)**2./(2.*1.34**2)) # standard
        
    else:
        
        if ln:                                      # logarithmic approximation
            q = np.log(np.sqrt(np.pi/2.)*amp*1.34/tau)
            w = (2.*tau*(x - mu) + 1.34**2)/(2.*tau**2)  # logarithmic term
            e = ln_erfc((tau*(x - mu) + 1.34**2)/(np.sqrt(2.)*1.34*tau))
            edx = np.where(np.isfinite(e) == 0)     # set inf values to first
                                                    # non-inf value
            e[edx] = e[np.where(np.isfinite(e) == 1)[0][0]]

            val = inc + slp*x + np.exp(q + w + e)
            
    
        else:                                       # standard degraded Gauss
            val = inc + slp*x +\
                  np.sqrt(np.pi/2.)*amp*1.34*np.exp((2.*tau*(x - mu) +\
                  1.34**2)/(2.*tau**2))*(erfc((tau*(x - mu) +\
                  1.34**2)/(np.sqrt(2.)*1.34**2*tau)))/tau
        
    return(val)

def conv_line_shape_cont_noslp(x, inc, amp, mu, sig, tau, ln=True):
    """
    Calculate differential spectrum of a degraded (asymmetric) Gaussian line
    shape
    """
    x = np.array(x)
    
    if tau < 1e-3:                                  # if tau too small
        val = inc + amp*np.exp(-(x - mu)**2./(2.*sig**2)) # standard
        
    else:
        
        if ln:                                      # logarithmic approximation
            q = np.log(np.sqrt(np.pi/2.)*amp*sig/tau)
            w = (2.*tau*(x - mu) + sig**2)/(2.*tau**2)  # logarithmic term
            e = ln_erfc((tau*(x - mu) + sig**2)/(np.sqrt(2.)*sig*tau))
            edx = np.where(np.isfinite(e) == 0)     # set inf values to first
                                                    # non-inf value
            e[edx] = e[np.where(np.isfinite(e) == 1)[0][0]]

            val = inc + np.exp(q + w + e)
            
    
        else:                                       # standard degraded Gauss
            val = inc +\
                  np.sqrt(np.pi/2.)*amp*sig*np.exp((2.*tau*(x - mu) +\
                  sig**2)/(2.*tau**2))*(erfc((tau*(x - mu) +\
                  sig**2)/(np.sqrt(2.)*sig**2*tau)))/tau
        
    return(val)

def conv_line_shape_cont_noslp_nosig(x, inc, amp, mu, tau, ln=True):
    """
    Calculate differential spectrum of a degraded (asymmetric) Gaussian line
    shape
    """
    x = np.array(x)
    
    if tau < 1e-3:                                  # if tau too small
        val = inc + amp*np.exp(-(x - mu)**2./(2.*1.34**2)) # standard
        
    else:
        
        if ln:                                      # logarithmic approximation
            q = np.log(np.sqrt(np.pi/2.)*amp*1.34/tau)
            w = (2.*tau*(x - mu) + 1.34**2)/(2.*tau**2)  # logarithmic term
            e = ln_erfc((tau*(x - mu) + 1.34**2)/(np.sqrt(2.)*1.34*tau))
            edx = np.where(np.isfinite(e) == 0)     # set inf values to first
                                                    # non-inf value
            e[edx] = e[np.where(np.isfinite(e) == 1)[0][0]]

            val = inc + np.exp(q + w + e)
            
    
        else:                                       # standard degraded Gauss
            val = inc +\
                  np.sqrt(np.pi/2.)*amp*1.34*np.exp((2.*tau*(x - mu) +\
                  1.34**2)/(2.*tau**2))*(erfc((tau*(x - mu) +\
                  1.34**2)/(np.sqrt(2.)*1.34**2*tau)))/tau
        
    return(val)

def deriv_conv_line(x, amp, mu, sig, tau):
    """
    Returns the derivative of an asymmetric Gaussian line shape 
    (CONV_LINE_SHAPE) at an energy X exactly
    """
    p = np.array([amp, mu, sig, tau])
    v = - 1. / p[3]**2 * (np.exp(((2.*x - 2.*p[1])*p[3] + p[2]**2)/(2.*p[3]**2))\
          * ( - 0.5 * np.sqrt(2. * np.pi) * p[2] * np.exp(ln_erfc((p[3]*(x-p[1]) +\
          p[2]**2)/(np.sqrt(2.)*p[2]*p[3]))) + np.exp(- (p[3]*(x-p[1]) + p[2]**2)**2\
          /(2. * p[2]**2 * p[3]**2)) * p[3]) * p[0])
    return(v)

def fwhm(sig, sig_err):
    """
    Calculate the full width at half maximum of a Gaussian
    """
    v = 2.*np.sqrt(2.*np.log(2.))*sig
    verr = 2.*np.sqrt(2.*np.log(2.))*sig_err
    return(v, verr)
    
def deg_fwhm(sig, tau, sig_err, tau_err):
    """
    Approximate the full width at half maximum of a degraded Gaussian
    (Frim Kretschmer 2011 or Siegert PhD 2016)
    """
    v = 2*np.sqrt(2.*np.log(2.))
    gamma = v*sig
    a0 = 0.913735
    a1 = 0.710648
    
    fwhm = gamma*(a0 + np.sqrt((1 - a0)**2 + (a1*tau/gamma)**2))
    
    dtau = ((1 - a0)**2 + (a1*tau/gamma)**2)**-0.5 * a1**2*tau/gamma
    dsig = v*a0 + v*((1 - a0)**2 + (a1*tau/gamma)**2)**0.5 \
           + ((1 - a0)**2 + (a1*tau/gamma)**2)**-0.5 * -(a1**2*tau**2/sig**2/v)
    
    fwhm_err = np.sqrt((dsig*sig_err)**2 + (dtau*tau_err)**2)
    
    return(fwhm, fwhm_err)

def fwhm_line(sig, tau):
    """
    Approximate the full width at half maximum of a degraded Gaussian
    (From Kretschmer 2011 or Siegert PhD 2017) without uncertainty estimate
    """
    v = 2*np.sqrt(2.*np.log(2.))
    gamma = v*sig
    a0 = 0.913735
    a1 = 0.710648
    
    fwhm = gamma*(a0 + np.sqrt((1 - a0)**2 + (a1*tau/gamma)**2))
    
    return(fwhm)

def flux_line(amp,sig):
    """
    Integrated flux of a Gaussian line (degraded and not degraded)
    """
    flux = np.sqrt(2*np.pi)*amp*sig
    return(flux)


def conv_peak_val(mu, sig, tau, exact=True):
    """
    Calculate Peak position of asymmetric Gaussians as defined in
    conv_line_shape
    """
            
    if exact:
        
        if tau >= 0.05:
            xx = np.arange(1000)/1000. * tau - tau + mu
            temp = np.exp(conv_line_shape(xx, 1., mu, sig, tau))
            p0 = xx[np.where(temp == np.max(temp))][0]

            val = newton(deriv_conv_line, p0, args=(1., mu, sig, tau),
                         maxiter=int(1e5))

        else:
            val = mu - tau
        
        return(val)


def cls_plaw_function(energy,p):
    """
    Convolved (degraded) gaussian line shape of a single(!) line on top of a power-law shaped continuum.
    
    Parameters:
    :param energy: 1D-array of energies where function is evaluated (in keV)
    :param p[0]:   Normalisation (in ph/cm2/s/keV)
    :param p[1]:   Power-law index (unitless)
    :param p[2]:   Integrated line flux (in ph/cm2/s), NOT amplitude
    :param p[3]:   Peak position of symmetric Gaussian (in keV), NOT of degraded Gaussian
                   (use conv_peak_val(mu,sig,tau) to get peak of degraded Gaussian) 
    :param p[4]:   Instrumental resolution (FWHM, in keV), NOT 1-sigma value
    :param p[5]:   Instrumental degradation (tau, in keV)
    :param p[6]:   Doppler-broadening of signal (FWHM, keV), NOT 1-sigma value

    Misc:
    Ec = 1000 keV  Pivotal energy, i.e. where the power-law is normalised at (in keV)

    """
    plaw = powerlaw(energy,[p[0],p[1]])
    sigma_tot = np.sqrt(p[4]**2 + p[6]**2)/(2*np.sqrt(2*np.log(2)))
    line = p[2]/(np.sqrt(2*np.pi)*sigma_tot)*conv_line_shape(energy, 1., p[3], sigma_tot, p[5], ln=True)
    model = plaw+line
    return(model)
    