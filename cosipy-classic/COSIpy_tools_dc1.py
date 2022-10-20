import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import astropy.units as u
import sys
import time
import scipy.interpolate as interpol
import matplotlib
from matplotlib import ticker, cm
from scipy.ndimage import gaussian_filter
from tqdm import tqdm_notebook as tqdm
import os, glob 

import pandas as pd

import ROOT as M
# Load MEGAlib into ROOT
M.gSystem.Load("$(MEGAlib)/lib/libMEGAlib.so")
# Initialize MEGAlib
G = M.MGlobal()
G.Initialize()


    
def one_func(x,y,grid=False):
    """Returns an array of ones for any given input array (or scalar) x.
    :param: x       Input array (or scalar, tuple, list)
    :param: y       Optional second array or value
    :option: grid   Standard keyword to work with RectBivarSpline
    """
    if isinstance(x, (list, tuple, np.ndarray)):
        return np.ones(len(x))
    else:
        return 1.
    
    
    
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



def GreatCircleGrid(l1,b1,l2,b2,deg=True):
    """
    Calculate the Great Circle length on a sphere from longitude/latitude pairs to others
    in units of rad on a unit sphere
    :param: l1    longitude of points Ai
    :param: b1    latitude of points Ai
    :param: l2    longitude of point Bj
    :param: b2    latitude of point Bj
    :option: deg  Default True to convert degree input to radians for trigonometric function use
                  If False, radian input is assumed
    """
    if deg == True:
        l1,b1,l2,b2 = np.deg2rad(l1),np.deg2rad(b1),np.deg2rad(l2),np.deg2rad(b2)

    L1, L2 = np.meshgrid(l1,l2)
    B1, B2 = np.meshgrid(b1,b2)
        
    return np.sin(B1)*np.sin(B2) + np.cos(B1)*np.cos(B2)*np.cos(L1-L2)



def zenaziGrid(scx_l, scx_b, scy_l, scy_b, scz_l, scz_b, src_l, src_b):
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
    :param: src_l      SOURCEgrid longitudes
    :param: src_b      SOURCEgrid latitudes
    
    """
#     # make matrices for response calculation on a pre-defined grid
    SCZ_L, SRC_L = np.meshgrid(scz_l,src_l)
    SCZ_B, SRC_B = np.meshgrid(scz_b,src_b)
    # Zenith is the distance from the optical axis (here z)
    costheta = GreatCircleGrid(scz_l,scz_b,src_l,src_b)                                                                        
    # Azimuth is the combination of the remaining two
    
    SCX_L, SRC_L = np.meshgrid(scx_l,src_l)
    SCX_B, SRC_B = np.meshgrid(scx_b,src_b)    
    cosx = GreatCircle(SCX_L,SCX_B,SRC_L,SRC_B)
    SCY_L, SRC_L = np.meshgrid(scy_l,src_l)
    SCY_B, SRC_B = np.meshgrid(scy_b,src_b)  
    cosy = GreatCircle(SCY_L,SCY_B,SRC_L,SRC_B)
    
    # theta = zenith
    theta = np.rad2deg(np.arccos(costheta))
    # phi = azimuth
    phi = np.rad2deg(np.arctan2(cosx,cosy))
    
    # make azimuth going from 0 to 360 deg
    if phi.size == 1:
        if (phi < 0):
            phi += 360
    else:
        phi[phi < 0] += 360
    
    return theta,phi    
                                                                      


def cashstat(data,model):
    return -2*np.sum(data*np.log(model)-model-(data*np.nan_to_num(np.log(data))-data))


