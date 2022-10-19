from COSIpy_tools import *

# Fitting model:
# define response calculation threshold globally
cut = 60

# Calculate prior
def lnprior_PSsearchFlux(theta,zs):
    """
    Return log-prior probability of a parameter set theta,
    setting normal, weakly informative, priors on the flux
    and background scaling, plus a uniform circular region
    inside the field of view, and normal priors on the
    position (also weakly informative)
    :param: theta     Array of fitted parameters with
                      theta[0] = F_sky (sky flux)
                      theta[1] = A_bg (BG scaling)
                      theta[2] = Gal. Lon. of source (deg)
                      theta[3] = Gal. Lat. of source (deg)
    :param: zs        Array of pointing directions to calculate
                      threshold for outside field of view
    
    """
    # define fit parameters
    F_sky,A_bg,l,b = theta
    # normal prior on flux (typically 1e-6..1e-1)
    mu_F_sky = 0
    sigma_F_sky = 10
    # normal prior on BG
    mu_A_bg = 0
    sigma_A_bg = 10000

    # calculate minimum and maximum observation (and mean)
    lzm = minmax(zs[:,0])
    bzm = minmax(zs[:,1])
    mean_zl = np.mean(zs[:,0])
    mean_zb = np.mean(zs[:,1])
    
    # normal prior on source position (mean of observation) as centre
    mu_l = mean_zl
    sigma_l = 20 # within 3 sigma, (= 60 deg), everything can be captured
    mu_b = mean_zb
    sigma_b = 20
    
    # get distance to outermost observations
    dist = angular_distance(lzm,bzm,l,b,deg=True)
    
    """if l > 180:
        l -= 360
    if l <= -180:
        l += 360"""
    
    # if distance of source is larger than 60 deg to outermost points, return -infinity
    # also check for inside "one sky"
    if ((dist[1] < cut) & (dist[0] < cut)) & (-90 <= b <= 90) & (-180 < l <= 180):
        return -0.5*(F_sky-mu_F_sky)**2/sigma_F_sky**2 -0.5*(A_bg-mu_A_bg)**2/sigma_A_bg**2 -0.5*(l-mu_l)**2/sigma_l**2 -0.5*(b-mu_b)**2/sigma_b**2
    return -np.inf
    
# Calculate likelihood
def lnlike_PSsearchFlux(theta, data, response_sky, xs, ys, zs, response_bg, Texp, area_flag):
    """
    Return the Poisson log-likelihood of the fitting model
    with parameters theta, given the data. Using instrument
    orientation, to calculate sky response on the fly and
    background model response to fit for flux, BG, l, b.
    :param: theta           Array of fitted parameters with
                            theta[0] = F_sky (sky flux)
                            theta[1] = A_bg (BG scaling)
                            theta[2] = Gal. Lon. of source (deg)
                            theta[3] = Gal. Lat. of source (deg)
    :param: response_sky    Response grid with regular sky dimension
    :param: xs              Array of pointing x-directions
    :param: ys              Array of pointing y-direction
    :param: zs              Array of pointing directions to calculate
                            threshold for outside field of view
    :param: response_bg     Background response, valid for all times
    :param: Texp            Exposure/observation time in seconds
    :param: area_flag       Option to calculate response by area
                            (which is buggy) instead of by distance
    """
    # define fit parameters
    F_sky,A_bg,l,b = theta
    
    # define zeniths and azimuths for given source position l/b
    zens,azis = zenazi(xs[:,0],xs[:,1],
                       ys[:,0],ys[:,1],
                       zs[:,0],zs[:,1],
                       l,b)
    
    # check which response calculation should be used
    if area_flag == True:
        sky_response = np.mean(get_response_with_weights_area(response_sky,zens,azis,cut=cut),axis=0)
    else:
        sky_response = np.mean(get_response_with_weights(response_sky,zens,azis,cut=cut),axis=0)
    
    # calculate sky model count expectaion
    model_sky = F_sky*Texp*sky_response
    # calculate BG model count expectation
    model_bg = A_bg*response_bg
    # add together
    model_tot = model_sky + model_bg

    # check for any negative model counts to return -infinity (zero counts model predictions are allowed)
    if np.any(model_tot < 0):
        return -np.inf
    else:
        # calculate Poisson likelihood
        stat = -2*np.sum(model_tot - data*np.nan_to_num(np.log(model_tot)))
        return stat

# Calculate posterior
def lnprob_PSsearchFlux(theta, data, response_sky, xs, ys, zs, response_bg,Texp,area_flag):
    """
    Return the log-posterior of the fitting model using above-
    defined likelihood and prior.
    :param: theta           Array of fitted parameters with
                            theta[0] = F_sky (sky flux)
                            theta[1] = A_bg (BG scaling)
                            theta[2] = Gal. Lon. of source (deg)
                            theta[3] = Gal. Lat. of source (deg)
    :param: response_sky    Response grid with regular sky dimension
    :param: xs              Array of pointing x-directions
    :param: ys              Array of pointing y-direction
    :param: zs              Array of pointing directions to calculate
                            threshold for outside field of view
    :param: response_bg     Background response, valid for all times
    :param: Texp            Exposure/observation time in seconds
    :param: area_flag       Option to calculate response by area
                            (which is buggy) instead of by distance
    """
    # get prior
    lp = lnprior_PSsearchFlux(theta,zs)
    if not np.isfinite(lp):
        return -np.inf
    # add prior to likelihood and return posterior
    return lp + lnlike_PSsearchFlux(theta, data, response_sky, xs, ys, zs, response_bg,Texp,area_flag)