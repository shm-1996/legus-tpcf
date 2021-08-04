"""
Tools for computing two-point correlation functions. 
Adopted from astromL two-point correlation module. 

Modified to account for edge/chip effects in field of view by 
masking out regions in the random catalog. The masks are prepared using
StScI's footprint finder. 

AUTHOR (modified)
Shyam Harimohan Menon (2020)
"""

import warnings
import numpy as np
from sklearn.neighbors import BallTree
from astroML.utils import check_random_state
from regions import read_ds9
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import regions
from header import *
# Check if scikit-learn's two-point functionality is available.
# This was added in scikit-learn version 0.14
try:
    from sklearn.neighbors import KDTree
    sklearn_has_two_point = True
except ImportError:
    import warnings
    sklearn_has_two_point = False


def exponential_sample(galaxy_class,r_c=50.,age=None,len_random=10,age_cut=1.e7):
    from GalaxyModels import Create_Galaxy
    from PaperPlots import get_separations
    from Plot_Class import myPlot
    name = galaxy_class.name.split('_')[0] + ' ' + galaxy_class.name.split('_')[1]

    #Center coordinates
    ra_dec = SkyCoord.from_name(name)
    ra_centre,dec_centre = ra_dec.ra.value,ra_dec.dec.value
    #Deproject FOV to prepare mask
    hdu = fits.open(galaxy_class.fits_file)[0]
    wcs = WCS(hdu.header)
    xpix_c,ypix_c = wcs.all_world2pix(ra_dec.ra.value,ra_dec.dec.value,0)

    galaxy_class.get_ra_dec(age=age,age_cut=age_cut)
    galaxy_ra,galaxy_dec = galaxy_class.ra,galaxy_class.dec
    
    xyz_centre = np.array(ra_dec_to_xyz(ra_centre,dec_centre))
    xyz_data = np.array(ra_dec_to_xyz(galaxy_ra,galaxy_dec))
    
    # Radial distances from centre in deg
    distance_data = np.sqrt((xyz_data[0]-xyz_centre[0])**2 + (xyz_data[1]-\
            xyz_centre[1])**2 + (xyz_data[2]-xyz_centre[2])**2)
    
    #Convert arcsec to degree and set R_max = L_FOV/2.
    R_max = np.max(distance_data)*u.radian.to(u.deg)
    r_min = np.min(distance_data)*u.radian.to(u.deg)
    #Parameters of exponential disk
    r_c = r_c*u.arcsec.to(u.deg)
    z_h = 2.0*r_c        

    factor = 2.0
    size_r = len_random
    while(size_r<=len_random):
        size = int(factor*len_random)

        ra_r_raw,dec_r_raw = Create_Galaxy(r_c=r_c,z_h=z_h,i=galaxy_class.inclination,R_max=R_max,r_min=r_min,size=size)
        dec_r_raw = dec_r_raw/np.cos(np.deg2rad(galaxy_class.inclination))
        ra_r_raw = ra_r_raw + ra_centre
        dec_r_raw = dec_r_raw + dec_centre
        
        #Deproject this observed galaxy

        random_x,random_y = wcs.wcs_world2pix(ra_r_raw,dec_r_raw,0)
        
        random_xy = np.vstack((random_x,random_y)).T

        #Mask with FOV   
        region = read_ds9(galaxy_class.region_file)
        mask_xy = None
        for i in range(0,np.size(region)):
            region_dep = deproject_region(region,i,galaxy_class.pa,galaxy_class.inclination,xpix_c,ypix_c)
            if(region_dep is not None):
                if(mask_xy is None):
                    mask_xy = region_dep.contains(regions.PixCoord(random_x,random_y))
                else:
                    mask_xy += region_dep.contains(regions.PixCoord(random_x,random_y))

        random_xy_masked = random_xy[mask_xy]
        ra_masked,dec_masked = random_xy_masked[:,0], random_xy_masked[:,1]
        ra_r,dec_r = wcs.all_pix2world(ra_masked,dec_masked,0)
        size_r = np.size(ra_r)
        if(size_r<len_random):
            factor *=2
    
    #Randomly cull to set len_random points in random sample
    indices = np.random.randint(0,size_r,len_random)
    ra_r = ra_r[indices]
    dec_r = dec_r[indices]
    
    return ra_r,dec_r


def random_distribution_probability(RA,DEC,distance_hist,xyz_centre,galaxy_class) :
    """
    Returns the probability that a random RA,DEC combination should exist, 
    by taking into account the distribution of distances for the clusters in the
    data from the galactic centre.  
    Parameters
    ----------
    RA  : ndarrays
        RA values of the point sources
        units are degrees
    DEC : ndarrays
        DEC values of the point sources
    distance_hist : ndarray
        The histogram ndarray returned by np.histogram for the distribution of
        distances from the galactic centre for the data, to which the random
        sample should be matched.
    xyz_centre : tuple
        x,y,z coordinates of the galactic centre 
    region_file: string
        Region file of HST camera footprint
    fits_file : string
        Fits file for WCS information
    Returns
    -------
    rand_prob : ndarray of size len(RA)
        the probability of the combination of RA/DEC
    """

    rand_prob = np.zeros(len(RA))
    prob_bins = distance_hist[0]
    bins_radial = distance_hist[1]
    
    #Check whether point lies in masked region
    region = read_ds9(galaxy_class.region_file)
    # For converting pixel region to sky region read the fits file
    hdu = fits.open(galaxy_class.fits_file)[0]
    wcs = WCS(hdu.header)

    random_x,random_y = wcs.wcs_world2pix(RA,DEC,0)

    xyz_array = np.array(ra_dec_to_xyz(RA,DEC))
    distance = np.sqrt((xyz_array[0]-xyz_centre[0])**2 + (xyz_array[1]-\
            xyz_centre[1])**2 + (xyz_array[2]-xyz_centre[2])**2)
    
    #separation in kpc
    distance = distance*galaxy_class.distance*1.e6*const.Parsec/(1.e3*const.Parsec)
    
    #Put into bins
    distance_bins = np.digitize(distance,bins_radial)
    distance_bins = distance_bins-1
    # Set prob of points above max distance = 0.0
    prob_bins = np.append(prob_bins,0.0)
    
    #Get probability of sample associated with distance
    rand_prob = prob_bins[distance_bins]
    
    #Compute FOV mask
    name = galaxy_class.name.split('_')[0] + ' ' + galaxy_class.name.split('_')[1]
    ra_dec = SkyCoord.from_name(name)
    xpix_c,ypix_c = wcs.all_world2pix(ra_dec.ra.value,ra_dec.dec.value,0)
    mask_xy = None
    for i in range(0,np.size(region)):
        region_dep = deproject_region(region,i,galaxy_class.pa,galaxy_class.inclination,xpix_c,ypix_c)
        if(region_dep is not None):
            if(mask_xy is None):
                mask_xy = region_dep.contains(regions.PixCoord(random_x,random_y))
            else:
                mask_xy += region_dep.contains(regions.PixCoord(random_x,random_y))


    #Set points in mask with 0.0 probability
    rand_prob[np.where(mask_xy == False)] = 0.0
    
    return rand_prob

def masked_radial_random_sample(galaxy_class,len_random=100) :
    """
    Prepare a random catalog that takes into account the FOV edge effects
    and ACS chip gaps, as well as match the distribution of the distances of 
    the clusters in the random sample to that of the data. 
    Parameters
    ----------
    galaxy_class : class Galaxy
        Instance of class Galaxy that contains metadata of the galaxy_class.
    len_random : integer
        Number of attempts of points in the random array. 
        Keep at least > 100*len(data)
    Returns
    -------
    RA, DEC : ndarray
        the random sample with len_random items. 
    """

    #Compute distance to galactic centre
    name = galaxy_class.name.split('_')[0] + ' ' + galaxy_class.name.split('_')[1]
    ra_dec = SkyCoord.from_name(name)
    ra_centre,dec_centre = ra_dec.ra.value,ra_dec.dec.value

    #RA/DEC
    RA = galaxy_class.ra
    DEC = galaxy_class.dec
    
    #xyz coordinates 
    xyz_data = np.array(ra_dec_to_xyz(RA,DEC))
    xyz_centre = np.array(ra_dec_to_xyz(ra_centre,dec_centre))

    # Relative separation from centre of galaxy
    distance_data = np.sqrt((xyz_data[0]-xyz_centre[0])**2 + (xyz_data[1]-\
            xyz_centre[1])**2 + (xyz_data[2]-xyz_centre[2])**2)

    #separation in kpc
    #Galaxy distance is in Mpc: convert to kpc
    distance_data = distance_data*galaxy_class.distance*1.e6*const.Parsec/(1.e3*const.Parsec)
    
    # 50 Bins b/w 0.01 to 11 kpc
    bins_radial = np.linspace(0.01,11,50)
    hist,bins = np.histogram(distance_data,bins=bins_radial)
    distance_hist = (hist/(np.size(distance_data)),bins)
    #Random RA/DEC
    ra_rand, dec_rand = uniform_sphere(galaxy_class,len_random,ignore_deproject=True)

    #Compute probabilities of existing
    rand_prob = random_distribution_probability(ra_rand,dec_rand,distance_hist,
        xyz_centre,galaxy_class)

    #Random nos b/w 0 & 1
    rand_values = np.random.random(len(ra_rand))
    
    #Metropolis-Hastings Algorithm condition check
    # i.e. if(random_no<prob_of_point) -> accept else discard
    accepted_indices = np.where(rand_prob>rand_values)

    # Pick out accepted combinations of RA/DEC
    ra_accepted = ra_rand[accepted_indices]
    dec_accepted = dec_rand[accepted_indices]

    return ra_accepted,dec_accepted

def masked_random_sample(galaxy_class,len_random=100):
    """
    Prepare a random catalog that takes into account the FOV edge effects
    and ACS chip gaps. 

    Parameters
    ----------
    galaxy : class Galaxy
        Instance of class Galaxy that contains metadata of the galaxy.
    len_random : integer
        Number of attempts of points in the random array. 
        Keep at least > 100*len(data)

    Returns
    -------
    RA, DEC : ndarray
        the random sample with len_random items. 
    """

    # DS9 Region file to select footprint region
    region = read_ds9(galaxy_class.region_file)
    fits_file = galaxy_class.fits_file
    RA = galaxy_class.ra
    DEC = galaxy_class.dec

    # For converting pixel region to sky region read the fits file
    hdu = fits.open(fits_file)[0]
    wcs = WCS(hdu.header)

    #Compute corresponding pixel limits 
    min_x,min_y = wcs.wcs_world2pix(min(RA),min(DEC),0)
    max_x,max_y = wcs.wcs_world2pix(max(RA),max(DEC),0)
    min_x,max_x,min_y,max_y = np.int(min_x.item()),np.int(max_x.item()),\
    np.int(min_y.item()),np.int(max_y.item())

    if(max_x<min_x) :
        min_x,max_x = max_x,min_x
    if(max_y<min_y) :
        min_y,max_y = max_y,min_y
        
    random_x = min_x+ np.random.random(len_random)*(max_x-min_x)
    random_y = min_y+ np.random.random(len_random)*(max_y-min_y)
    random_xy = np.vstack((random_x,random_y)).T

    name = galaxy_class.name.split('_')[0] + ' ' + galaxy_class.name.split('_')[1]
    ra_dec = SkyCoord.from_name(name)
    xpix_c,ypix_c = wcs.all_world2pix(ra_dec.ra.value,ra_dec.dec.value,0)
    mask_xy = None
    for i in range(0,np.size(region)):
        region_dep = deproject_region(region,i,galaxy_class.pa,galaxy_class.inclination,xpix_c,ypix_c)
        if(region_dep is not None):
            if(mask_xy is None):
                mask_xy = region_dep.contains(regions.PixCoord(random_x,random_y))
            else:
                mask_xy += region_dep.contains(regions.PixCoord(random_x,random_y))

    random_xy_masked = random_xy[mask_xy]
    ra_masked,dec_masked = random_xy_masked[:,0], random_xy_masked[:,1]
    ra_random,dec_random = wcs.all_pix2world(ra_masked,dec_masked,0)
                
    return ra_random, dec_random


def deproject_region(region,i,pa,inclination,xpix_c,ypix_c):
    #rotate clockwise by angle PA about centre of galaxy
    region_rotated = region[i].rotate(regions.PixCoord(xpix_c,ypix_c),-pa*u.deg)
    try:
        size = np.size(region_rotated.vertices)
    except:
        return None
    x = np.zeros(size)
    y = np.zeros(size)
    for i in range(0,size):
        x[i] = region_rotated.vertices[i].x/np.cos(np.deg2rad(inclination)) -\
        xpix_c/np.cos(np.deg2rad(inclination)) + xpix_c

        y[i] = region_rotated.vertices[i].y
    regions_dep = regions.PolygonPixelRegion(vertices=regions.PixCoord(x=x,y=y))    
    return regions_dep

def uniform_sphere(galaxy, len_random=100,ignore_deproject=False):
    """Draw a uniform sample on a sphere

    Parameters
    ----------
    galaxy : class Galaxy
        Instance of class Galaxy that contains metadata of the galaxy.
    len_random : integer
        Number of attempts of points in the random array. 
        Keep at least > 100*len(data)
    ignore_deproject : bool
        Flag to ignore the deprojection of the random sample to account
        for inclination. This is used by the masked_radial sample to just obtain 
        a rectangular region in RA-DEC space. 

    Returns
    -------
    RA, DEC : ndarray
        the random sample on the sphere within the given limits.
        arrays have shape equal to size.
    """

    # RA/DEC Limits
    RAlim = min(galaxy.ra),max(galaxy.ra)
    DEClim = min(galaxy.dec),max(galaxy.dec)

    zlim = np.sin(np.pi * np.asarray(DEClim) / 180.)

    z = zlim[0] + (zlim[1] - zlim[0]) * np.random.random(len_random)
    DEC = (180. / np.pi) * np.arcsin(z)
    RA = RAlim[0] + (RAlim[1] - RAlim[0]) * np.random.random(len_random)

    if(galaxy.deproject_galaxy == True and ignore_deproject == False):
        ra_random = RA*np.cos(galaxy.pa) - DEC*np.sin(galaxy.pa)
        dec_random = RA*np.sin(galaxy.pa) + DEC*np.cos(galaxy.pa)
        ra_random = ra_random/(np.cos(galaxy.inclination))
    else : 
        ra_random, dec_random = RA, DEC


    return ra_random,dec_random


def ra_dec_to_xyz(ra, dec):
    """Convert ra & dec to Euclidean points

    Parameters
    ----------
    ra, dec : ndarrays

    Returns
    x, y, z : ndarrays
    """
    sin_ra = np.sin(ra * np.pi / 180.)
    cos_ra = np.cos(ra * np.pi / 180.)

    sin_dec = np.sin(np.pi / 2 - dec * np.pi / 180.)
    cos_dec = np.cos(np.pi / 2 - dec * np.pi / 180.)

    return (cos_ra * sin_dec,
            sin_ra * sin_dec,
            cos_dec)


def angular_dist_to_euclidean_dist(D, r=1):
    """convert angular distances to euclidean distances"""
    return 2 * r * np.sin(0.5 * D * np.pi / 180.)


def two_point(data, bins, method='standard',
              data_R=None, random_state=None):
    """Two-point correlation function

    Parameters
    ----------
    data : array_like
        input data, shape = [n_samples, n_features]
    bins : array_like
        bins within which to compute the 2-point correlation.
        shape = Nbins + 1
    method : string
        "standard" or "landy-szalay".
    data_R : array_like (optional)
        if specified, use this as the random comparison sample
    random_state : integer, np.random.RandomState, or None
        specify the random state to use for generating background

    Returns
    -------
    corr : ndarray
        the estimate of the correlation function within each bin
        shape = Nbins
    """
    data = np.asarray(data)
    bins = np.asarray(bins)
    rng = check_random_state(random_state)

    if method not in ['standard', 'landy-szalay']:
        raise ValueError("method must be 'standard' or 'landy-szalay'")

    if bins.ndim != 1:
        raise ValueError("bins must be a 1D array")

    if data.ndim == 1:
        data = data[:, np.newaxis]
    elif data.ndim != 2:
        raise ValueError("data should be 1D or 2D")

    n_samples, n_features = data.shape
    Nbins = len(bins) - 1

    # shuffle all but one axis to get background distribution
    if data_R is None:
        data_R = data.copy()
        for i in range(n_features - 1):
            rng.shuffle(data_R[:, i])
    else:
        data_R = np.asarray(data_R)
        if (data_R.ndim != 2) or (data_R.shape[-1] != n_features):
            raise ValueError('data_R must have same n_features as data')

    factor = len(data_R) * 1. / len(data)

    if sklearn_has_two_point:
        # Fast two-point correlation functions added in scikit-learn v. 0.14
        KDT_D = KDTree(data)
        KDT_R = KDTree(data_R)

        counts_DD = KDT_D.two_point_correlation(data, bins)
        counts_RR = KDT_R.two_point_correlation(data_R, bins)

    else:
        warnings.warn("Version 0.3 of astroML will require scikit-learn "
                      "version 0.14 or higher for correlation function "
                      "calculations. Upgrade to sklearn 0.14+ now for much "
                      "faster correlation function calculations.")

        BT_D = BallTree(data)
        BT_R = BallTree(data_R)

        counts_DD = np.zeros(Nbins + 1)
        counts_RR = np.zeros(Nbins + 1)

        for i in range(Nbins + 1):
            counts_DD[i] = np.sum(BT_D.query_radius(data, bins[i],
                                                    count_only=True))
            counts_RR[i] = np.sum(BT_R.query_radius(data_R, bins[i],
                                                    count_only=True))

    DD = np.diff(counts_DD)
    RR = np.diff(counts_RR)

    # check for zero in the denominator
    RR_zero = (RR == 0)
    RR[RR_zero] = 1

    if method == 'standard':
        corr = factor ** 2 * DD / RR - 1
    elif method == 'landy-szalay':
        if sklearn_has_two_point:
            counts_DR = KDT_R.two_point_correlation(data, bins)
        else:
            counts_DR = np.zeros(Nbins + 1)
            for i in range(Nbins + 1):
                counts_DR[i] = np.sum(BT_R.query_radius(data, bins[i],
                                                        count_only=True))
        DR = np.diff(counts_DR)

        corr = (factor ** 2 * DD - 2 * factor * DR + RR) / RR

    corr[RR_zero] = np.nan

    return corr



def bootstrap_two_point(data, bins, Nbootstrap=10,
                        method='standard', return_bootstraps=False,
                        random_state=None):
    """Bootstrapped two-point correlation function

    Parameters
    ----------
    data : array_like
        input data, shape = [n_samples, n_features]
    bins : array_like
        bins within which to compute the 2-point correlation.
        shape = Nbins + 1
    Nbootstrap : integer
        number of bootstrap resamples to perform (default = 10)
    method : string
        "standard" or "landy-szalay".
    return_bootstraps: bool
        if True, return full bootstrapped samples
    random_state : integer, np.random.RandomState, or None
        specify the random state to use for generating background

    Returns
    -------
    corr, corr_err : ndarrays
        the estimate of the correlation function and the bootstrap
        error within each bin. shape = Nbins
    """
    data = np.asarray(data)
    bins = np.asarray(bins)
    rng = check_random_state(random_state)

    if method not in ['standard', 'landy-szalay']:
        raise ValueError("method must be 'standard' or 'landy-szalay'")

    if bins.ndim != 1:
        raise ValueError("bins must be a 1D array")

    if data.ndim == 1:
        data = data[:, np.newaxis]
    elif data.ndim != 2:
        raise ValueError("data should be 1D or 2D")

    if Nbootstrap < 2:
        raise ValueError("Nbootstrap must be greater than 1")

    n_samples, n_features = data.shape

    # get the baseline estimate
    corr = two_point(data, bins, method=method, random_state=rng)

    bootstraps = np.zeros((Nbootstrap, len(corr)))

    for i in range(Nbootstrap):
        indices = rng.randint(0, n_samples, n_samples)
        bootstraps[i] = two_point(data[indices, :], bins, method=method,
                                  random_state=rng)

    # use masked std dev in case of NaNs
    corr_err = np.asarray(np.ma.masked_invalid(bootstraps).std(0, ddof=1))

    if return_bootstraps:
        return corr, corr_err, bootstraps
    else:
        return corr, corr_err



def two_point_angular(ra, dec, bins, method='standard', random_state=None):
    """Angular two-point correlation function

    A separate function is needed because angular distances are not
    euclidean, and random sampling needs to take into account the
    spherical volume element.

    Parameters
    ----------
    ra : array_like
        input right ascention, shape = (n_samples,)
    dec : array_like
        input declination
    bins : array_like
        bins within which to compute the 2-point correlation.
        shape = Nbins + 1
    method : string
        "standard" or "landy-szalay".
    random_state : integer, np.random.RandomState, or None
        specify the random state to use for generating background

    Returns
    -------
    corr : ndarray
        the estimate of the correlation function within each bin
        shape = Nbins
    """
    ra = np.asarray(ra)
    dec = np.asarray(dec)
    rng = check_random_state(random_state)

    if method not in ['standard', 'landy-szalay']:
        raise ValueError("method must be 'standard' or 'landy-szalay'")

    if bins.ndim != 1:
        raise ValueError("bins must be a 1D array")

    if (ra.ndim != 1) or (dec.ndim != 1) or (ra.shape != dec.shape):
        raise ValueError('ra and dec must be 1-dimensional '
                         'arrays of the same length')

    n_features = len(ra)
    Nbins = len(bins) - 1

    # draw a random sample with N points
    ra_R, dec_R = uniform_sphere((min(ra), max(ra)),
                                 (min(dec), max(dec)),
                                 2 * len(ra))

    data = np.asarray(ra_dec_to_xyz(ra, dec), order='F').T
    data_R = np.asarray(ra_dec_to_xyz(ra_R, dec_R), order='F').T

    # convert spherical bins to cartesian bins
    bins_transform = angular_dist_to_euclidean_dist(bins)

    return two_point(data, bins_transform, method=method,
                     data_R=data_R, random_state=rng)



def bootstrap_two_point_angular(galaxy, method='standard',
                                Nbootstraps=10, random_state=None,
                                random_method='masked',r_c=50.):
    """Angular two-point correlation function

    A separate function is needed because angular distances are not
    euclidean, and random sampling needs to take into account the
    spherical volume element.

    Parameters
    ----------
    galaxy : class Galaxy
        Instance of class Galaxy describing the properties of the galaxy
        for which TPCF to be calculated.
    method : string
        "standard" or "landy-szalay".
    Nbootstraps : int
        number of bootstrap resamples
    random_state : integer, np.random.RandomState, or None
        specify the random state to use for generating background
    random_method : string
        "uniform", "masked" if random sample should sample FOV, or "masked_radial"
        if random sample should also mirror data distribution in distance to 
        galactic centre. 

    Returns
    -------
    corr : ndarray
        the estimate of the correlation function within each bin
        shape = Nbins
    dcorr : ndarray
        error estimate on dcorr (sample standard deviation of
        bootstrap resamples)
    bootstraps : ndarray
        The full sample of bootstraps used to compute corr and dcorr
    """

    #Initialise
    ra = galaxy.ra
    dec = galaxy.dec
    bins = galaxy.bins
    rng = check_random_state(random_state)

    #Safety checks
    if method not in ['standard', 'landy-szalay']:
        raise ValueError("method must be 'standard' or 'landy-szalay'")
    if random_method not in ['uniform', 'masked','masked_radial','exponential']:
        raise ValueError("method must be 'uniform' or 'masked' or 'masked_radial' or 'exponential' ")

    if bins.ndim != 1:
        raise ValueError("bins must be a 1D array")

    if (ra.ndim != 1) or (dec.ndim != 1) or (ra.shape != dec.shape):
        raise ValueError('ra and dec must be 1-dimensional '
                         'arrays of the same length')

    n_features = len(ra)
    Nbins = len(bins) - 1
    data = np.asarray(ra_dec_to_xyz(ra, dec), order='F').T

    # convert spherical bins to cartesian bins
    bins_transform = angular_dist_to_euclidean_dist(bins)

    bootstraps = []

    for i in range(Nbootstraps):
        if(random_method == 'uniform') :
            #Draw uniform random sample with N points
            ra_R, dec_R = uniform_sphere(galaxy,2 * len(ra))
                                      
        elif(random_method == 'masked') :
            #Draw a masked random sample with N points
            ra_R, dec_R = masked_random_sample(galaxy,10*len(ra))

        elif(random_method == 'masked_radial') :
            #Draw a random sample 
            ra_R, dec_R = masked_radial_random_sample(galaxy,len_random=100*len(ra))
        elif(random_method == 'exponential'):
            ra_R, dec_R = exponential_sample(galaxy,r_c=r_c,len_random=len(ra))

        data_R = np.asarray(ra_dec_to_xyz(ra_R, dec_R), order='F').T

        if i > 0:
            # random sample of the data
            ind = np.random.randint(0, data.shape[0], data.shape[0])
            data_b = data[ind]
        else:
            data_b = data

        bootstraps.append(two_point(data_b, bins_transform, method=method,
                                    data_R=data_R, random_state=rng))

    bootstraps = np.asarray(bootstraps)
    corr = np.mean(bootstraps, 0)
    corr_err = np.std(bootstraps, 0, ddof=1)

    return corr, corr_err, bootstraps


def linear_function(theta,A_1,alpha_1,alpha_2,beta) :
    A_2 = A_1 
    function = np.piecewise(theta,[np.log(theta)<beta],[lambda theta :  A_1 + alpha_1 * np.log(theta), 
        lambda theta : A_2 + (alpha_1-alpha_2)*beta + alpha_2*np.log(theta)])
    return function

def smooth_function(theta,A1,alpha_1,alpha_2,beta):
    function = A1*((theta/beta)**alpha_1 + (theta/beta)**alpha_2)
    return function

def onepowerlaw_function(theta,A1,alpha_1) :
    function = A1 + alpha_1*np.log(theta)
    return function

def linear_truncation(theta,A1,alpha_1,theta_c) :
    function = A1 + alpha_1*np.log(theta) - theta/theta_c
    return function    

def piecewise_truncation(theta,A_1,alpha_1,alpha_2,beta,theta_c):
    A_2 = A_1 
    function = np.piecewise(theta,[np.log(theta)<beta],[lambda theta :  A_1 + alpha_1 * np.log(theta), 
        lambda theta : A_2 + (alpha_1-alpha_2)*beta + alpha_2*np.log(theta) - theta/theta_c])
    return function

        
