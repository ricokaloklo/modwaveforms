import numpy as np
import mpmath as mp
from scipy.special import gamma
from bilby.gw.source import *
from .utils import point_mass as pm
from .utils.constants import *

""" Point lens waveform """

"Diffraction integral"
#Amplification function
def laguerre(n,a,x):
    """ Laguerre polynomial of degree n with parameter a at point x 
    
    Parameters
    ----------
    n : int
        Degree of the Laguerre polynomial
    a : float
        Parameter of the Laguerre polynomial
    x : float
        Point at which to evaluate the Laguerre polynomial

    Returns
    -------
    float
        Value of the Laguerre polynomial at x
    """

    nplaguerre = np.frompyfunc(mp.laguerre,3,1)
    return nplaguerre(n,a,x)
vlaguerre = np.vectorize(laguerre)

def F_pointlens(f,ML,y):
    """ Amplification function of a point lens

    Follows definition of Eq.16 https://arxiv.org/abs/2005.10702

    Parameters
    ----------
    fs : float
        Frequency of the gravitational wave in Hz
    ML : float
        Redshifted mass of the lens in solar masses
    y : float
        Dimensionless impact parameter of the lens

    Returns
    -------
    complex
        Amplification function of the point lens
    """

    w = 2.*np.pi * (4*TSUN*ML) * f

    lague = np.array(vlaguerre(-0.5j*w, 0, 0.5j*w*y**2.),dtype='complex')
    ampl = np.power(-0.5j,1.+0.5j*w)*np.power(w,1.+0.5j*w)*gamma(-0.5j*w)*lague
    ampl[w==0. + 0.j] = 1. #gamma of 0 is inf but F(0,b)=1.

    t_plus_pm = pm.t_delay_geom_plus(y)* pm.t_ref(ML) #+ image time delay in seconds
    ampl *= np.exp(-1j*t_plus_pm*2*np.pi*f) #eliminate global time delay of the + image
    return ampl

def pointlens(frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
    phi_12, a_2, tilt_2, phi_jl, theta_jn, phase,ML,y, **kwargs
    ):

    """ Waveform generation for a binary black hole merger lensed by a point mass lens

    Generate the frequency domain waveform for a binary black hole merger lensed by a point mass lens

    Parameters
    ----------
    frequency_array : array
        Array of frequencies in Hz
    mass_1 : float
        Mass of the primary black hole in solar masses
    mass_2 : float
        Mass of the secondary black hole in solar masses
    luminosity_distance : float
        Luminosity distance to the source in Mpc
    a_1 : float
        Dimensionless spin magnitude of the primary black hole
    tilt_1 : float
        Tilt angle of the primary black hole spin in radians
    phi_12 : float
        Azimuthal angle between the two black hole spins in radians
    a_2 : float
        Dimensionless spin magnitude of the secondary black hole
    tilt_2 : float
        Tilt angle of the secondary black hole spin in radians
    phi_jl : float
        Azimuthal angle between the total angular momentum and the line of sight in radians
    theta_jn : float
        Angle between the total angular momentum and the line of sight in radians
    phase : float
        Reference phase of the binary in radians
    ML : float
        Redshifted mass of the lens in solar masses
    y : float
        Dimensionless impact parameter of the lens
    **kwargs : dict
        Additional keyword arguments to pass to the waveform generator

    Returns
    -------
    wf : dict
        Dictionary containing the plus and cross polarizations of the lensed waveform
    """

    frequency_domain_source_model = lal_binary_black_hole

    # Actually generate the waveform by calling the generator
    wf = frequency_domain_source_model(
        frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
        phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, **kwargs
    )

    # Calculate the amplification factor
    # complex conjugate because of different sign convention in the Fourier transform
    F = F_pointlens(frequency_array,ML,y)

    wf["plus"] = wf["plus"]*F
    wf["cross"] = wf["cross"]*F

    return wf