import numpy as np
import mpmath as mp
from scipy.special import gamma
from bilby.gw.source import *

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
    tMsun = 4.925490947e-6 #solar mass in seconds
    w = 2.*np.pi * (4*tMsun*ML) * f

    lague = np.array(vlaguerre(-0.5j*w, 0, 0.5j*w*y**2.),dtype='complex')
    ampl = np.power(-0.5j,1.+0.5j*w)*np.power(w,1.+0.5j*w)*gamma(-0.5j*w)*lague
    ampl[w==0. + 0.j] = 1. #gamma of 0 is inf but F(0,b)=1.
    return ampl

def pointlens(frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
    phi_12, a_2, tilt_2, phi_jl, theta_jn, phase,ML,y, **kwargs
    ):
    
    """
    Generate the frequency domain waveform for lensed image of a binary black hole
    
    Delta_phase is the phase difference w.r.t. the unlensed waveform
    """

    frequency_domain_source_model = lal_binary_black_hole

    # Actually generate the waveform by calling the generator
    wf = frequency_domain_source_model(
        frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
        phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, **kwargs
    )

    F = F_pointlens(frequency_array,ML,y)

    wf["plus"] = wf["plus"]*F
    wf["cross"] = wf["cross"]*F

    return wf