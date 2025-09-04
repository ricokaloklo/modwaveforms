"Python module with lensing function for a point mass"
import numpy as np
from ..utils.constants import *

"Time delays"
#Geometric optics    
def t_delay_geom_plus(y):
    return (y**2. + 2. - y*np.sqrt(y**2 +4.))/4.-np.log(np.abs(y+np.sqrt(y**2+4.))/2.)
def t_delay_geom_minus(y):
    return (y**2 + 2. + y*np.sqrt(y**2 +4.))/4.-np.log(np.abs(y-np.sqrt(y**2+4.))/2.)
def DeltaT(y):
    return t_delay_geom_minus(y)-t_delay_geom_plus(y)

def t_ref(M_Lz):
    """Reference time delay for a point mass lens
    
    Parameters
    ----------
    M_L : float
        Lens mass in solar masses
    z_L : float
        Lens redshift
    """
    return 4 *  Gnewton * M_Lz*MSUN / (Clight**3)

def Delta_t(M_Lz,y):
    return DeltaT(y) * t_ref(M_Lz)

"Magnification"
def mu_plus(y):
    return 0.5 + (y**2. + 2.)/(2.*y*np.sqrt(y**2 + 4.))
def mu_minus(y):
    return 0.5 - (y**2. + 2.)/(2.*y*np.sqrt(y**2 + 4.))

def mu_rel(y):
    return abs(mu_minus(y)/mu_plus(y))