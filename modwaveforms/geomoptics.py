import numpy as np
from bilby.gw.source import *

def one_image_BBH(frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
    phi_12, a_2, tilt_2, phi_jl, theta_jn, phase,Delta_phase, **kwargs
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

    F_one_images = np.exp(-1j*np.sign(frequency_array)*Delta_phase)

    wf["plus"] = wf["plus"]*F_one_images
    wf["cross"] = wf["cross"]*F_one_images

    return wf


def two_images_BBH(frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
    phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, mu_rel,Delta_t,Delta_phase, **kwargs
    ):
    
    """
    Generate the frequency domain waveform for a binary black hole system with two images

    mu_rel is the relative amplitude of the second image
    Delta_t is the time delay between the two images
    Delta_phase is the phase difference between the two images
    """

    frequency_domain_source_model = lal_binary_black_hole

    # Actually generate the waveform by calling the generator
    wf = frequency_domain_source_model(
        frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
        phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, **kwargs
    )

    F_two_images = (1.+np.sqrt(mu_rel)*np.exp(1j*(2*np.pi*frequency_array*Delta_t-np.sign(frequency_array)*Delta_phase)))

    wf["plus"] = wf["plus"]*F_two_images
    wf["cross"] = wf["cross"]*F_two_images

    return wf