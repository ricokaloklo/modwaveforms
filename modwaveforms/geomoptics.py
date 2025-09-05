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

    # complex conjugate because of different sign convention in the Fourier transform of 
    # LALSuite (engineering) compared to the one used in wave optics (physics)
    F_one_images = np.conjugate(F_one_images)

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

    # complex conjugate because of different sign convention in the Fourier transform of 
    # LALSuite (engineering) compared to the one used in wave optics (physics)
    F_two_images = np.conjugate(F_two_images)

    wf["plus"] = wf["plus"]*F_two_images
    wf["cross"] = wf["cross"]*F_two_images

    return wf

def fold_caustic(frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
    phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, Delta_t, positive_phase, **kwargs
    ):
    
    """
    Generate the frequency domain waveform for a binary black hole system with two images

    Delta_t is the time delay between the two images
    positive_phase is the phase of the positive parity image
    """

    frequency_domain_source_model = lal_binary_black_hole

    # Actually generate the waveform by calling the generator
    wf = frequency_domain_source_model(
        frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
        phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, **kwargs
    )

    F_fold_two_images = positive_phase*(1. + np.exp(1j*(2*np.pi*frequency_array*Delta_t-positive_phase*np.sign(frequency_array)*np.pi/2)))

    # complex conjugate because of different sign convention in the Fourier transform of 
    # LALSuite (engineering) compared to the one used in wave optics (physics)
    F_fold_two_images = np.conjugate(F_fold_two_images)

    wf["plus"] = wf["plus"]*F_fold_two_images
    wf["cross"] = wf["cross"]*F_fold_two_images

    return wf

def cusp_caustic(frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
    phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, Delta_t_10, Delta_t_20, mu_rel, positive_phase, **kwargs
    ):
    
    """
    Generate the frequency domain waveform for a binary black hole system with two images

    Delta_t is the time delay between the two images
    positive_phase is the phase of the positive parity image
    """

    frequency_domain_source_model = lal_binary_black_hole

    # Actually generate the waveform by calling the generator
    wf = frequency_domain_source_model(
        frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
        phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, **kwargs
    )

    F_cusp_three_images = positive_phase*(1. + np.sqrt(abs(mu_rel))*np.exp(1j*(2*np.pi*frequency_array*Delta_t_10-positive_phase*np.sign(frequency_array)*np.pi/2))
                                          + np.sqrt(1-abs(mu_rel))*np.exp(1j*(2*np.pi*frequency_array*Delta_t_20-positive_phase*np.sign(frequency_array)*np.pi/2)))

    # complex conjugate because of different sign convention in the Fourier transform of 
    # LALSuite (engineering) compared to the one used in wave optics (physics)
    F_cusp_three_images = np.conjugate(F_cusp_three_images)
    
    wf["plus"] = wf["plus"]*F_cusp_three_images
    wf["cross"] = wf["cross"]*F_cusp_three_images

    return wf