# coding=utf-8

import math
from mpmath import asec
import numpy as np

wl_seeing = np.arange(1000, 25000, 250)

# wl_seeing:angstrom, Slit_*:arcsec, Angstrom, ref_wl:angstrom (usually 5000 A) with seeing (ref_seeing), Telescope D in meters, Lzero:  wave-front outer-scale in meters


def get_slit_efficiency(wavelength, airmass, slit_width, slit_length, ref_seeing, fwhm_ins, tel_diameter, l_zero, index=-1./5, ref_wl=5000):
    seeing_lambda_vector = []
    fwhm_iq = []
    fwhm_iq_total = []

    for i in range(0, len(wl_seeing)):
        total_FWHM_ins = math.sqrt(fwhm_atm(
            ref_seeing, airmass, wl_seeing[i], ref_wl, index, tel_diameter, l_zero)**2 + fwhm_tel(wl_seeing[i], tel_diameter)**2 + fwhm_ins**2)
        total_FWHM = math.sqrt(fwhm_atm(
            ref_seeing, airmass, wl_seeing[i], ref_wl, index, tel_diameter, l_zero)**2 + fwhm_tel(wl_seeing[i], tel_diameter)**2)
        # total_FWHM = ref_seeing*(wl_seeing[i]/ref_wl)**index #OLD
        coeffA = (slit_width * (math.log(4) ** 0.5)) / \
            (total_FWHM)  # adimensionale
        coeffB = (slit_length * (math.log(4) ** 0.5)) / \
            (total_FWHM)  # adimensionale
        seeing_lambda = math.erf(0.5 * (2 ** (1. / 2)) * coeffA) * \
            math.erf(0.5 * (2 ** (1. / 2)) * coeffB)

        fwhm_iq = fwhm_iq + [total_FWHM]
        fwhm_iq_total = fwhm_iq_total + [total_FWHM_ins]
        seeing_lambda_vector = seeing_lambda_vector + [seeing_lambda]

    # Interpolare for wavelength range
    seeing_lambda_vector = np.interp(
        wavelength, wl_seeing, seeing_lambda_vector)
    fwhm_iq_total = np.interp(wavelength, wl_seeing, fwhm_iq_total)

    return seeing_lambda_vector, fwhm_iq_total


def get_pinhole_efficiency(wl_seeing, airmass, d_pinhole, ref_wl, ref_seeing, tel_diameter, l_zero, index):
    seeing_lambda_vector = []

    for i in range(0, len(wl_seeing)):
        total_FWHM = math.sqrt(fwhm_atm(
            ref_seeing, airmass, wl_seeing[i], ref_wl, index, tel_diameter, l_zero)**2 + fwhm_tel(wl_seeing[i], tel_diameter)**2)
        # total_FWHM = ref_seeing*(wl_seeing[i]/ref_wl)**index #OLD
        seeing_lambda = 1 - \
            (1 / math.exp((math.sqrt(np.log(2)) * (d_pinhole / total_FWHM))**2))
        seeing_lambda_vector = seeing_lambda_vector + [seeing_lambda]

    return seeing_lambda_vector  # Vector of pinhole efficiency


def fwhm_atm(seeing, airmass, lmbd, ref_wl, index, tel_diameter, l_zero):
    return seeing * airmass**0.6 * ((lmbd / ref_wl) ** index) * math.sqrt(1 + fKolb(tel_diameter, l_zero) * 2.183 * (r_zero(seeing, lmbd, ref_wl, airmass) / l_zero)**0.356)


def fKolb(tel_diameter, l_zero):
    return 1 / (1 + 300 * tel_diameter / l_zero) - 1


def r_zero(seeing, lmbd, ref_wl, airmass):  # Return Rzero in meters
    return 0.100 * seeing**-1 * (lmbd / ref_wl)**1.2 * airmass**-0.6


def fwhm_tel(lmbd, tel_diameter):  # FWHM returned in arcsec
    # lmbd in nm and tel_diameter in meters
    return 0.000212 * (lmbd / 10) / tel_diameter


def fwhm_cam(PSF_sigma, plate_scale):
    return 2.35 * PSF_sigma * plate_scale


def fwhm_iq_camera(wl_seeing, airmass, ref_wl, ref_seeing, tel_diameter, l_zero, index, PSF_sigma, plate_scale):
    FWHM_iq_total = []
    for i in range(0, len(wl_seeing)):
        FWHM_iq = math.sqrt(fwhm_atm(ref_seeing, airmass, wl_seeing[i], ref_wl, index, tel_diameter, l_zero)**2
                            + fwhm_tel(wl_seeing[i], tel_diameter)**2
                            + fwhm_cam(PSF_sigma, plate_scale)**2)
        FWHM_iq_total.append(FWHM_iq)

    return FWHM_iq_total

# P = 760 mm Hg
# T = 15°C


def sea_level_refractive_index(ll):
    return (64.328 + (29498.1 / (146 - ((1 / ll)**2))) + (255.4 / (41 - ((1 / ll)**2)))) * 10**-6

    #       La silla meteo monitor default values
    # Lambda, Pressure in mmHg, Temperature in C°, Water Vapor


def high_altitude_refractive_index(ll, P=575.6, T=7, f=None):
    # print(sea_level_refractive_index(ll))
    ha_ref_index = sea_level_refractive_index(
        ll) * (P * (1 + ((1.049 - 0.0157 * T) * (10**-6)) * P)) / (720.883 * (1 + 0.003661 * T))
    # print(ha_ref_index)
    if f is not None:
        ha_ref_index = ha_ref_index * 10**6 - \
            (((0.0624 - 0.000680 / (ll**2)) / (1 + 0.003661 * T)) * f)
    return ha_ref_index * 10**-6


def differential_Atmospheric_refraction(ll, ref_ll, airmass, T, P, RH):
    return 206264.8 * abs(high_altitude_refractive_index(ll, P, T, RH) - high_altitude_refractive_index(ref_ll, P, T, RH)) * math.tan(asec(airmass))


def fwhm_diff_atm(cut_off_limits, airmass, fwhm):
    ref_wl = 1.1  # um
    T = 10  # C°
    P = 577.547  # mmHg
    WP = 9.5  # mmHg
    # Ang to um
    lowest_wl = differential_Atmospheric_refraction(
        cut_off_limits[0] / 10000, ref_wl, airmass, T, P, WP)
    highest_wl = differential_Atmospheric_refraction(
        cut_off_limits[1] / 10000, ref_wl, airmass, T, P, WP)
    distortion = abs(lowest_wl - highest_wl)
    print("Distortion:")
    print(distortion)
    return math.sqrt(fwhm**2 + distortion**2)
