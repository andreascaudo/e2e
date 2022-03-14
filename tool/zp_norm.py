# coding=utf-8

import os
from os.path import dirname, join
from scipy import integrate
import numpy as np

filters_folder = "filters"
filters_folder = join(
    dirname(dirname(os.path.realpath(__file__))), "data", filters_folder)

band_filters = {}

for filename in os.listdir(filters_folder):
    if filename.endswith(".txt"):
        sed_txt = np.loadtxt(os.path.join(filters_folder, filename),
                             comments='#')
        band_filters[filename] = sed_txt


def get_zp_norm(wave_strumento, nfotoni, filer_name, zeropoint, lambda_0, mag):
    ind = np.argwhere(wave_strumento >= lambda_0)
    ind1 = ind[0][0]
    F0 = nfotoni[ind1]
    Zp_scalato = zeropoint * (2.512) ** (-mag)
    K = Zp_scalato / F0

    return nfotoni * K  # [phot/s/cm2/A]


def get_filter_norm(wave_strumento, nfotoni, filer_name, zeropoint, mag):
    try:
        bf = band_filters["phot_" + filer_name + ".txt"]
    except Exception as e:
        print(e)
        print("Band filter file not found!")

    bf_w = bf.T[0] * 10  # Å
    bf_f = bf.T[1]  # Unitless

    bf_f_scaled = np.interp(bf_w, wave_strumento, nfotoni)
    bf_f_scaled = bf_f_scaled * bf_f  # phot/s/cm^2/ang

    integral_flux = integrate.trapz(bf_f_scaled, bf_w)
    integral_filt = integrate.trapz(bf_f, bf_w)

    Zp_scalato = (zeropoint*integral_filt) * (2.512) ** (-mag)
    K = Zp_scalato/integral_flux

    return nfotoni * K  # [phot/s/cm2/A]
