from scipy import interpolate
from . import unit_converter
import numpy as np


def interp(x, y, new_x, fill_value, kind="linear"):
    function = interpolate.interp1d(x, y, kind=kind, fill_value=fill_value)
    return function(new_x)


def integration(lam, delta_lam, flux):
    l1 = lam  # [m]
    l2 = lam + delta_lam  # [m]
    dl = delta_lam / 100  # [m]

    lam_vect = np.arange(l1, l2, dl)

    # print(len(lam_vect))

    spec_flux_norm_interp = interp(unit_converter.wavelength(
        flux.wavelength, "A", "m"), flux.flux, lam_vect, fill_value="extrapolate")  # flux in [phot/s/cm^2/ang] and lambda in [m]

    # E=N*h*nu -> N=E*lambda/h*c
    dl = unit_converter.wavelength(dl, "m", "A")
    counts = 0

    for i in range(0, len(spec_flux_norm_interp)):
        counts = counts + \
            spec_flux_norm_interp[i] * dl  # Nfot/(s*cm^2)

    return counts
