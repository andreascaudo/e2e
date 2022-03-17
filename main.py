from .simulation import Configuration
from .tool import efficiency
from .tool import generic as plt
import numpy as np


def run(configuration: Configuration):
    # Unpack configuration
    telescope = configuration.telescope
    spectrograph = configuration.spectrograph
    acquisition = configuration.acquisition

    # TBI: Implement a function TO CHECK if acquisition reflects the spectrograph parameters

    # First step: generate flux
    # Angstrom, Ph/s/cm^2/A
    sed_wavelength, sed_flux = acquisition.sed.get_flux()

    plt("SED", sed_wavelength, sed_flux, [
        "wavelength [$\AA$]", "flux " + r"[Ph/s/cm$^2$/$\AA$]"])

    sed_wavelength, sed_flux = acquisition.sed.normalize()

    plt("SED Normalized", sed_wavelength, sed_flux, [
        "wavelength [$\AA$]", "flux " + r"[Ph/s/cm$^2$/$\AA$]"])

    # Get Slit Efficiency and Image Quality
    slit_efficiency, fwhm_iq = efficiency.get_slit_efficiency(sed_wavelength, acquisition.sky.airmass,
                                                              acquisition.characteristics.slit_size_x, acquisition.characteristics.slit_size_y,
                                                              acquisition.sky.seeing, spectrograph.fwhm_instrument, telescope.diameter, telescope.l_zero)

    plt("Slit efficiecny", sed_wavelength,
        slit_efficiency, ["wavelength [$\AA$]", "[-]"])
    plt("FWHM Image Quality", sed_wavelength,
        fwhm_iq,  ["wavelength [$\AA$]", "[arcsec]"])

    sky_wavelength, sky_flux = acquisition.sky.get_sky()
