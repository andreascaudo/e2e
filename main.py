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

    # 1st step: generate flux
    # [Angstrom], [Ph/s/cm^2/A]
    sed_wavelength, sed_flux = acquisition.sed.get_flux()

    plt("SED", sed_wavelength, sed_flux, [
        "wavelength [$\AA$]", "flux " + r"[Ph/s/cm$^2$/$\AA$]"])

    sed_wavelength, sed_flux = acquisition.sed.normalize()

    plt("SED Normalized", sed_wavelength, sed_flux, [
        "wavelength [$\AA$]", "flux " + r"[Ph/s/cm$^2$/$\AA$]"])

    # 2nd step: get sky radiance and transission
    # [-] ,[ph/s/cm2/A]
    transimission, radiance = acquisition.sky.get_sky(
        acquisition.characteristics.slit_size_x, acquisition.characteristics.slit_size_y)

    plt("Sky transmission", transimission.wavelength,
        transimission.transmission, ["wavelength [$\AA$]", "[-]"])

    plt("Sky radiance", radiance.wavelength,
        radiance.flux, ["wavelength [$\AA$]", "[ph/s/cm2/A]"])

    # 3th step: Get Slit Efficiency and Image Quality
    slit_efficiency, fwhm_iq = efficiency.get_slit_efficiency(spectrograph.wavelength_band, acquisition.sky.airmass,
                                                              acquisition.characteristics.slit_size_x, acquisition.characteristics.slit_size_y,
                                                              acquisition.sky.seeing, spectrograph.fwhm_instrument, (telescope.diameter/100), telescope.l_zero)

    plt("Slit efficiecny", spectrograph.wavelength_band,
        slit_efficiency, ["wavelength [$\AA$]", "[-]"])
    plt("FWHM Image Quality", spectrograph.wavelength_band,
        fwhm_iq,  ["wavelength [$\AA$]", "[arcsec]"])

    # 4th step: Get Instrument-Spectrograph Efficiency
