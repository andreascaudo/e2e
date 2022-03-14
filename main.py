from .simulation import Configuration
from .tool import generic as plt
import numpy as np


def run(configuration: Configuration):
    # Unpack configuration
    acquisition = configuration.acquisition
    spectrograph = configuration.spectrograph

    # First step: generate flux
    #Angstrom, Ph/s/cm^2/A
    sed_wavelength, sed_flux = acquisition.sed.get_flux()

    plt("SED", sed_wavelength, sed_flux, [
        "wavelength [$\AA$]", "flux " + r"[Ph/s/cm$^2$/$\AA$]"])

    sed_wavelength, sed_flux = acquisition.sed.normalize()

    plt("SED Normalized", sed_wavelength, sed_flux, [
        "wavelength [$\AA$]", "flux " + r"[Ph/s/cm$^2$/$\AA$]"])

    sky_wavelength, sky_flux = acquisition.sky.get_sky()
