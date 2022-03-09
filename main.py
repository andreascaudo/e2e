from .simulation import Configuration
from .tool import generic as plt
import numpy as np


def run(configuration: Configuration):
    # Unpack configuration
    acquisition = configuration.acquisition
    spectrograph = configuration.spectrograph

    # First step: generate flux
    #Angstrom, Ph/s/cm^2/A
    wavelength, flux = acquisition.sed.get_flux()

    plt("SED", wavelength, flux, [
        "wavelength [$\AA$]", "flux " + r"[Ph/s/cm$^2$/$\AA$]"])

    wavelength, flux = acquisition.sed.normalize()

    plt("SED Normalized", wavelength, flux, [
        "wavelength [$\AA$]", "flux " + r"[Ph/s/cm$^2$/$\AA$]"])
