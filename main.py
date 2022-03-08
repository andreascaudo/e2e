from .simulation import Configuration
import matplotlib.pyplot as plt
from scipy import constants
import numpy as np


def run(configuration: Configuration):
    # Unpack configuration
    output = configuration.output
    acquisition = configuration.acquisition
    telescope = configuration.telescope
    spectrograph = configuration.spectrograph

    # First step: generate flux
    wavelength, flux = acquisition.sed.get_flux(spectrograph.wavelength_band)

    plt.plot(wavelength, flux)
    plt.show()
