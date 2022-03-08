import numpy as np
from scipy import constants
from ..instrument import spectrograph
from ..tool import unit_converter


class Sed:
    def __init__(
        self,
        sed_type: str,
        magnitude: float,               # [-]
        magnitude_system: str           # ["Vega" or "AB"]
    ) -> None:
        self.sed_type = sed_type
        self.magnitude = magnitude
        self.magnitude_system = magnitude_system


class Blackbody(Sed):
    def __init__(
        self,
        sed_type: str,
        magnitude: float,               # [-]
        magnitude_system: str,          # ["Vega" or "AB"]
        temperature: float              # K
    ) -> None:
        super().__init__(sed_type, magnitude, magnitude_system)
        self.temperature = temperature

    def get_flux(self, wavelength_band):
        self.wavelength = wavelength_band
        self.flux = 8 * constants.pi * constants.h * constants.c ** 2 / ((wavelength_band / 10 ** 10) ** 5 * (
            np.exp(constants.h * constants.c / (constants.k * self.temperature * (wavelength_band / 10 ** 10))) - 1))  # J/(s * m2 * m) or W / m^2 * m

        self.flux = unit_converter.flux(self.flux, "J/(s * m^3)",
                                        "photons/cm^2/s/A", self.wavelength)

        return self.wavelength, self.flux


class Powerlaw(Sed):
    def __init__(
        self,
        sed_type: str,
        magnitude: float,               # [-]
        magnitude_system: str,          # ["Vega" or "AB"]
        index: float
    ) -> None:
        super().__init__(sed_type, magnitude, magnitude_system)
        self.index = index

    def get_flux(self, wavelength_band):
        self.wavelength = wavelength_band
        self.flux = np.power(wavelength_band, self.index)  # J/(s * m2 * m)

        self.flux = unit_converter.flux(self.flux, "J/(s * m^3)",
                                        "photons/cm^2/s/A", self.wavelength)

        return self.wavelength, self.flux


class Spectrum(Sed):
    def __init__(
        self,
        sed_type: str,
        magnitude: float,               # [-]
        magnitude_system: str,          # ["Vega" or "AB"]
        spectrum_file: float
    ) -> None:
        super().__init__(sed_type, magnitude, magnitude_system)
        self.spectrum_file = spectrum_file

    def get_flux(self, wavelength_band):
        try:
            sed_file = np.loadtxt(self.spectrum_file, comments='#')
        except Exception as e:
            print(e)

        self.wavelength = wavelength_band
        self.flux = np.interp(wavelength_band, sed_file.T[0], sed_file.T[1])

        self.flux = unit_converter.flux(self.flux, "ergs/cm^2/s/A",
                                        "photons/cm^2/s/A", self.wavelength)

        return self.wavelength, self.flux
