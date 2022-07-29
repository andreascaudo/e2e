import numpy as np
from scipy import constants
from ..instrument import spectrograph
from ..tool import unit_converter
from ..tool import zp_norm
from ..tool import magnitude
from ..tool import tools


class Sed:
    def __init__(
        self,
        sed_type: str
    ) -> None:
        self.sed_type = sed_type

    def normalize(self):
        if self.magnitude_system == "AB":
            self.magnitude = magnitude.ab_to_vega_converter(
                self.magnitude, self.band)
            self.magnitude_system = "Vega"

        lambda_0, zeropoint = magnitude.get_vega_flux_zeropoints(
            self.band, "PHll")

        if self.bandpass_normalization:
            self.flux = zp_norm.get_filter_norm(
                self.wavelength, self.flux, self.band, zeropoint, self.magnitude)  # [phot/s/cm2/A]
        else:
            self.flux = zp_norm.get_zp_norm(
                self.wavelength, self.flux, self.band, zeropoint, lambda_0, self.magnitude)  # [phot/s/cm2/A]

        return self.wavelength, self.flux

    def set_efficiency(self, transmission_matrix, wavematrix, efficiency):
        self.wavelength_matrix = wavematrix
        self.sed_total_efficincy = transmission_matrix * efficiency


class Blackbody(Sed):
    def __init__(
        self,
        sed_type: str,
        band: str,
        magnitude: float,               # [-]
        magnitude_system: str,          # ["Vega" or "AB"]
        temperature: float,              # K
        bandpass_normalization: bool = True
    ) -> None:
        super().__init__(sed_type)
        self.band = band
        self.magnitude = magnitude
        self.magnitude_system = magnitude_system
        self.bandpass_normalization = bandpass_normalization
        self.temperature = temperature
        self.calibration = False
        self.wavelength = np.arange(1000, 25000, 1)  # Angstrom

    def get_flux(self):
        self.flux = 8 * constants.pi * constants.h * constants.c ** 2 / ((self.wavelength / 10 ** 10) ** 5 * (
            np.exp(constants.h * constants.c / (constants.k * self.temperature * (self.wavelength / 10 ** 10))) - 1))  # J/(s * m2 * m) or W / m^2 * m

        self.flux = unit_converter.flux(self.flux, "J/(s * m^3)",
                                        "photons/cm^2/s/A", self.wavelength)

        return self.wavelength, self.flux


class Powerlaw(Sed):
    def __init__(
        self,
        sed_type: str,
        band: str,
        magnitude: float,               # [-]
        magnitude_system: str,          # ["Vega" or "AB"]
        index: float,
        bandpass_normalization: bool = True
    ) -> None:
        super().__init__(sed_type)
        self.band = band
        self.magnitude = magnitude
        self.magnitude_system = magnitude_system
        self.bandpass_normalization = bandpass_normalization
        self.index = index
        self.calibration = False
        self.wavelength = np.arange(1000, 25000, 1)  # Angstrom

    def get_flux(self):
        self.flux = np.power(self.wavelength, self.index)  # J/(s * m2 * m)

        self.flux = unit_converter.flux(self.flux, "J/(s * m^3)",
                                        "photons/cm^2/s/A", self.wavelength)

        return self.wavelength, self.flux


class Flat(Sed):
    def __init__(
        self,
        sed_type: str,
        energy: float
    ) -> None:
        super().__init__(sed_type)
        self.energy = energy
        self.calibration = True
        self.wavelength = np.arange(1000, 25000, 1)  # Angstrom

    def get_flux(self):
        self.flux = np.full(len(self.wavelength), self.energy)

        self.flux = unit_converter.flux(self.flux, "ergs/cm^2/s/A",
                                        "photons/cm^2/s/A", self.wavelength)

        return self.wavelength, self.flux


class Lamp(Sed):
    def __init__(
        self,
        sed_type: str,
        spectrum_file: float
    ) -> None:
        super().__init__(sed_type)
        self.spectrum_file = spectrum_file
        self.calibration = True

    def get_flux(self):
        try:
            sed_file = np.loadtxt(self.spectrum_file, comments='#')
        except Exception as e:
            print(e)

        delta_lambda = tools.myround(sed_file[1, 0] - sed_file[0, 0])
        self.wavelength = np.arange(1000, 25000, delta_lambda)  # Angstrom

        self.flux = np.interp(
            self.wavelength, sed_file.T[0], sed_file.T[1])  # 0: Wavelength [A], 1: Flux [ergs/s/cm^2/A]

        self.flux = unit_converter.flux(self.flux, "ergs/cm^2/s/A",
                                        "photons/cm^2/s/A", self.wavelength)

        return self.wavelength, self.flux


class Spectrum(Sed):
    def __init__(
        self,
        sed_type: str,
        band: str,
        magnitude: float,               # [-]
        magnitude_system: str,          # ["Vega" or "AB"]
        spectrum_file: float,
        bandpass_normalization: bool = True
    ) -> None:
        super().__init__(sed_type)
        self.band = band
        self.magnitude = magnitude
        self.magnitude_system = magnitude_system
        self.bandpass_normalization = bandpass_normalization
        self.spectrum_file = spectrum_file
        self.calibration = False

    def get_flux(self):
        try:
            sed_file = np.loadtxt(self.spectrum_file, comments='#')
        except Exception as e:
            print(e)

        delta_lambda = tools.myround(sed_file[1, 0] - sed_file[0, 0])
        self.wavelength = np.arange(1000, 25000, delta_lambda)  # Angstrom

        self.flux = np.interp(
            self.wavelength, sed_file.T[0], sed_file.T[1])  # 0: Wavelength [A], 1: Flux [ergs/s/cm^2/A]

        self.flux = unit_converter.flux(self.flux, "ergs/cm^2/s/A",
                                        "photons/cm^2/s/A", self.wavelength)

        return self.wavelength, self.flux
