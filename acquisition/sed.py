from scipy import constants
from ..instrument import spectrograph


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
        temperature: float            # K
    ) -> None:
        super().__init__(sed_type, magnitude, magnitude_system)
        self.temperature = temperature

    def get_flux(self):
        print("BLACK BODY FLUX")
        # return 8 * constants.pi * constants.h * constants.c ** 2 / ((source / 10 ** 10) ** 5 * (
        #    np.exp(constants.h * constants.c / (constants.k * temperature * (source / 10 ** 10))) - 1))  # J/(s * m2 * m) or W / m^2 * m


class Powerlaw(Sed):
    def __init__(
        self,
        sed_type: str,
        magnitude: float,               # [-]
        magnitude_system: str,          # ["Vega" or "AB"]
        index: float            # K
    ) -> None:
        super().__init__(sed_type, magnitude, magnitude_system)
        self.index = index

    def get_flux(self):
        print("POWER LAW FLUX")


class Spectrum(Sed):
    def __init__(
        self,
        sed_type: str,
        magnitude: float,               # [-]
        magnitude_system: str,          # ["Vega" or "AB"]
        spectrum_file: float            # K
    ) -> None:
        super().__init__(sed_type, magnitude, magnitude_system)
        self.spectrum_file = spectrum_file

    def get_flux(self):
        print("SPECTRUM FLUX")
