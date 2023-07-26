import pathlib
from os import system, remove
from os.path import join, dirname, realpath
import uuid

import numpy as np
import pandas as pd
from astropy.table import Table

import matplotlib.pyplot as plt

from ..tool import unit_converter
from ..tool import tools

transmission_folder = "transmission"
radiance_folder = "radiance"


def write_skycalc_file(airmass, pwv, moon_sun_sep):

    id = str(uuid.uuid4())

    input = join(dirname(dirname(realpath(__file__))), "data",
                 "sky", "skycalc", 'SkyCalc_input_' + id + '.txt')
    output = join(dirname(dirname(realpath(__file__))), "data",
                  "sky", "skycalc", 'SkyCalc_input_' + id + '_Out.fits')
    template = join(dirname(dirname(realpath(__file__))),
                    "data", "sky", "skycalc", 'SkyCalc_input.txt')

    # From template file, setup new input
    with open(input, 'w') as new:
        with open(template) as template:
            cnt = 1
            line = template.readline()
            while line:
                if cnt == 1:
                    # AIRMASS
                    # Convert the string to a list
                    line_new = "airmass         :  " + airmass + "\n"

                elif cnt == 5:
                    # PWV
                    # Convert the string to a list
                    line_new = "pwv             :  " + pwv + "\n"

                elif cnt == 8:
                    # MOON
                    # Convert the string to a list
                    line_new = "moon_sun_sep    :  " + moon_sun_sep + "\n"

                else:
                    line_new = line

                new.write(line_new)
                cnt = 1 + cnt
                line = template.readline()

        template.close()
    new.close()

    # Call to SkyCalc -i input -o output
    system("skycalc_cli " +
           "-i '" + str(input) +
           "' -o '" + str(output) + "'")

    remove(input)  # Remove input file
    return output


# Write radiance.dat and transmission.dat files
def write_sky_file(output, transmission_file, radiance_file):

    len_output = len(output)

    wavelength = np.zeros(len_output)
    flux = np.zeros(len_output)
    transmission = np.zeros(len_output)

    for i in range(0, len_output):
        wavelength[i] = output[i][0]  # [nm]
        flux[i] = output[i][1]  # Radiance
        transmission[i] = output[i][4]  # Transmission

    output_file_radiance = pd.DataFrame(
        {'lambda': wavelength, 'radiance': flux})
    output_file_transmission = pd.DataFrame(
        {'lambda': wavelength, 'transmission': transmission})

    output_file_radiance.to_csv(
        radiance_file, header=False, index=False, sep=" ")
    output_file_transmission.to_csv(
        transmission_file, header=False, index=False, sep=" ")


class Radiance:
    def __init__(
        self,
        radiance_file: np.ndarray = None,
        slit_dimension: list = None
    ) -> None:
        if radiance_file is not None:
            self.wavelength = unit_converter.wavelength(
                radiance_file.T[0], "nm", "A")  # [nm] -> [A]
        else:
            # I Changed before was np.arange(300, 24000.1, 0.1)
            self.wavelength = None

        if slit_dimension is not None:
            self.flux = unit_converter.flux(
                radiance_file.T[1], "photons/m^2/s/µm/asec^2", "photons/cm^2/s/A", slit_dimension)  # [Ph/s/m2/µm/asec2] -> [Ph/s/cm2/A]
        else:
            # I Changed before was np.zeros(len(self.wavelength))
            self.flux = None


class Transimission:
    def __init__(
        self,
        transmission_file: np.ndarray = None
    ) -> None:
        if transmission_file is not None:
            self.wavelength = unit_converter.wavelength(
                transmission_file.T[0], "nm", "A")  # [nm] -> [A]
            self.transmission = transmission_file.T[1]  # [-]
        else:
            self.wavelength = np.arange(300, 24000.1, 0.1)
            self.transmission = np.ones(len(self.wavelength))


class Sky:
    def __init__(
        self,
        airmass: str,           # [-]
        moon_fli: float,        # [days]
        pwv: float,             # [mm]
        seeing: str,             # [arcsec]
        active: bool = True
    ) -> None:
        self.airmass = airmass
        self.moon_fli = moon_fli
        self.pwv = pwv
        self.seeing = seeing
        self.active = active
        # Load spectrum

    def get_sky(self, slit_x, slit_y):
        if self.active:
            transmission_filename = "Sky_tra_MFLI_" + \
                str(self.moon_fli) + "_Airm_" + \
                str(self.airmass) + "_PWV_" + \
                str(self.pwv) + ".dat"  # float ???

            transmission_file = pathlib.Path(
                join(dirname(dirname(realpath(__file__))), "data", "sky", transmission_folder), transmission_filename)

            radiance_filename = "Sky_rad_MFLI_" + \
                str(self.moon_fli) + "_Airm_" + \
                str(self.airmass) + "_PWV_" + str(self.pwv) + ".dat"

            radiance_file = pathlib.Path(
                join(dirname(dirname(realpath(__file__))), "data", "sky", radiance_folder), radiance_filename)

            print(transmission_file)
            print(radiance_file)

            # If files do not exist -> call skycal
            if not (transmission_file.is_file() and radiance_file.is_file()):
                output = write_skycalc_file(str(self.airmass), str(self.pwv), str(
                    np.round(self.moon_fli * (180 / 1), 2)))
                output_skycalc = Table.read(output)
                remove(output)
                write_sky_file(
                    output_skycalc, transmission_file, radiance_file)

            self.transmission = Transimission(np.loadtxt(transmission_file))
            self.radiance = Radiance(np.loadtxt(
                radiance_file), [slit_x, slit_y])

            return self.transmission, self.radiance
        else:
            self.transmission = Transimission()
            self.radiance = Radiance()

            return self.transmission, self.radiance

    def set_efficiency(self, wavematrix, efficiency):
        self.wavelength_matrix = wavematrix
        size_matrix = self.wavelength_matrix.shape
        self.transmission.transmission_matrix = np.zeros(size_matrix)

        for i in range(0, size_matrix[0]):
            self.transmission.transmission_matrix[i] = tools.interp(
                self.transmission.wavelength, self.transmission.transmission, self.wavelength_matrix[i], fill_value="extrapolate")

        self.sky_efficiency_matrix = efficiency
