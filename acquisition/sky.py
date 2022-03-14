import pathlib
from os import system, remove
from os.path import join, dirname, realpath
import uuid

import numpy as np
import pandas as pd
from astropy.table import Table

from ..tool import unit_converter

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

    system("skycalc_cli " +
           "-i '" + str(input) +
           "' -o '" + str(output) + "'")

    remove(input)
    return output


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
        {'lambda': transmission, 'radiance': flux})
    output_file_transmission = pd.DataFrame(
        {'lambda': wavelength, 'transmission': transmission})

    output_file_radiance.to_csv(
        radiance_file, header=False, index=False, sep=" ")
    output_file_transmission.to_csv(
        transmission_file, header=False, index=False, sep=" ")


class Radiance:
    def __init__(
        self,
        radiance_file: np.ndarray
    ) -> None:
        self.wavelength = unit_converter.wavelength(
            radiance_file[0], "nm", "A")  # [nm] -> [A]
        self.flux = radiance_file[1]


class Transimission:
    def __init__(
        self,
        transmission_file: np.ndarray
    ) -> None:
        self.wavelength = unit_converter.wavelength(
            transmission_file[0], "nm", "A")  # [nm] -> [A]
        self.transmission = transmission_file[1]  # [-]


class Sky:
    def __init__(
        self,
        airmass: str,           # [-]
        moon_fli: float,        # [days]
        pwv: float,             # [mm]
        seeing: str             # [arcsec]
    ) -> None:
        self.airmass = airmass
        self.moon_fli = moon_fli
        self.pwv = pwv
        self.seeing = seeing
        # Load spectrum

    def get_sky(self):
        transmission_filename = "Sky_tra_MFLI_" + \
            str(int(self.moon_fli)) + "_Airm_" + \
            str(self.airmass) + "_PWV_" + str(self.pwv) + ".dat"  # float ???

        transmission_file = pathlib.Path(
            join(dirname(dirname(realpath(__file__))), "data", "sky", transmission_folder), transmission_filename)

        radiance_filename = "Sky_rad_MFLI_" + \
            str(int(self.moon_fli)) + "_Airm_" + \
            str(self.airmass) + "_PWV_" + str(self.pwv) + ".dat"

        radiance_file = pathlib.Path(
            join(dirname(dirname(realpath(__file__))), "data", "sky", radiance_folder), radiance_filename)

        print(transmission_file)
        print(radiance_file)

        if transmission_file.is_file() and radiance_file.is_file():
            transmission = Transimission(np.loadtxt(transmission_file))
            radiance = Radiance(np.loadtxt(radiance_file))
        else:
            output = write_skycalc_file(str(self.airmass), str(self.pwv), str(
                int(self.moon_fli * (180 / 1))))
            output_skycalc = Table.read(output)
            remove(output)
            write_sky_file(output_skycalc, transmission_file, radiance_file)

        return transmission_file, radiance_file
