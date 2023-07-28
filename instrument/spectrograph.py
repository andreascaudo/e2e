import numpy as np
import scipy.io
import math
from . import grating as grating_obj
from .grating import Grating
from .image_slicer import Slice
from ..tool import unit_converter
from ..tool import tools


class Spectrograph:
    def __init__(
        self,
        name: str,                      # Spectrograph Name
        arm: str,                       # Spectrograph Arm

        # Slicing
        slit_size_x: list,              # [arcsec]
        slit_size_y: list,              # [arcsec]
        slit_size_x_calibration: list,  # cm
        slit_size_y_calibration: list,  # cm

        # Grating
        grating: Grating,

        # Detector
        n_pixels: int,                  # [-]
        dimension_pixel: float,         # [um]

        # Efficiency
        telescope_fdr_file: str,
        instrument_fdr_file: str,
        psf_map_file: str,

        # Non-mandatory parameters
        image_slicer: dict = None

    ) -> None:
        self.name = name
        self.arm = arm

        grating_type = getattr(grating_obj, grating["type"])
        self.grating = grating_type(**grating)

        self.telescope_fdr_file = load_fdr(telescope_fdr_file)
        self.instrument_fdr_file = load_fdr(instrument_fdr_file)
        self.psf_map_file = load_fdr(psf_map_file)

        self.slit_size_x = slit_size_x
        self.slit_size_y = slit_size_y
        self.slit_size_x_calibration = slit_size_x_calibration
        self.slit_size_y_calibration = slit_size_y_calibration

        self.wavematrix, self.len_n_orders, self.n_p = self.grating.get_efficiency()

        # GRATING
        self.telescope_fdr = np.zeros((self.len_n_orders, self.n_p))
        self.instrument_fdr = np.zeros((self.len_n_orders, self.n_p))

        for i in range(0, self.len_n_orders):
            self.telescope_fdr[i] = tools.interp(
                self.telescope_fdr_file.T[0], self.telescope_fdr_file.T[1], self.wavematrix[i], fill_value="extrapolate")

            self.instrument_fdr[i] = tools.interp(
                self.instrument_fdr_file.T[0], self.instrument_fdr_file.T[1], self.wavematrix[i], fill_value="extrapolate")

        # TBD: When implementing CUBES remember to check psf map vs MATLAB
        if self.psf_map_file[0][0] < self.psf_map_file[-1][0]:
            self.psf_map = np.flip(self.psf_map_file)
        else:
            self.psf_map = self.psf_map_file

        self.wavematrix = unit_converter.wavelength(self.wavematrix, "um", "A")

        self.telescope_spectrograph_efficiency_fdr = self.telescope_fdr * self.instrument_fdr

        self.n_pixels = n_pixels
        self.dimension_pixel = dimension_pixel

        # Genrate a set of slices based on two cases:
        # 1. If the user provides an image slicer, then use those slices
        # 2. If the user does not provide an image slicer, then use one slice per order
        self.slices = []
        if image_slicer != None:
            for i in range(0, image_slicer["n_slices"]):
                self.slices.append(Slice(i+1, image_slicer["shift_arc"][i]))
        else:
            self.slices.append(Slice(1, None))

        self.len_n_slices = len(self.slices)

    def set_subpixels(self, pixel_oversampling, psf_map_pixel_number):
        self.n_pixels_subpixel = self.n_pixels * pixel_oversampling
        self.psf_map_pixel_number_subpixel = psf_map_pixel_number * pixel_oversampling
        self.subpixel_edge = self.psf_map_pixel_number_subpixel * 2
        self.detector_subpixel = np.zeros(
            (310*pixel_oversampling, self.n_pixels_subpixel + self.subpixel_edge, self.len_n_orders,  self.len_n_slices))


def load_fdr(file):
    try:
        if file.endswith(".mat"):
            return scipy.io.loadmat(file)
        elif file.endswith(".npy"):
            return np.load(file, allow_pickle=True)
        else:
            return np.loadtxt(file, delimiter=" ")
    except Exception as e:
        print(e)
        raise "ERROR: " + e
