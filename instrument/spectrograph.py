from re import sub
import numpy as np
import scipy.io
import math
from . import grating as grating_obj
from .grating import Grating
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

        wavelength_min: float,          # [A]
        wavelength_max: float,          # [A]

        # Efficiency
        telescope_fdr_file: str,
        commonpath_ir_fdr_file: str,
        commonpath_vis_fdr_file: str,
        collimator_fdr_file: str,
        field_mirror_fdr_file: str,
        cross_disperser_fdr_file: str,
        fold_mirror_fdr_file: str,
        camera_fdr_file: str,
        qe_detector_file: str,
        psf_map_file: str


    ) -> None:
        self.name = name
        self.arm = arm

        grating_type = getattr(grating_obj, grating["type"])
        self.grating = grating_type(**grating)

        self.telescope_fdr_file = load_fdr(telescope_fdr_file)
        self.commonpath_ir_fdr_file = load_fdr(commonpath_ir_fdr_file)
        self.commonpath_vis_fdr_file = load_fdr(commonpath_vis_fdr_file)
        self.collimator_fdr_file = load_fdr(collimator_fdr_file)
        self.field_mirror_fdr_file = load_fdr(field_mirror_fdr_file)
        self.cross_disperser_fdr_file = load_fdr(cross_disperser_fdr_file)
        self.fold_mirror_fdr_file = load_fdr(fold_mirror_fdr_file)
        self.camera_fdr_file = load_fdr(camera_fdr_file)
        self.qe_detector_file = load_fdr(qe_detector_file)
        self.psf_map_file = load_fdr(psf_map_file)

        self.slit_size_x = slit_size_x
        self.slit_size_y = slit_size_y
        self.slit_size_x_calibration = slit_size_x_calibration
        self.slit_size_y_calibration = slit_size_y_calibration

        self.wavematrix, self.b_phase, self.b, self.len_n_orders, self.n_p = self.grating.get_efficiency()

        # CommonPath IR
        wavetemp = math.floor(self.wavematrix[0][0]*1000)
        index_start = np.where(self.commonpath_vis_fdr_file == wavetemp)[0]
        index_end = np.where(self.commonpath_vis_fdr_file == 849)[0]
        # Extrap of data from UV-VIS
        CPVISnofilt_2_800 = self.commonpath_vis_fdr_file[index_start[0]
            :index_end[0]+1].T[1]
        CPNIRnofilt_2_800 = 1 - (CPVISnofilt_2_800/100)
        # Re definition of CPIR_fdr
        vect_800_849 = np.arange((wavetemp/1000), 0.850, 0.001)
        CPIR_fdr_800 = np.c_[vect_800_849, CPNIRnofilt_2_800]
        CPIR_fdr = np.concatenate(
            (CPIR_fdr_800, self.commonpath_ir_fdr_file))

        # GRATING
        self.grating_fdr = np.zeros((self.len_n_orders, self.n_p))
        self.grating_rsc = np.zeros((self.len_n_orders, self.n_p))

        self.telescope_fdr = np.zeros((self.len_n_orders, self.n_p))
        self.commonpath_ir_fdr = np.zeros((self.len_n_orders, self.n_p))
        self.collimator_fdr = np.zeros((self.len_n_orders, self.n_p))
        self.field_mirror_fdr = np.zeros((self.len_n_orders, self.n_p))
        self.cross_disperser_fdr = np.zeros((self.len_n_orders, self.n_p))
        self.fold_mirror_fdr = np.zeros((self.len_n_orders, self.n_p))
        self.camera_fdr = np.zeros((self.len_n_orders, self.n_p))
        self.qe_detector = np.zeros((self.len_n_orders, self.n_p))

        for i in range(0, self.len_n_orders):
            # Rescale grating efficiency
            self.grating_fdr[i] = tools.interp(
                self.grating.grating_lambda, self.grating.grating_efficiency, self.wavematrix[i], fill_value="extrapolate")
            self.grating_fdr[i] = self.grating_fdr[i] * self.b[i]

            self.telescope_fdr[i] = tools.interp(
                self.telescope_fdr_file.T[0], self.telescope_fdr_file.T[1], self.wavematrix[i], fill_value="extrapolate")
            self.commonpath_ir_fdr[i] = tools.interp(
                CPIR_fdr.T[0], CPIR_fdr.T[1], self.wavematrix[i], fill_value="extrapolate", kind="nearest")
            self.collimator_fdr[i] = tools.interp(
                self.collimator_fdr_file.T[0], self.collimator_fdr_file.T[1], self.wavematrix[i], fill_value="extrapolate", kind="nearest")
            self.field_mirror_fdr[i] = tools.interp(
                self.field_mirror_fdr_file.T[0], self.field_mirror_fdr_file.T[1], self.wavematrix[i], fill_value="extrapolate", kind="nearest")
            self.cross_disperser_fdr[i] = tools.interp(
                self.cross_disperser_fdr_file.T[0], self.cross_disperser_fdr_file.T[1], self.wavematrix[i], fill_value="extrapolate", kind="nearest")
            self.fold_mirror_fdr[i] = tools.interp(
                self.fold_mirror_fdr_file.T[0], self.fold_mirror_fdr_file.T[1], self.wavematrix[i], fill_value="extrapolate", kind="nearest")
            self.camera_fdr[i] = tools.interp(self.camera_fdr_file.T[0], self.camera_fdr_file.T[1],
                                              self.wavematrix[i], fill_value="extrapolate", kind="nearest")
            self.qe_detector[i] = tools.interp(self.qe_detector_file.T[0], self.qe_detector_file.T[1],
                                               self.wavematrix[i], fill_value="extrapolate", kind="nearest")

        self.psf_map = np.flip(self.psf_map_file["PSFmap_Struct"][0])

        '''
        for len_i in range(0, self.len_n_orders):
            plt.plot(self.wavematrix[len_i], self.b[len_i])

        plt.show()

        for len_i in range(0, self.len_n_orders):
            plt.plot(self.wavematrix[len_i], self.grating_fdr[len_i])

        plt.show()
        '''

        self.wavematrix = unit_converter.wavelength(self.wavematrix, "um", "A")

        self.telescope_spectrograph_efficiency_fdr = self.telescope_fdr * self.commonpath_ir_fdr * self.collimator_fdr * \
            self.field_mirror_fdr * self.grating_fdr * self.cross_disperser_fdr * \
            self.fold_mirror_fdr * self.camera_fdr * self.qe_detector

        self.spectrograph_efficiency_fdr = self.commonpath_ir_fdr * self.collimator_fdr * self.field_mirror_fdr * \
            self.grating_fdr * self.cross_disperser_fdr * \
            self.fold_mirror_fdr * self.camera_fdr * self.qe_detector

        '''
        for len_i in range(0, self.len_n_orders):
            plt.plot(self.wavematrix[len_i],
                     self.telescope_spectrograph_efficiency_fdr[len_i])

        plt.show()

        for len_i in range(0, self.len_n_orders):
            plt.plot(self.wavematrix[len_i],
                     self.spectrograph_efficiency_fdr[len_i])

        plt.show()

        '''
        self.n_pixels = n_pixels
        self.dimension_pixel = dimension_pixel

        self.wavelength_min = wavelength_min
        self.wavelength_max = wavelength_max
        self.wavelength_band = np.arange(wavelength_min, wavelength_max, 1)

    def set_subpixels(self, pixel_oversampling, psf_map_pixel_number):
        self.n_pixels_subpixel = self.n_pixels * pixel_oversampling
        self.psf_map_pixel_number_subpixel = psf_map_pixel_number * pixel_oversampling
        self.subpixel_edge = self.psf_map_pixel_number_subpixel * 2
        self.detector_subpixel = np.zeros(
            (310*pixel_oversampling, self.n_pixels_subpixel + self.subpixel_edge, self.len_n_orders))


def load_fdr(file):
    try:
        if file.endswith(".mat"):
            return scipy.io.loadmat(file)
        else:
            return np.loadtxt(file, delimiter=" ")
    except Exception as e:
        print(e)
        raise "ERROR: " + e
