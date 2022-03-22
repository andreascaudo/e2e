import numpy as np
from scipy import interpolate
import math

import matplotlib.pyplot as plt


class Spectrograph:
    def __init__(
        self,
        name: str,                      # Spectrograph Name
        arm: str,                       # Spectrograph Arm
        type: str,                      # Spectrograph Type

        optical_wavelength_file: str,
        optical_rtc_psf_map_file: str,
        orders_table_file: str,

        # Full Width Half Maximum
        fwhm_instrument: float,

        # Common Path f numbers (Pre-Slit unit)

        common_path_f_in: float,        # [-]
        common_path_f_out: float,       # [-]

        # Common Path anamorphic factor (Pre-Slit unit)

        common_path_AFxy: float,

        # Slicing

        n_slice: int,                   # [-]

        slit_size_x: list,              # [arcsec]
        slit_size_y: list,              # [arcsec]
        slit_size_x_calibration: list,  # cm
        slit_size_y_calibration: list,  # cm

        f_collimator: float,            # [-]

        n_order_start: int,             # [-]
        n_order_end: int,               # [-]
        line_density: float,            # [l/um] #rho
        blaze_angle: float,             # [deg]
        eps_angle: float,               # [deg]

        # Camera f length
        f_camera: float,                # [mm]

        n_pixels: int,                  # [-]
        dimension_pixel: float,         # [um]

        wavelength_min: float,          # [A]
        wavelength_max: float,          # [A]

        resolving_power: float,          # [-]

        # Efficiency
        ilg: float,
        n_p: float,
        grating_lambda: list,
        grating_efficiency: list,
        telescope_fdr_file: str,
        commonpath_ir_fdr_file: str,
        commonpath_vis_fdr_file: str,
        collimator_fdr_file: str,
        field_mirror_fdr_file: str,
        cross_disperser_fdr_file: str,
        fold_mirror_fdr_file: str,
        camera_fdr_file: str,
        qe_detector_file: str


    ) -> None:
        self.name = name
        self.arm = arm
        self.type = type

        self.optical_wavelength_file = optical_wavelength_file
        self.optical_rtc_psf_map_file = optical_rtc_psf_map_file
        self.orders_table_file = orders_table_file

        if n_order_start > n_order_end:
            self.n_orders = [*range(n_order_start, n_order_end-1, -1)]
        elif n_order_start < n_order_end:
            self.n_orders = [*range(n_order_start, n_order_end+1, 1)]
        else:
            self.n_orders = [n_order_start]

        self.len_n_orders = len(self.n_orders)

        try:
            self.order_table = np.loadtxt(
                self.orders_table_file, delimiter=" ", skiprows=1)

            self.sx_wavelegnth_per_order = []
            self.dx_wavelegnth_per_order = []

            # Create from order table two list containing the wavelength min and max for each order
            for n_ord in self.n_orders:
                self.sx_wavelegnth_per_order.append(
                    self.order_table[self.order_table[:, 0] == n_ord, :][0][2])
                self.dx_wavelegnth_per_order.append(
                    self.order_table[self.order_table[:, 0] == n_ord, :][-1][2])

        except Exception as e:
            print(e)

        self.telescope_fdr_file = load_fdr(telescope_fdr_file)
        self.commonpath_ir_fdr_file = load_fdr(commonpath_ir_fdr_file)
        self.commonpath_vis_fdr_file = load_fdr(commonpath_vis_fdr_file)
        self.collimator_fdr_file = load_fdr(collimator_fdr_file)
        self.field_mirror_fdr_file = load_fdr(field_mirror_fdr_file)
        self.cross_disperser_fdr_file = load_fdr(cross_disperser_fdr_file)
        self.fold_mirror_fdr_file = load_fdr(fold_mirror_fdr_file)
        self.camera_fdr_file = load_fdr(camera_fdr_file)
        self.qe_detector_file = load_fdr(qe_detector_file)

        self.fwhm_instrument = fwhm_instrument

        self.cp_f_in = common_path_f_in
        self.cp_f_out = common_path_f_out
        self.n_slice = n_slice

        self.slit_size_x = slit_size_x
        self.slit_size_y = slit_size_y
        self.slit_size_x_calibration = slit_size_x_calibration
        self.slit_size_y_calibration = slit_size_y_calibration

        self.f_collimator = f_collimator

        self.line_density = line_density
        self.blaze_angle = blaze_angle
        self.eps_angle = eps_angle
        self.f_camera = f_camera

        self.ilg = ilg
        self.n_p = n_p * 3

        self.grating_lambda = grating_lambda
        self.grating_efficiency = grating_efficiency

        if self.type == "echelle":
            self.wavematrix, self.b_phase, self.b = echelle_grating_efficiency(
                self.line_density, self.blaze_angle, self.eps_angle, self.sx_wavelegnth_per_order, self.dx_wavelegnth_per_order, self.n_orders, self.len_n_orders, self.ilg, self.n_p)

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
                self.grating_fdr[i] = interp(
                    self.grating_lambda, self.grating_efficiency, self.wavematrix[i], fill_value="extrapolate")
                self.grating_fdr[i] = self.grating_fdr[i] * self.b[i]

                self.telescope_fdr[i] = interp(
                    self.telescope_fdr_file.T[0], self.telescope_fdr_file.T[1], self.wavematrix[i], fill_value="extrapolate")
                self.commonpath_ir_fdr[i] = interp(
                    CPIR_fdr.T[0], CPIR_fdr.T[1], self.wavematrix[i], fill_value="extrapolate", kind="nearest")
                self.collimator_fdr[i] = interp(
                    self.collimator_fdr_file.T[0], self.collimator_fdr_file.T[1], self.wavematrix[i], fill_value="extrapolate", kind="nearest")
                self.field_mirror_fdr[i] = interp(
                    self.field_mirror_fdr_file.T[0], self.field_mirror_fdr_file.T[1], self.wavematrix[i], fill_value="extrapolate", kind="nearest")
                self.cross_disperser_fdr[i] = interp(
                    self.cross_disperser_fdr_file.T[0], self.cross_disperser_fdr_file.T[1], self.wavematrix[i], fill_value="extrapolate", kind="nearest")
                self.fold_mirror_fdr[i] = interp(
                    self.fold_mirror_fdr_file.T[0], self.fold_mirror_fdr_file.T[1], self.wavematrix[i], fill_value="extrapolate", kind="nearest")
                self.camera_fdr[i] = interp(self.camera_fdr_file.T[0], self.camera_fdr_file.T[1],
                                            self.wavematrix[i], fill_value="extrapolate", kind="nearest")
                self.qe_detector[i] = interp(self.qe_detector_file.T[0], self.qe_detector_file.T[1],
                                             self.wavematrix[i], fill_value="extrapolate", kind="nearest")

        for len_i in range(0, self.len_n_orders):
            plt.plot(self.wavematrix[len_i], self.b[len_i])

        plt.show()

        for len_i in range(0, self.len_n_orders):
            plt.plot(self.wavematrix[len_i], self.grating_fdr[len_i])

        plt.show()

        self.telescope_efficiency_total_fdr = self.telescope_fdr * self.commonpath_ir_fdr * self.collimator_fdr * \
            self.field_mirror_fdr * self.grating_fdr * self.cross_disperser_fdr * \
            self.fold_mirror_fdr * self.camera_fdr * self.qe_detector

        self.spectograph_efficiency_total_fdr = self.commonpath_ir_fdr * self.collimator_fdr * self.field_mirror_fdr * \
            self.grating_fdr * self.cross_disperser_fdr * \
            self.fold_mirror_fdr * self.camera_fdr * self.qe_detector

        for len_i in range(0, self.len_n_orders):
            plt.plot(self.wavematrix[len_i],
                     self.commonpath_ir_fdr[len_i])

        plt.show()

        for len_i in range(0, self.len_n_orders):
            plt.plot(self.wavematrix[len_i],
                     self.telescope_efficiency_total_fdr[len_i])

        plt.show()

        for len_i in range(0, self.len_n_orders):
            plt.plot(self.wavematrix[len_i],
                     self.spectograph_efficiency_total_fdr[len_i])

        plt.show()

        self.n_pixels = n_pixels
        self.dimension_pixel = dimension_pixel

        self.wavelength_min = wavelength_min
        self.wavelength_max = wavelength_max
        self.wavelength_band = np.arange(wavelength_min, wavelength_max, 1)

        self.resolving_power = resolving_power


def cut_spurius_efficiency(b, len_n_orders):
    for i in range(0, len_n_orders):
        max_idx = np.where(b[i] == np.amax(b[i]))[0]
        for k in range(0, max_idx[0]):
            if b[i][k] > b[i][k+1]:
                b[i][k] = 0
            else:
                break
        for k in reversed(range(max_idx[0], len(b[i]))):
            if b[i][k] > b[i][k-1]:
                b[i][k] = 0
            else:
                break
    return b


def echelle_grating_efficiency(line_density, blaze_angle, eps_angle, sx_wavelegnth_per_order, dx_wavelegnth_per_order, n_orders, len_n_orders, ilg, n_p):
    wave_matrix = np.zeros((len_n_orders, n_p))
    b_phase = np.zeros((len_n_orders, n_p))
    b = np.zeros((len_n_orders, n_p))

    # d=grating const and N=number of lines
    d = 1/line_density  # [um]
    N = round(ilg*line_density*1000)  # [-]

    for i in range(0, len_n_orders):
        wave_matrix[i] = np.linspace(
            sx_wavelegnth_per_order[i], dx_wavelegnth_per_order[i], n_p)
        for index, lmbd in enumerate(wave_matrix[i]):
            beta = math.asin(((n_orders[i] * line_density * lmbd) / math.cos(
                math.radians(eps_angle))) - math.sin(math.radians(blaze_angle)))

            b_phase[i, index] = (math.pi*d*math.cos(math.radians(blaze_angle))/lmbd) * (math.sin(
                math.radians(blaze_angle-blaze_angle)) + math.sin(beta-math.radians(blaze_angle)))

            b[i, index] = np.sinc(b_phase[i][index]/math.pi)**2

    b = cut_spurius_efficiency(b, len_n_orders)

    return wave_matrix, b_phase, b


def load_fdr(file):
    try:
        return np.loadtxt(file, delimiter=" ")
    except Exception as e:
        print(e)
        raise "ERROR: " + e


def interp(x, y, new_x, fill_value, kind="linear"):
    function = interpolate.interp1d(x, y, kind=kind, fill_value=fill_value)
    return function(new_x)
