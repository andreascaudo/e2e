import numpy as np


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
        line_density: float,            # [l/um]
        blaze_angle: float,             # [deg]
        eps_angle: float,               # [deg]

        # Camera f length
        f_camera: float,                # [mm]

        n_pixels: int,                  # [-]
        dimension_pixel: float,         # [um]

        wavelength_min: float,          # [A]
        wavelength_max: float,          # [A]

        resolving_power: float          # [-]

    ) -> None:
        self.name = name
        self.arm = arm
        self.type = type

        self.optical_wavelength_file = optical_wavelength_file
        self.optical_rtc_psf_map_file = optical_rtc_psf_map_file
        self.orders_table_file = orders_table_file
        try:
            self.order_table = np.loadtxt(
                self.orders_table_file, delimiter=" ", skiprows=1)
        except Exception as e:
            print(e)

        self.fwhm_instrument = fwhm_instrument

        self.cp_f_in = common_path_f_in
        self.cp_f_out = common_path_f_out
        self.n_slice = n_slice

        self.slit_size_x = slit_size_x
        self.slit_size_y = slit_size_y
        self.slit_size_x_calibration = slit_size_x_calibration
        self.slit_size_y_calibration = slit_size_y_calibration

        self.f_collimator = f_collimator

        if n_order_start > n_order_end:
            self.n_orders = [*range(n_order_start, n_order_end-1, -1)]
        elif n_order_start < n_order_end:
            self.n_orders = [*range(n_order_start, n_order_end+1, 1)]
        else:
            self.n_orders = [n_order_start]

        self.line_density = line_density
        self.blaze_angle = blaze_angle
        self.eps_angle = eps_angle
        self.f_camera = f_camera

        self.n_pixels = n_pixels
        self.dimension_pixel = dimension_pixel

        self.wavelength_min = wavelength_min
        self.wavelength_max = wavelength_max
        self.wavelength_band = np.arange(wavelength_min, wavelength_max, 1)

        self.resolving_power = resolving_power
