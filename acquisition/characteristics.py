class Characteristics:
    def __init__(
        self,
        spectrograph_obj,
        slit_size_x_index: int,                     # [arcsec]
        slit_size_y_index: int,                     # [arcsec]
        detector_integration_time: float,       # [s]
        number_of_dit: int,                     # [-]
        number_of_integration: int,             # [-]
        bin_x: int,                             # [-]
        bin_y: int                              # [-]
    ) -> None:
        if slit_size_x_index > len(spectrograph_obj.slit_size_x):
            raise ValueError("Slit size x index out of range")
        if slit_size_y_index > len(spectrograph_obj.slit_size_y):
            raise ValueError("Slit size y index out of range")

        self.slit_size_x = spectrograph_obj.slit_size_x[slit_size_x_index]
        self.slit_size_y = spectrograph_obj.slit_size_y[slit_size_y_index]
        self.slit_size_x_calibration = spectrograph_obj.slit_size_x_calibration[
            slit_size_x_index]
        self.slit_size_y_calibration = spectrograph_obj.slit_size_y_calibration[
            slit_size_y_index]
        self.slit_area_calibration = self.slit_size_x_calibration * \
            self.slit_size_y_calibration

        self.detector_integration_time = detector_integration_time
        self.number_of_dit = number_of_dit
        self.number_of_integration = number_of_integration

        self.bin_x = bin_x
        self.bin_y = bin_y
