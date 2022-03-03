class Characteristics:
    def __init__(

        self,
        slit_size_x: float,                     # [arcsec]
        slit_size_y: float,                     # [arcsec]
        detector_integration_time: float,       # [s]
        number_of_dit: int,                     # [-]
        number_of_integration: int,             # [-]
        bin_x: int,                             # [-]
        bin_y: int                              # [-]

    ) -> None:

        self.slit_size_x = slit_size_x
        self.slit_size_y = slit_size_y
        self.detector_integration_time = detector_integration_time
        self.number_of_dit = number_of_dit
        self.number_of_integration = number_of_integration
        self.bin_x = bin_x
        self.bin_y = bin_y
