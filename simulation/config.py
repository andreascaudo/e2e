from ..instrument import Spectrograph
from ..instrument import Telescope
from ..acquisition import Acquisition
from .output import Output


class Parameter:
    def __init__(
        self,
        # pixel wavelength
        pixel_oversampling: int,
        wavelength_overscanning: float,
        mask_oversampling: float,
        # psf
        psf_pupil_sampling: float,
        psf_field_sampling: int,
        psf_map_pixel_number: int
    ):
        self.pixel_oversampling = pixel_oversampling
        self.wavelength_overscanning = wavelength_overscanning
        self.mask_oversampling = mask_oversampling
        self.psf_pupil_sampling = psf_pupil_sampling
        self.psf_field_sampling = psf_field_sampling
        self.psf_map_pixel_number = psf_map_pixel_number
        self.psf_box_subpix = self.psf_map_pixel_number * self.pixel_oversampling


class Configuration:
    def __init__(
        self,
        output: Output,                     # Output Specification
        acquisition: Acquisition,           # Acquisition Parameters
        telescope: Telescope,               # Telescope Parameters
        spectrograph: Spectrograph,         # Spectrograph Parameters
        parameters: Parameter                    # Simulation Paramenters
    ):
        self.output = output,
        self.acquisition = acquisition
        self.telescope = telescope
        self.spectrograph = spectrograph
        self.parameters = parameters
