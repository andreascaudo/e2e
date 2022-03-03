from distutils.fancy_getopt import OptionDummy

from ..instrument import Spectrograph
from ..instrument import Telescope
from ..acquisition import Acquisition
from .output import Output


class Configuration:
    def __init__(
        self,
        output: Output,                     # Output Specification
        acquisition: Acquisition,           # Acquisition Parameters
        telescope: Telescope,               # Telescope Parameters
        spectrograph: Spectrograph,         # Spectrograph Parameters
        parameters: dict                    # Simulation Paramenters
    ):
        self.output = output,
        self.acquisition = acquisition
        self.telescope = telescope
        self.spectrograph = spectrograph
        self.parameters = parameters
