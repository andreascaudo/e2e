from distutils.fancy_getopt import OptionDummy

from ..spectrograph import Spectrograph
from ..acquisition import Acquisition
from .output import Output


class Configuration:
    def __init__(
        self,
        output: Output,
        acquisition: Acquisition,
        spectrograph: Spectrograph,
        # telescope
        parameters: dict
    ):
        self.output = output,
        self.acquisition = acquisition
        self.spectrograph = spectrograph
        self.parameters = parameters
