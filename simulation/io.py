import yaml

from .output import Output
from ..acquisition import Acquisition
from ..instrument import Spectrograph
from ..instrument import Telescope
from .config import Configuration, Parameter


def to_output(dct: dict):
    return Output(**dct)


def to_acquisition(dct: dict):
    return Acquisition(**dct)


def to_telescope(dct: dict):
    return Telescope(**dct)


def to_spectrograph(dct: dict):
    return Spectrograph(**dct)


def to_parameter(dct: dict):
    return Parameter(**dct)


def build_config(configuration_file: dict) -> Configuration:
    # Dictionary -> Object of specific Class
    output_obj = to_output(configuration_file["output"])
    acquisition_obj = to_acquisition(configuration_file["acquisition"])
    telescope_obj = to_telescope(configuration_file["telescope"])
    spectrograph_obj = to_spectrograph(configuration_file["spectrograph"])
    parameter_obj = to_parameter(configuration_file["simulation"])

    # Dictionary -> Dictionary

    # Setup Environment Configuration
    config = Configuration(
        output=output_obj,
        acquisition=acquisition_obj,
        telescope=telescope_obj,
        spectrograph=spectrograph_obj,
        parameters=parameter_obj
    )

    return config


# Load YAML file configuration
def load_config(path: str):
    with open(path, "r") as stream:
        try:
            configuration_file = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

        config = build_config(configuration_file)

    return config
