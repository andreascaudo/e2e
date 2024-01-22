import yaml

from .output import Output
from ..acquisition import Acquisition
from ..instrument import Spectrograph
from ..instrument import Zemax
from ..instrument import calibration
from ..instrument import Telescope
from .config import Configuration, Parameter


def to_output(dct: dict):
    return Output(**dct)


def to_acquisition(dct: dict, spectrograph_obj):
    return Acquisition(spectrograph_obj, **dct)


def to_calibration(dct: dict):
    calibration_mode = getattr(calibration, dct["mode"])
    return calibration_mode(**dct)


def to_telescope(dct: dict):
    return Telescope(**dct)


def to_spectrograph(dct: dict):
    return Spectrograph(**dct)


def to_zemax(dct: dict):
    return Zemax(**dct)


def to_parameter(dct: dict):
    return Parameter(**dct)


def build_config(configuration_file: dict) -> Configuration:
    # Dictionary -> Object of specific Class
    output_obj = to_output(configuration_file["output"])
    calibration_obj = to_calibration(
        configuration_file["calibration"]) if "calibration" in configuration_file else None
    telescope_obj = to_telescope(configuration_file["telescope"])
    zemax_obj = to_zemax(
        configuration_file["zemax"]) if "zemax" in configuration_file else None
    spectrograph_obj = to_spectrograph(configuration_file["spectrograph"])
    acquisition_obj = to_acquisition(
        configuration_file["acquisition"], spectrograph_obj)
    parameter_obj = to_parameter(configuration_file["simulation"])

    # Dictionary -> Dictionary

    # Setup Environment Configuration
    config = Configuration(
        output=output_obj,
        acquisition=acquisition_obj,
        calibration=calibration_obj,
        telescope=telescope_obj,
        spectrograph=spectrograph_obj,
        zemax=zemax_obj,
        parameters=parameter_obj
    )

    return config


# Load YAML file configuration
def load_config(path: str):
    with open(path, "r") as stream:
        try:
            configuration_file = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            if hasattr(exc, 'problem_mark'):
                mark = exc.problem_mark
                print("Error position: (%s:%s)" % (mark.line+1, mark.column+1))
        config = build_config(configuration_file)

    return config
