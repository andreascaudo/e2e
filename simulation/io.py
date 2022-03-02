import yaml

from .output import Output
from ..acquisition import Acquisition
from ..spectrograph import Spectrograph
from .config import Configuration


def to_output(dct: dict):
    return Output(**dct)


def to_acquisition(dct: dict):
    return Acquisition(**dct)


def to_spectrograph(dct: dict):
    return Spectrograph(**dct)


def build_config(configuration_file: dict):
    output_obj = to_output(configuration_file["output"])
    acquisition_obj = to_acquisition(configuration_file["acquisition"])
    spectrograph_obj = to_spectrograph(configuration_file["spectrograph"])
    parameters_dict = configuration_file["simulation"]

    config = Configuration(
        output=output_obj,
        acquisition=acquisition_obj,
        spectrograph=spectrograph_obj,
        parameters=parameters_dict
    )

    return config


def load_config(path: str):
    with open(path, "r") as stream:
        try:
            configuration_file = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

        config = build_config(configuration_file)

    return config
