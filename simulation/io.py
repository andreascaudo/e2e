# This File should be moved outside the folder simulation and into "io" subpackages where output object is specified
# (or not really it could be the the IO of simulation)

import typing as t
import yaml

from .output import Output
from .config import Configuration

#from .instrument import spectograph
#from .acquisition import sed


def to_output(dct: dict):
    return Output(**dct)


'''
def to_instrument():
    # TBW
    # create class containing instrument specification
    # Class instrument - SubClass Spectograph


def to_acquisition():
    # TBW
    # create class containing sed & sky specification
'''


def build_config(configuration_file: dict):
    output_obj = to_output(configuration_file["output"])

    config = Configuration(
        output=output_obj
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
