from .simulation import Configuration


def run(configuration: Configuration):
    # Unpack configuration
    output = configuration.output
    acquisition = configuration.acquisition
    telescope = configuration.telescope
    spectrograph = configuration.spectrograph

    # First step: generate flux
    flux = acquisition.sed.get_flux()
