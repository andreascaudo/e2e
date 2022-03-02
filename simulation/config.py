from .output import Output


class Configuration:
    output: Output
    #instrument: Instrument
    #acquisition: Acquisition

    def __init__(self, output):
        self.output = output

    def get_configuration(self):
        return self.output
