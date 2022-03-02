class Sed:
    def __init__(
        self,
        spectrum_file: str,
        temperature: float,
        magnitude: float,
        magnitude_system: str
    ) -> None:
        self.spectrum_file = spectrum_file
        self.temperature = temperature
        self.magnitude = magnitude
        self.magnitude_system = magnitude_system
        # Load spectrum
