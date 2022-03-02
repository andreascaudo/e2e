class Sky:
    def __init__(
        self,
        airmass: str,
        moon_fli: float,
        pwv: float,
        seeing: str
    ) -> None:
        self.airmass = airmass
        self.moon_fli = moon_fli
        self.pwv = pwv
        self.seeing = seeing
        # Load spectrum
