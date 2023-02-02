
class Telescope:
    def __init__(
        self,
        name: str,                          # Telescope Name
        unit: int,                          # Unit Telescope
        diameter: float,                    # [cm]
        l_zero: float,                      # [m]
        pupil_equivalent_diameter: float,   # [cm]
        f_number: float,                    # focal_length/diameter
    ) -> None:
        self.name = name

        self.unit = unit
        self.diameter = diameter
        self.l_zero = l_zero
        self.pupil_equivalent_diameter = pupil_equivalent_diameter

        self.f_number = f_number
