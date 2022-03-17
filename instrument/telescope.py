
class Telescope:
    def __init__(
        self,
        name: str,                          # Telescope Name
        diameter: float,                    # [cm]
        l_zero: float,                      # [m]
        pupil_equivalent_diameter: float,   # [cm]
        f_number: float,                    # focal_length/diameter
        image_scale: float                  # [arcsec/mm]
    ) -> None:
        self.name = name

        self.diameter = diameter
        self.l_zero = l_zero
        self.pupil_equivalent_diameter = pupil_equivalent_diameter

        self.f_number = f_number

        self.image_scale = image_scale
