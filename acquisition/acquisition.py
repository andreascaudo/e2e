from . import sed as sd
from .sky import Sky
from .characteristics import Characteristics


class Acquisition:
    def __init__(self,
                 sed: sd.Sed,
                 sky: Sky,
                 characteristics: Characteristics
                 ) -> None:
        sed_class = getattr(sd, sed["sed_type"])
        self.sed = sed_class(**sed)
        self.sky = Sky(**sky)
        self.characteristics = Characteristics(**characteristics)
