from . import sed as sd
from .sky import Sky
from .characteristics import Characteristics


class Acquisition:
    def __init__(self,
                 sed: sd.Sed,
                 characteristics: Characteristics,
                 sky: Sky = None
                 ) -> None:
        sed_class = getattr(sd, sed["sed_type"])
        self.sed = sed_class(**sed)
        if sky is not None:
            self.sky = Sky(**sky)
        else:
            self.sky = Sky(-1, -1, -1, -1, active=False)
        self.characteristics = Characteristics(**characteristics)
