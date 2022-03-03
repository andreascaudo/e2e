from .sed import Sed
from .sky import Sky
from .characteristics import Characteristics


class Acquisition:
    def __init__(self,
                 sed: Sed,
                 sky: Sky,
                 characteristics: Characteristics
                 ) -> None:
        self.sed = Sed(**sed)
        self.sky = Sky(**sky)
        self.characteristics = Characteristics(**characteristics)
