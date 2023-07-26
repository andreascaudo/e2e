import numpy as np


class Slice:
    def __init__(
        self,
        slice_id: int,
        shift_arc: float
    ) -> None:
        self.slice_id = slice_id
        self.shift_arc = shift_arc

    def to_pix_arcsec(self, plate_scale):
        self.shift_pix = self.shift_arc * plate_scale

    def to_subpix(self, oversampling):
        self.shift_subpix = self.shift_pix * oversampling
