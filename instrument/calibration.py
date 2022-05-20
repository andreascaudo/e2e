import math
import numpy as np
from ..tool import tools
from numba import njit


class Calibration:
    def __init__(
        self,
        mode: str
    ) -> None:
        self.mode = mode


class Multipinhole(Calibration):
    def __init__(
        self,
        mode: str,
        pinhole_radius: float,
        pinhole_radius_pixel: float,
        holes_number: int
    ) -> None:
        super().__init__(mode)
        self.pinhole_radius = pinhole_radius
        self.pinhole_radius_pixel = pinhole_radius_pixel
        self.area = math.pi * (pinhole_radius**2)
        self.holes_number = holes_number

    # TO USE JIT THIS SHOULD BE A STATIC METHOD
    def get_mask(self, ps_y_fact, image_size, pixel_oversampling, mask_oversampling):
        radius_sub = self.pinhole_radius_pixel * pixel_oversampling * mask_oversampling
        d_xy = 4 * radius_sub
        # central hole: position and size (r=1 pixel --> 1*PPP*FF sub-sub-pixes)
        central_hole_xy = (d_xy+1)/2  # +1 only with odd numbers

        radius_y = radius_sub * ps_y_fact

        x = np.arange(0, d_xy, 1)
        y = np.arange(0, d_xy, 1)
        x, y = np.meshgrid(x, y)

        # Init MASK
        mask_ph_m = np.zeros(
            [image_size[0]*mask_oversampling, image_size[1]*mask_oversampling])

        # Where submask equals 1
        mask_ph_m1 = (((x-central_hole_xy)**2)/(radius_sub**2)) + \
            (((y-central_hole_xy)**2)/((radius_y**2))) < 1

        # position the sub-mask in the mask
        # PLACING THE MASK ON ALL OTHER POSITIONS --2
        holes_number_2 = np.floor(self.holes_number/2)
        dist_holes = 158  # um

        dist_holes_pix = dist_holes*(4/110)*(36/48)  # to go in pixel scale
        dist_holes_pix1 = dist_holes_pix * pixel_oversampling
        dist_holes_pix2 = dist_holes_pix1 * mask_oversampling
        displ_y = np.round(np.arange(-holes_number_2,
                                     holes_number_2+1, 1) * dist_holes_pix2)

        for i in range(0, self.holes_number):
            y_submask_01 = int((
                (image_size[0]*mask_oversampling)/2) - ((np.shape(mask_ph_m1)[0])/2) - displ_y[i])

            x_submask_01 = int((
                (image_size[1]*mask_oversampling)/2) - ((np.shape(mask_ph_m1)[1])/2))

            y_submask_11 = int((
                (image_size[0]*mask_oversampling)/2) + ((np.shape(mask_ph_m1)[0])/2) - displ_y[i])

            x_submask_11 = int((
                (image_size[1]*mask_oversampling)/2) + ((np.shape(mask_ph_m1)[1])/2))

            mask_ph_m[y_submask_01:y_submask_11,
                      x_submask_01:x_submask_11] = mask_ph_m1

        # rebin to sub-pixel scale and re-normalize to 1
        mask = tools.rebin_image(
            mask_ph_m, [mask_oversampling, mask_oversampling])
        mask = mask/np.max(mask)
        return mask
