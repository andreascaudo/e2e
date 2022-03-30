import math
from scipy import interpolate
import numpy as np

from . import unit_converter


def gaussian_distribution(x, y, counts, refx, refy, sx, sy):
    peak = counts / (2*sx*sy*math.pi)
    return peak * np.exp(-((x - refx)**2/(2*(sx**2)) + (y - refy)**2/(2*(sy**2))))


def interp(x, y, new_x, fill_value, kind="linear"):
    function = interpolate.interp1d(x, y, kind=kind, fill_value=fill_value)
    return function(new_x)


def integration(lam, delta_lam, flux):
    l1 = lam  # [m]
    l2 = lam + delta_lam  # [m]
    dl = delta_lam / 100  # [m]

    lam_vect = np.arange(l1, l2, dl)

    # print(len(lam_vect))

    spec_flux_norm_interp = interp(unit_converter.wavelength(
        flux.wavelength, "A", "m"), flux.flux, lam_vect, fill_value="extrapolate")  # flux in [phot/s/cm^2/ang] and lambda in [m]

    # E=N*h*nu -> N=E*lambda/h*c
    dl = unit_converter.wavelength(dl, "m", "A")
    counts = 0

    for i in range(0, len(spec_flux_norm_interp)):
        counts = counts + \
            spec_flux_norm_interp[i] * dl  # Nfot/(s*cm^2)

    return counts


def object_slit(counts, telescope, spectograph, acquisition, parameter, image_size, ps_y_fact):
    pixelsz = spectograph.dimension_pixel / parameter.pixel_oversampling

    refx = -((image_size[0]/2) - np.round(image_size[0]/2))
    refy = ((image_size[1]/2) - np.round(image_size[1]/2))

    f = telescope.f_number * telescope.pupil_equivalent_diameter * 10
    plate_scale = (f/206265)*1000
    FWHM_x = plate_scale * acquisition.sky.seeing  # in um
    FWHM_y = plate_scale * acquisition.sky.seeing  # in um

    sx = (FWHM_x/(2.35*pixelsz))*ps_y_fact
    sy = (FWHM_y/(2.35*pixelsz))

    detector = np.zeros((image_size[0], image_size[1]))
    detector_centre = [np.round(image_size[0]/2),
                       np.round(image_size[1]/2)]  # 180, 60

    # image_size is equal to:

    # x = column
    # y = row

    for i in range(0, image_size[1]):
        for j in range(0, image_size[0]):

            displ_x = detector_centre[0] - j - 1
            displ_y = detector_centre[1] - i - 1

            if (displ_x < 0):
                xmin = 0.5 + (abs(displ_x)-1)
                xmax = 0.5 + abs(displ_x)
            elif (displ_x > 0):
                xmin = -0.5 + (displ_x)
                xmax = -0.5 + (displ_x-1)
            elif (displ_x == 0):
                xmin = -0.5
                xmax = 0.5

            if (displ_y < 0):
                ymin = -0.5 - (abs(displ_y))
                ymax = -0.5 - (abs(displ_y)-1)
            elif (displ_y > 0):
                ymin = 0.5 + (displ_y-1)
                ymax = 0.5 + (displ_y)
            elif (displ_y == 0):
                ymin = -0.5
                ymax = 0.5

            # Integrale di volume
            phi_vect = np.zeros(4)
            phi_vect[0] = gaussian_distribution(
                xmin, ymin, counts, refx, refy, sx, sy)
            phi_vect[1] = gaussian_distribution(
                xmax, ymin, counts, refx, refy, sx, sy)
            phi_vect[2] = gaussian_distribution(
                xmax, ymax, counts, refx, refy, sx, sy)
            phi_vect[3] = gaussian_distribution(
                xmin, ymax, counts, refx, refy, sx, sy)

            max_val = max(phi_vect)
            max_index = np.argmax(phi_vect)
            min_val = min(phi_vect)
            min_index = np.argmin(phi_vect)

            # trovo restanti indici
            val_rest = np.zeros(2)
            k_search_rest = 0

            for s_rest in range(0, 4):
                if (s_rest != max_index and s_rest != min_index):
                    val_rest[k_search_rest] = phi_vect[s_rest]
                    k_search_rest = + 1

            vol1 = max_val/2

            vol_b = (((max_val - min_val) + (max_val - val_rest[1]))/2)/3
            vol_d = (((max_val - min_val) + (max_val - val_rest[0]))/2)/3

            detector[j][i] = (vol1 - vol_b) + (vol1 - vol_d)

    return detector


def mask_ideal_slit(image_size, sy_m, sx_m):

    x = np.arange(0, image_size[1], 1)
    y = np.arange(0, image_size[0], 1)

    x, y = np.meshgrid(x, y)

    mask = (x > ((image_size[1]-sx_m)/2)) & (x <= ((image_size[1]+sx_m) /
                                                   2)) & (y > ((image_size[0]-sy_m)/2)) & (y <= ((image_size[0]+sy_m)/2))

    return mask


def rebin_image(image, factor):
    p = int(factor[0])
    print(p)
    q = int(factor[1])
    print(q)
    size_img = image.shape
    print(size_img)

    new_image = np.sum(np.reshape(image, (p, -1), order="F"), axis=0)

    new_image = np.reshape(new_image, (int(size_img[0]/p), -1), order="F").T

    new_image = np.sum(np.reshape(new_image, (q, -1), order="F"), axis=0)
    new_image = np.reshape(new_image, (int(size_img[1]/q), -1), order="F").T
    # print(new_image.shape)

    return new_image
