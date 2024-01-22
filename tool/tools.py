from astropy.io import fits
import matplotlib.pyplot as plt
import math
from charset_normalizer import detect
from scipy import ndimage
from scipy import interpolate
from scipy import signal
import numpy as np
from numba import njit

from . import unit_converter


@njit()
def gaussian_distribution(x, y, counts, refx, refy, sx, sy):
    peak = counts / (2*sx*sy*math.pi)
    return peak * np.exp(-((x - refx)**2/(2*(sx**2)) + (y - refy)**2/(2*(sy**2))))

# Round floats down to keep one non-zero decimal only


@njit()
def myround(n):
    if n == 0:
        return 0
    sgn = -1 if n < 0 else 1
    scale = int(-math.floor(math.log10(abs(n))))
    if scale <= 0:
        scale = 1
    factor = 10**scale
    return sgn*math.floor(abs(n)*factor)/factor


def interp(x, y, new_x, fill_value, kind="linear"):

    function = interpolate.interp1d(x, y, kind=kind, fill_value=fill_value)

    return function(new_x)

# Return Plate Scale in um/arcsec
# ----------------------> Pupil in cm


def get_plate_scale(f_number, pupil_equivalent_diameter):
    # focal length [mm] = f/N * Pupil diameter [cm -> mm]
    f = f_number * pupil_equivalent_diameter * 10
    plate_scale = (f/206265)*1000  # um/arcsec
    return plate_scale


def integration(lam, delta_lam, flux):
    l1 = lam  # [m]
    l2 = lam + delta_lam  # [m]
    dl = delta_lam / 100  # [m]

    lam_vect = np.arange(l1, l2, dl)

    # get lam in [A]
    lam_temp = unit_converter.wavelength(lam_vect, "m", "A")
    # Find the index of the closest wavelength in the spectrum
    index = np.searchsorted(flux.wavelength, lam_temp)
    # If there is no flux, return 0
    if not flux.flux[index].any():
        return 0

    spec_flux_norm_interp = interp(unit_converter.wavelength(
        flux.wavelength, "A", "m"), flux.flux, lam_vect, fill_value="extrapolate")  # flux in [phot/s/cm^2/ang] and lambda in [m]

    # If there is no flux, return 0 [Just to be sure]
    if not spec_flux_norm_interp.any():
        return 0

    # E=N*h*nu -> N=E*lambda/h*c
    dl = unit_converter.wavelength(dl, "m", "A")
    counts = 0

    for i in range(0, len(spec_flux_norm_interp)):
        counts = counts + \
            spec_flux_norm_interp[i] * dl  # Nfot/(s*cm^2)

    return counts


@njit()
def calibration_slit(counts, image_size, area, ps_y_fact, slit_x, pixel_oversampling):
    d_arcsec = (slit_x/(slit_x*4))/pixel_oversampling
    h_mask_mu = np.ceil(slit_x/d_arcsec)

    detector = np.ones((image_size[0], image_size[1])) * (counts * area) / \
        ((math.pi/4) * ps_y_fact * h_mask_mu * h_mask_mu)
    return detector


@njit()
def object_slit(counts, pixel_platescale, dimension_pixel, seeing, pixel_oversampling, image_size, ps_y_fact):
    pixelsz = dimension_pixel / pixel_oversampling

    refx = -((image_size[0]/2) - np.round(image_size[0]/2))
    refy = ((image_size[1]/2) - np.round(image_size[1]/2))

    # focal length [mm] = f/N * Pupil diameter [cm -> mm]

    plate_scale = pixel_platescale * dimension_pixel  # um/arcsec

    # Seeing in arcsec
    FWHM_x = plate_scale * seeing  # in um
    FWHM_y = plate_scale * seeing  # in um

    sx = (FWHM_x/(2.35*pixelsz))
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
            phi_vect = np.empty(4)

            phi_vect[0] = gaussian_distribution(
                xmin, ymin, counts, refx, refy, sx, sy)[0]
            phi_vect[1] = gaussian_distribution(
                xmax, ymin, counts, refx, refy, sx, sy)[0]
            phi_vect[2] = gaussian_distribution(
                xmax, ymax, counts, refx, refy, sx, sy)[0]
            phi_vect[3] = gaussian_distribution(
                xmin, ymax, counts, refx, refy, sx, sy)[0]

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


@njit()
def mask_maker(x, y, image_size, sx_m, sy_m):
    return (x > int((image_size[1]-1-sx_m)/2)) & (x <= int((image_size[1]-1+sx_m) / 2)) & (y > ((image_size[0]-sy_m)/2)) & (y <= ((image_size[0]+sy_m)/2))
    # OR return (x >= int((image_size[1]-1-sx_m)/2)) & (x <= int((image_size[1]-1+sx_m) / 2)) & (y > ((image_size[0]-sy_m)/2)) & (y <= ((image_size[0]+sy_m)/2)) 50


def mask_ideal_slit(image_size, sy_m, sx_m, mask_oversampling):

    image_size = image_size * mask_oversampling

    x = np.arange(0, image_size[1], 1)
    y = np.arange(0, image_size[0], 1)

    x, y = np.meshgrid(x, y)

    mask = mask_maker(x, y, image_size, sx_m, sy_m)

    return mask


def multi_pinhole():
    return


def rebin_image(image, factor):
    p = int(factor[0])
    q = int(factor[1])
    size_img = image.shape

    new_image = np.sum(np.reshape(image, (p, -1), order="F"), axis=0)
    new_image = np.reshape(new_image, (int(size_img[0]/p), -1), order="F").T

    new_image = np.sum(np.reshape(new_image, (q, -1), order="F"), axis=0)
    new_image = np.reshape(new_image, (int(size_img[1]/q), -1), order="F").T

    return new_image


def interpolate_psf_map(psf_map_shape, order_psf_map_cube, wavelength, order_wavelength_subpix):
    psf_map = np.zeros(psf_map_shape)

    # Iter for each grid-point to interpolate in lambda

    for x in range(0, psf_map_shape[0]):
        for y in range(0, psf_map_shape[1]):
            psf_data_xy = order_psf_map_cube[x][y]
            psf_map[x][y] = interp(
                wavelength, psf_data_xy, order_wavelength_subpix, "extrapolate")  # in subpix

    return psf_map


def interpolate_griddata_psf_map(psf_map_j_norm, v1, v2):
    # Indexing is specified to be 'ij' to be consistent with the output of np.meshgrid and interpn
    x2, y2 = np.meshgrid(v2, v2, indexing='ij')
    # Interpolate the PSF map

    # GRidata is SLOW and interp2d is DEPRECATED
    # It took 0.2 seconds in contrast to 0.0005 seconds with interpn
    # The difference in accuracy is negligible, the maximum difference is 0.00010 while the average is 2.3e-8
    # The difference in speed is huge, interp2d is 400 times faster
    # So I decided to use interpn instead of griddata

    # Interpolate with interpn
    psf_interp = interpolate.interpn(
        (v1, v1), psf_map_j_norm, (x2, y2), method="linear", bounds_error=False, fill_value=None)
    psf_interp = psf_interp/((v1[1]-v1[0])**2)

    return psf_interp, np.sum(psf_interp)


def convolve(image, kernel):
    # return ndimage.convolve(image, kernel)
    return signal.fftconvolve(image, kernel, mode='same')


@njit(cache=True)
def convolve2D(image, kernel, strides=1):
    # Cross Correlation
    kernel = np.flipud(np.fliplr(kernel))

    # Gather Shapes of Kernel + Image + Padding
    xKernShape = kernel.shape[0]
    yKernShape = kernel.shape[1]
    xImgShape = image.shape[0]
    yImgShape = image.shape[1]

    # Shape of Output Convolution
    xOutput = int(((xImgShape - xKernShape) / strides) + 1)
    yOutput = int(((yImgShape - yKernShape) / strides) + 1)
    output = np.zeros((xOutput, yOutput))

    # Iterate through image
    for y in range(image.shape[1]):
        # Exit Convolution
        if y > image.shape[1] - yKernShape:
            break
        # Only Convolve if y has gone down by the specified Strides
        if y % strides == 0:
            for x in range(image.shape[0]):
                # Go to next row once kernel is out of bounds
                if x > image.shape[0] - xKernShape:
                    break
                try:
                    # Only Convolve if x has moved by the specified Strides
                    if x % strides == 0:
                        output[x, y] = (
                            kernel * image[x: x + xKernShape, y: y + yKernShape]).sum()
                except:
                    break

    return output


def flip_slit(detector, to_tilt):
    tilt = -to_tilt
    # Rotate the detector by tilt
    return ndimage.rotate(detector, tilt, reshape=False)


@njit()
def normalize_psf_map(psf_map_j):
    psf_map_j_norm = psf_map_j / np.sum(psf_map_j)
    psf_box_z = psf_map_j_norm.shape[0]

    return psf_map_j_norm, psf_box_z


@njit()
def init_conv_matrix(psf_map_pixel_number, dimension_pixel, psf_box_z):
    v1 = np.linspace(0, (psf_map_pixel_number *
                         dimension_pixel), psf_box_z)
    v2 = np.linspace(0, (psf_map_pixel_number * dimension_pixel),
                     (psf_map_pixel_number * dimension_pixel))
    return v1, v2

# CUBES


def add_pre_over_scan(detector):
    # Adding the 16 Pixels to get 9232 in spatial direction
    detecor_final = np.zeros(
        (detector.shape[0]+16, detector.shape[1]))
    detecor_final[8:-8, :] = detector

    detector = detecor_final

    # Rotate detector by 90 degrees
    detector = np.rot90(detector)

    new_detector = np.zeros((9920, 9296))
    new_detector[24:1176, 0:4616] = detector[0:1152, 0:4616]
    new_detector[24:1176, 4680:-1] = detector[0:1152, 4616:-1]
    new_detector[1264:2416, 0:4616] = detector[1152:(2*1152), 0:4616]
    new_detector[1264:2416, 4680:-1] = detector[1152:(2*1152), 4616:-1]
    new_detector[2504:3656, 0:4616] = detector[(2*1152):(3*1152), 0:4616]
    new_detector[2504:3656, 4680:-1] = detector[(2*1152):(3*1152), 4616:-1]
    new_detector[3744:4896, 0:4616] = detector[(3*1152):(4*1152), 0:4616]
    new_detector[3744:4896, 4680:-1] = detector[(3*1152):(4*1152), 4616:-1]

    new_detector[5024:(5024+1152), 0:4616] = detector[(4*1152):(5*1152), 0:4616]
    new_detector[5024:(5024+1152),
                 4680:-1] = detector[(4*1152):(5*1152), 4616:-1]
    new_detector[6264:(6264+1152), 0:4616] = detector[(5*1152):(6*1152), 0:4616]
    new_detector[6264:(6264+1152),
                 4680:-1] = detector[(5*1152):(6*1152), 4616:-1]
    new_detector[7504:(7504+1152), 0:4616] = detector[(6*1152):(7*1152), 0:4616]
    new_detector[7504:(7504+1152),
                 4680:-1] = detector[(6*1152):(7*1152), 4616:-1]
    new_detector[8744:(8744+1152), 0:4616] = detector[(7*1152):(8*1152), 0:4616]
    new_detector[8744:(8744+1152), 4680:-
                 1] = detector[(7*1152):(8*1152), 4616:-1]

    # Rotate back new_detector to the original orientation
    new_detector = np.rot90(new_detector, 3)
    return new_detector

# CUBES


def remove_pre_over_scan(detector):
    # Rotate detector by 90 degrees
    detector = np.rot90(detector)

    new_detector = np.zeros((9216, 9216+16))

    new_detector[0:1152, 0:4616] = detector[24:1176, 0:4616]
    new_detector[0:1152, 4616:-1] = detector[24:1176, 4680:-1]
    new_detector[1152:(2*1152), 0:4616] = detector[1264:2416, 0:4616]
    new_detector[1152:(2*1152), 4616:-1] = detector[1264:2416, 4680:-1]
    new_detector[(2*1152):(3*1152), 0:4616] = detector[2504:3656, 0:4616]
    new_detector[(2*1152):(3*1152), 4616:-1] = detector[2504:3656, 4680:-1]
    new_detector[(3*1152):(4*1152), 0:4616] = detector[3744:4896, 0:4616]
    new_detector[(3*1152):(4*1152), 4616:-1] = detector[3744:4896, 4680:-1]

    new_detector[(4*1152):(5*1152),
                 0:4616] = detector[5024:(5024+1152), 0:4616]
    new_detector[(4*1152):(5*1152), 4616:-
                 1] = detector[5024:(5024+1152), 4680:-1]
    new_detector[(5*1152):(6*1152),
                 0:4616] = detector[6264:(6264+1152), 0:4616]
    new_detector[(5*1152):(6*1152), 4616:-
                 1] = detector[6264:(6264+1152), 4680:-1]
    new_detector[(6*1152):(7*1152),
                 0:4616] = detector[7504:(7504+1152), 0:4616]
    new_detector[(6*1152):(7*1152), 4616:-
                 1] = detector[7504:(7504+1152), 4680:-1]
    new_detector[(7*1152):(8*1152),
                 0:4616] = detector[8744:(8744+1152), 0:4616]
    new_detector[(7*1152):(8*1152), 4616:-
                 1] = detector[8744:(8744+1152), 4680:-1]

    # Rotate back new_detector to the original orientation
    new_detector = np.rot90(new_detector, 3)

    new_data = np.zeros(
        (new_detector.shape[0]-16, new_detector.shape[1]))
    new_data = new_detector[8:-8, :]

    # (9216, 9216)
    return new_data
