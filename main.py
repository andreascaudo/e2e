from .simulation import Configuration
from .tool import tools
from .tool import efficiency
from .tool import unit_converter
from .tool import generic as plt
import numpy as np
import math
import time

DEBUG = False


def run(configuration: Configuration):
    # Unpack configuration
    telescope = configuration.telescope
    spectrograph = configuration.spectrograph
    acquisition = configuration.acquisition
    parameter = configuration.parameters

    # TBI: Implement a function TO CHECK if acquisition reflects the spectrograph parameters

    # 1st step: generate flux
    # [Angstrom], [Ph/s/cm^2/A]
    sed_wavelength, sed_flux = acquisition.sed.get_flux()

    if DEBUG:
        plt("SED", sed_wavelength, sed_flux, [
            "wavelength [$\AA$]", "flux " + r"[Ph/s/cm$^2$/$\AA$]"])

    sed_wavelength, sed_flux = acquisition.sed.normalize()

    if DEBUG:
        plt("SED Normalized", sed_wavelength, sed_flux, [
            "wavelength [$\AA$]", "flux " + r"[Ph/s/cm$^2$/$\AA$]"])

    # 2nd step: get sky radiance and transission
    # [-] ,[ph/s/cm2/A]
    transimission, radiance = acquisition.sky.get_sky(
        acquisition.characteristics.slit_size_x, acquisition.characteristics.slit_size_y)

    if DEBUG:
        plt("Sky transmission", transimission.wavelength,
            transimission.transmission, ["wavelength [$\AA$]", "[-]"])

        plt("Sky radiance", radiance.wavelength,
            radiance.flux, ["wavelength [$\AA$]", "[ph/s/cm2/A]"])

    # 3th step: Set Sky and Obj Efficiency

    acquisition.sky.set_efficiency(spectrograph.wavematrix,
                                   spectrograph.telescope_spectrograph_efficiency_fdr)
    acquisition.sed.set_efficiency(acquisition.sky.transmission.transmission_matrix,
                                   spectrograph.wavematrix, spectrograph.telescope_spectrograph_efficiency_fdr)

    # 4th step: Get Slit Efficiency and Image Quality
    slit_efficiency_matrix, fwhm_iq_matrix = efficiency.get_slit_efficiency(spectrograph.wavematrix, acquisition.sky.airmass,
                                                                            acquisition.characteristics.slit_size_x, acquisition.characteristics.slit_size_y,
                                                                            acquisition.sky.seeing, spectrograph.fwhm_instrument, (telescope.diameter/100), telescope.l_zero)

    spectrograph.set_subpixels(
        parameter.pixel_oversampling, parameter.psf_map_pixel_number)

    print("Start Calculations")
    start_time = time.time()

    calculation(configuration)

    print("End Calculations")
    print("--- %s seconds ---" % (time.time() - start_time))


def calculation(configuration):
    # Unpack configuration
    telescope = configuration.telescope
    spectrograph = configuration.spectrograph
    acquisition = configuration.acquisition
    parameter = configuration.parameters

    number_of_lambda = np.zeros((spectrograph.len_n_orders))
    ps_y_fact = np.zeros((spectrograph.len_n_orders))
    order_y_subpix_min = np.zeros((spectrograph.len_n_orders))
    psf_bin_mat = np.zeros(
        (spectrograph.n_pixels_subpixel, spectrograph.n_pixels_subpixel, spectrograph.len_n_orders))

    print("Start Orders")
    order_time = time.time()

    for i in [4]:  # range(0, spectrograph.len_n_orders)
        # ---------------------------------------------------------------------
        # order and single-order-table
        order_number = spectrograph.n_orders[i]
        order = spectrograph.order_table[spectrograph.order_table.T[0]
                                         == order_number]
        # wavlength vector in meters
        wavelength = unit_converter.wavelength(
            order.T[2], "um", "m")  # in meters
        delta_lambda = abs(wavelength[0]-wavelength[1])
        # ---------------------------------------------------------------------

        # ---------------------------------------------------------------------
        # ### Sub-Pixel scale: x, lambda, y, Sx
        # X
        order_x_subpix = order.T[3] * parameter.pixel_oversampling
        order_x_subpix_start = math.ceil(
            order.T[3][0] * parameter.pixel_oversampling)
        order_x_subpix_end = math.ceil(
            order.T[3][-1] * parameter.pixel_oversampling)
        order_x_subpix_new = np.arange(
            order_x_subpix_start, order_x_subpix_end+1, 1)
        # Lambda
        order_wavelength_subpix = tools.interp(
            order_x_subpix, wavelength, order_x_subpix_new, "extrapolate", "cubic")
        delta_lambda_subpix = np.diff(order_wavelength_subpix)

        # Y
        order_y_subpix = order.T[4] * parameter.pixel_oversampling
        order_y_subpix = tools.interp(
            order_x_subpix, order_y_subpix, order_x_subpix_new, "extrapolate", "cubic")

        # Sx
        order_sx_subpix = order.T[5]
        order_sx_subpix = tools.interp(
            order_x_subpix, order_sx_subpix, order_x_subpix_new, "extrapolate", "cubic")

        # Tilt Order
        order_tilt = order.T[13]
        order_tilt = tools.interp(
            order_x_subpix, order_tilt, order_x_subpix_new, "extrapolate", "cubic")

        # reference sub-pixel map
        order_y_subpix = np.round(order_y_subpix)
        order_len_wavelength_subpix = len(order_wavelength_subpix)
        number_of_lambda[i] = order_len_wavelength_subpix

        # Load wavematrix in meters x order
        spectrograph_wavematrix_order_i = unit_converter.wavelength(
            spectrograph.wavematrix[i], "A", "m")

        # Efficiency interpolation and application
        order_total_efficiency_sky = np.zeros((order_len_wavelength_subpix, 2))
        order_total_efficiency_sky.T[0] = order_wavelength_subpix
        order_total_efficiency_sky.T[1] = tools.interp(
            spectrograph_wavematrix_order_i, acquisition.sky.sky_efficiency_matrix[i], order_wavelength_subpix, "extrapolate", "cubic")

        order_total_efficiency_object = np.zeros(
            (order_len_wavelength_subpix, 2))
        order_total_efficiency_object.T[0] = order_wavelength_subpix
        order_total_efficiency_object.T[1] = tools.interp(
            spectrograph_wavematrix_order_i, acquisition.sed.sed_total_efficincy[i], order_wavelength_subpix, "extrapolate", "cubic")

        # ---------------------------------------------------------------------
        # Initializing PSF interpolated maps
        psf_map = np.zeros((parameter.psf_field_sampling,
                           parameter.psf_field_sampling, order_len_wavelength_subpix))

        # Iter for each grid-point to interpolate in lambda

        order_psf_map_cube = spectrograph.psf_map[i][1]  # [1] -> .Cubes

        for x in range(0, parameter.psf_field_sampling):
            for y in range(0, parameter.psf_field_sampling):
                psf_data_xy = order_psf_map_cube[x][y]
                psf_map[x][y] = tools.interp(
                    wavelength, psf_data_xy, order_wavelength_subpix, "extrapolate")  # in subpix

        # ---------------------------------------------------------------------
        # Effective slit length / height and width/sampling_x

        order_shape = order.shape
        order_center_index = int(np.round(order_shape[0])/2)
        order_efficiency_subpix = np.round(
            order[order_center_index][6]) * parameter.pixel_oversampling

        ps_y_fact[i] = order_efficiency_subpix / \
            (48 * parameter.pixel_oversampling)
        sx = np.round(order_sx_subpix * parameter.pixel_oversampling)

        # Slit Image Simulation

        # ---------------------------------------------------------------------
        # Vectors and matrices initializing
        object_counts = np.zeros((order_len_wavelength_subpix, 1))
        sky_counts = np.zeros((order_len_wavelength_subpix, 1))
        F2_mat = np.zeros((order_len_wavelength_subpix, 101))
        N_eff = np.zeros((order_len_wavelength_subpix, 1))
        N_eff2 = np.zeros((order_len_wavelength_subpix, 1))
        NC_eff = np.zeros((order_len_wavelength_subpix, 1))

        order_detector_subpixel = np.zeros(
            (310*parameter.pixel_oversampling, spectrograph.n_pixels_subpixel + spectrograph.subpixel_edge))

        order_y_subpix_min[i] = min(order_y_subpix)
        order_y_subpix = order_y_subpix - \
            (order_y_subpix_min[i]) + 1 + \
            (3 * spectrograph.psf_map_pixel_number_subpixel)

        # ---------------------------------------------------------------------
        # IMPLEMENT PROGRESS BAR
        # tic
        # TEMP to sim few pixel
        v1 = np.array([1, 500, 1000, 1500, 1700])*parameter.pixel_oversampling
        print(v1)
        for j in v1:  # in order_len_wavelength_subpix
            object_counts[j] = tools.integration(
                order_wavelength_subpix[j], delta_lambda_subpix[j], acquisition.sed)

            sky_counts[j] = tools.integration(
                order_wavelength_subpix[j], delta_lambda_subpix[j], acquisition.sky.radiance)

            print(sky_counts[j])

    print("End Order: " + str(i))
    print("--- %s seconds ---" % (time.time() - order_time))
