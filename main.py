from scipy import ndimage, interpolate

from .simulation import Configuration
from .tool import tools
from .tool import efficiency
from .tool import unit_converter
from .tool import plot
import numpy as np
import math
import time
from astropy.io import fits
from concurrent import futures
from itertools import repeat
from tqdm import tqdm
from os import path

DEBUG = False
TIME = False
PARALLEL = False


def run(configuration: Configuration):
    # Unpack configuration
    output = configuration.output
    telescope = configuration.telescope
    calibration = configuration.calibration
    spectrograph = configuration.spectrograph
    acquisition = configuration.acquisition
    parameter = configuration.parameters

    # TBI: Implement a function TO CHECK if acquisition reflects the spectrograph parameters

    # 1st step: generate flux
    # [Angstrom], [Ph/s/cm^2/A]
    sed_wavelength, sed_flux = acquisition.sed.get_flux()

    if DEBUG:
        plot.generic("SED", sed_wavelength, sed_flux, [
            "wavelength [$\AA$]", "flux " + r"[Ph/s/cm$^2$/$\AA$]"])

    if not acquisition.sed.calibration:
        sed_wavelength, sed_flux = acquisition.sed.normalize()

    if DEBUG:
        plot.generic("SED Normalized", sed_wavelength, sed_flux, [
            "wavelength [$\AA$]", "flux " + r"[Ph/s/cm$^2$/$\AA$]"])

    # 2nd step: get sky radiance and transission
    # [-] ,[ph/s/cm2/A]
    transimission, radiance = acquisition.sky.get_sky(
        acquisition.characteristics.slit_size_x, acquisition.characteristics.slit_size_y)
    if DEBUG:
        plot.generic("Sky transmission", transimission.wavelength,
                     transimission.transmission, ["wavelength [$\AA$]", "[-]"])

        plot.generic("Sky radiance", radiance.wavelength,
                     radiance.flux, ["wavelength [$\AA$]", "[ph/s/cm2/A]"])

        for i in range(0, spectrograph.len_n_orders):
            plot.generic("Spectrograph + Telesciope Efficiency",
                         spectrograph.wavematrix[i], spectrograph.telescope_spectrograph_efficiency_fdr[i], ["wavelength [$\AA$]", "[-]"])

    # 3th step: Set Sky and Obj Efficiency
    acquisition.sky.set_efficiency(spectrograph.wavematrix,
                                   spectrograph.telescope_spectrograph_efficiency_fdr)
    acquisition.sed.set_efficiency(acquisition.sky.transmission.transmission_matrix,
                                   spectrograph.wavematrix, spectrograph.telescope_spectrograph_efficiency_fdr)

    # 4th step: Get Slit Efficiency and Image Quality
    # slit_efficiency_matrix, fwhm_iq_matrix = efficiency.get_slit_efficiency(spectrograph.wavematrix, acquisition.sky.airmass,
    #                                                                        acquisition.characteristics.slit_size_x, acquisition.characteristics.slit_size_y,
    #                                                                        acquisition.sky.seeing, spectrograph.fwhm_instrument, (telescope.diameter/100), telescope.l_zero)

    spectrograph.set_subpixels(
        parameter.pixel_oversampling, parameter.psf_map_pixel_number)

    order_y_subpix_min = np.zeros((spectrograph.len_n_orders))

    if parameter.orders_index != None:
        orders = np.array(parameter.orders_index)
    else:
        orders = np.arange(0, spectrograph.len_n_orders)
    #orders = np.array([15, 14, 7, 6, 1, 0])
    #orders = np.array([15])

    print("Start Calculations")
    start_time = time.time()

    # order_y_subpix_min = do_in_parallel(orders, configuration)
    if PARALLEL:
        with futures.ProcessPoolExecutor() as executor:
            parallel_results = executor.map(
                calculation, orders, repeat(configuration))

            for i, result in enumerate(parallel_results):
                spectrograph.detector_subpixel[:, :, result[2]] = result[0]
                order_y_subpix_min[result[2]] = result[1]
    else:
        for order in orders:
            spectrograph.detector_subpixel[:, :, order], order_y_subpix_min[order], i = calculation(
                order, configuration)

 # -------------------------------------------------------------------------
    # Recombination of the orders

    detector_recombined = np.zeros(
        (spectrograph.n_pixels_subpixel, spectrograph.n_pixels_subpixel))

    for t in orders:
        yy_start = int(order_y_subpix_min[t] -
                       (3*spectrograph.psf_map_pixel_number_subpixel)) - 1
        yy_end = int(yy_start + (310*parameter.pixel_oversampling))

        detector_recombined[yy_start:yy_end, :] = detector_recombined[yy_start:yy_end, :] + \
            spectrograph.detector_subpixel[:, int(spectrograph.subpixel_edge/2):int(
                spectrograph.subpixel_edge/2)+spectrograph.n_pixels_subpixel, t]

    # Binnig to the final detector image with real size

    detector_recombined_binned = tools.rebin_image(
        detector_recombined, [parameter.pixel_oversampling, parameter.pixel_oversampling])

    # Display the final detector image [Binned]

    # Moltilpy for telescope are to obtain [/s]
    if calibration is None:
        detector_final = detector_recombined_binned * \
            ((telescope.diameter/2) ** 2) * math.pi
    else:
        detector_final = detector_recombined_binned

    # Moltiply for time exposure
    detector_final = detector_final * \
        acquisition.characteristics.detector_integration_time

    plot.detector('Detector - Real Pixel Scale - PSF Incl' + ' - DIT = ' +
                  str(acquisition.characteristics.detector_integration_time) + 's - [e-]', detector_final)

    # Get current time and use it to save the fits file
    if output.SAVE:
        current_time = time.strftime("%Y%m%d-%H%M%S")
        folder = output.output_folder
        filename = current_time + '.fits'
        fits.writeto(path.join(folder, filename),
                     detector_final, overwrite=True)

    # --- END OF THE SIMULATION ------------------------------------------------

    print("End Calculations")
    print("--- %s seconds ---" % (time.time() - start_time))


def calculation(i, configuration):
    # Unpack configuration
    telescope = configuration.telescope
    spectrograph = configuration.spectrograph
    calibration = configuration.calibration
    acquisition = configuration.acquisition
    parameter = configuration.parameters

    number_of_lambda = np.zeros((spectrograph.len_n_orders))
    ps_y_fact = np.zeros((spectrograph.len_n_orders))
    order_y_subpix_min = np.zeros((spectrograph.len_n_orders))
    psf_bin_matrix = np.zeros(
        (spectrograph.psf_map_pixel_number_subpixel, spectrograph.psf_map_pixel_number_subpixel, spectrograph.len_n_orders))

    # orders_loop = range(0, spectrograph.len_n_orders-1)
    # orders_loop = [4]
    # for i in tqdm(orders_loop):  # [4]
    print("Start Order: ", str(i))
    order_time = time.time()

    # ---------------------------------------------------------------------
    # order and single-order-table
    order_number = spectrograph.n_orders[i]
    print("Order Number: ", str(order_number))
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
    order_x_subpix = order_x_subpix_new
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
    # Initializing and interpolate PSF map
    psf_map_shape = ((parameter.psf_field_sampling,
                      parameter.psf_field_sampling, order_len_wavelength_subpix))
    psf_map = tools.interpolate_psf_map(
        psf_map_shape, spectrograph.psf_map[i][1], wavelength, order_wavelength_subpix)

    # ---------------------------------------------------------------------
    # Effective slit length / height and width/sampling_x

    order_shape = order.shape
    order_center_index = int(np.round(order_shape[0])/2)
    order_efficiency_subpix = np.round(
        order[order_center_index][6]) * parameter.pixel_oversampling  # Should be called as length

    ps_y_fact[i] = order_efficiency_subpix / \
        (48 * parameter.pixel_oversampling)
    sx = np.round(order_sx_subpix * parameter.pixel_oversampling)

    # Slit Image Simulation

    # ---------------------------------------------------------------------
    # Vectors and matrices initializing
    object_counts = np.zeros((order_len_wavelength_subpix, 1))
    sky_counts = np.zeros((order_len_wavelength_subpix, 1))
    object_efficiency = np.zeros((order_len_wavelength_subpix, 1))
    sky_efficiency = np.zeros((order_len_wavelength_subpix, 1))

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
    # v1 = np.array([1, 500, 1000, 1500, 1700]) * \
    #    parameter.pixel_oversampling
    # for j in tqdm(v1):

    #rng = range(0, order_len_wavelength_subpix-1)
    # rng = range(1050 * parameter.pixel_oversampling,
    #            1100 * parameter.pixel_oversampling)
    rng = np.array([1000])*parameter.pixel_oversampling
    #rng = np.arange(271, 288, 1)*parameter.pixel_oversampling

    for j in tqdm(rng):
        if TIME:
            print("Starting order: ", i, " sub-pixel: ", j)
            print("Start init")
            subpix_time = time.time()
            subpix_time_save = subpix_time

        # Obj Counts & Efficiency
        object_counts[j] = tools.integration(
            order_wavelength_subpix[j], delta_lambda_subpix[j], acquisition.sed)
        object_efficiency[j] = object_counts[j] * \
            order_total_efficiency_object.T[1][j]

        # Sky Counts & Efficiency
        sky_counts[j] = tools.integration(
            order_wavelength_subpix[j], delta_lambda_subpix[j], acquisition.sky.radiance)
        sky_efficiency[j] = sky_counts[j] * \
            order_total_efficiency_sky.T[1][j]

        # --- Image size --------------------------------------------------
        image_size = np.array([6*spectrograph.psf_map_pixel_number_subpixel,
                               2*spectrograph.psf_map_pixel_number_subpixel])
        if TIME:
            print("Ended Init in:")
            print("--- %s seconds ---" % (time.time() - subpix_time))
            print("Start obj slit")
            subpix_time = time.time()

        # --- OBJ + SKY slit ---------------------------------------------

        # OBJ
        if calibration is None:
            if acquisition.sed.calibration:
                object_slit = np.ones(image_size) * \
                    object_efficiency[j] / (order_efficiency_subpix * sx[j])
            else:
                object_slit = tools.object_slit(
                    object_efficiency[j],
                    telescope.f_number, telescope.pupil_equivalent_diameter,
                    spectrograph.dimension_pixel,
                    acquisition.sky.seeing,
                    parameter.pixel_oversampling,
                    image_size, ps_y_fact[i])
        else:
            calibration_slit = tools.calibration_slit(
                object_efficiency[j],
                image_size,
                calibration.area,
                ps_y_fact[i],
                acquisition.characteristics.slit_size_x,
                parameter.pixel_oversampling
            )

        if TIME:
            print("Ended obj slit in:")
            print("--- %s seconds ---" % (time.time() - subpix_time))
            print("Start sky slit")
            subpix_time = time.time()

        # SKY
        sky_slit = np.ones(image_size) * \
            sky_efficiency[j] / (order_efficiency_subpix * sx[j])

        if TIME:
            print("Ended sky slit in:")
            print("--- %s seconds ---" % (time.time() - subpix_time))

        # --- CCD --------------------------------------------------------
        if calibration is None:
            detector = object_slit + sky_slit
        else:
            detector = calibration_slit

        if TIME:
            print("Start mask")
            subpix_time = time.time()

        # --- MASK --------------------------------------------------------
        if calibration is None:
            sy_m = order_efficiency_subpix
            sx_m = sx[j]
            mask = tools.mask_ideal_slit(image_size, sy_m, sx_m)
        else:
            mask = calibration.get_mask(
                ps_y_fact[i], image_size, parameter.pixel_oversampling, parameter.mask_oversampling)
        # HERE
        detector = detector * mask

        if TIME:
            print("Ended mask in:")
            print("--- %s seconds ---" % (time.time() - subpix_time))
            print("Start tilting")
            subpix_time = time.time()

        # ----- FLIP SLIT(potrebbe essere che basta ruotare senza fare resampling).
        # detector_rotated = tools.flip_slit(detector, order_tilt[j])
        tilt = - order_tilt[j]
        # Rotate the detector by tilt
        detector_rotated = ndimage.rotate(detector, tilt, reshape=False)

        if TIME:
            print("Ended tilting in:")
            print("--- %s seconds ---" % (time.time() - subpix_time))
            print("Start interpolate for convoluting")
            subpix_time = time.time()

        # --- WHEN CONV TO BE DONE HERE ---------------------------------------
        # Extract the wave-interp map and normalization
        psf_map_j_norm, psf_box_z = tools.normalize_psf_map(
            psf_map[:, :, j])
        #  Initialize the convolution matrix
        v1, v2 = tools.init_conv_matrix(
            parameter.psf_map_pixel_number, spectrograph.dimension_pixel, psf_box_z)

        psf_interp, psf_interp_sum = tools.interpolate_griddata_psf_map(
            psf_map_j_norm, v1, v2)

        if TIME:
            print("Ended interpolate for convoluting in:")
            print("--- %s seconds ---" % (time.time() - subpix_time))
            print("Start normalize and rebin convoluting")
            subpix_time = time.time()

        # Normalize the PSF map
        psf_interp = psf_interp / psf_interp_sum

        # Resampling - rebinning factor : Num-SubPix of 1um size / Num-SubPix for PPP2
        rebin_factor = (parameter.psf_map_pixel_number * spectrograph.dimension_pixel) / (
            parameter.psf_map_pixel_number * parameter.pixel_oversampling)
        psf_bin = tools.rebin_image(
            psf_interp, [rebin_factor, rebin_factor])

        if TIME:
            print("Ended normalize and rebin convoluting in:")
            print("--- %s seconds ---" % (time.time() - subpix_time))
            print("Start convoluting")
            subpix_time = time.time()

        # Convolution
        # detector_conv = tools.convolve(
        #    detector_rotated, psf_bin)

        if DEBUG:
            plot.detector("Order: " + str(i) + ": Slit Image",
                          detector_rotated)

        detector_conv = tools.convolve(detector_rotated, psf_bin)

        # detector_conv = tools.convolve2D(detector_rotated, psf_bin) # DO NOT use "same" condition

        if TIME:
            print("Ended convoluting in:")
            print("--- %s seconds ---" % (time.time() - subpix_time))
            print("Start Positioning the image")
            subpix_time = time.time()

        # --- POSITION IMAGE SIMULATION ELEMENT ---------------------------

        ref_x = order_x_subpix[j] + (spectrograph.subpixel_edge/2)
        ref_y = order_y_subpix[j]
        ref_y_start = int(ref_y - (np.floor(detector_conv.shape[0]/2))) - 1
        ref_y_end = int(ref_y_start + detector_conv.shape[0])
        ref_x_start = int(ref_x - (np.floor(detector_conv.shape[1]/2))) - 1
        ref_x_end = int(ref_x_start + detector_conv.shape[1])

        order_detector_subpixel[ref_y_start:ref_y_end,
                                ref_x_start:ref_x_end] = order_detector_subpixel[ref_y_start:ref_y_end, ref_x_start:ref_x_end] + detector_conv

        if TIME:
            print("Ended Positioning the image in:")
            print("--- %s seconds ---" % (time.time() - subpix_time))
            print("Ended Subpixel: " + str(j))
            print("--- %s seconds ---" % (time.time() - subpix_time_save))

    spectrograph.detector_subpixel[:, :, i] = order_detector_subpixel
    psf_bin_matrix[:, :, i] = psf_bin

    print("End Order: " + str(i))
    print("--- %s seconds ---" % (time.time() - order_time))

    # Display the final detector image
    if DEBUG:
        plot.detector("Order: " + str(i),
                      spectrograph.detector_subpixel[:, :, i])

    return spectrograph.detector_subpixel[:, :, i], order_y_subpix_min[i], i
