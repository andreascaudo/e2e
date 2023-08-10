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
    start_time = time.time()

    # Unpack configuration
    output = configuration.output
    telescope = configuration.telescope
    calibration = configuration.calibration
    spectrograph = configuration.spectrograph
    zemax = configuration.zemax
    acquisition = configuration.acquisition
    parameter = configuration.parameters

    if zemax is not None:
        spectrograph.grating.get_orders_wavelegth_range = zemax.get_orders_wavelegth_range(
                spectrograph.grating.n_orders, spectrograph.n_pixels, spectrograph.dimension_pixel)
        #print("Order Wavelegth Range: ", spectrograph.grating.get_orders_wavelegth_range)
        if zemax.order_table_flag:
            spectrograph.grating.order_table = zemax.get_order_table(
                spectrograph.grating.n_orders, spectrograph.n_pixels, spectrograph.dimension_pixel, parameter.psf_map_pixel_number)
        #print("Order Table: ", spectrograph.grating.order_table[0])
        #It should already take in account the new order table
        if zemax.PSF_map_flag:
            spectrograph.psf_map = zemax.get_PSF_map(spectrograph.grating.n_orders, parameter.psf_field_sampling)

    # TBI: Implement a function TO CHECK if acquisition reflects the spectrograph parameters

    if TIME:
        print("Start generating SED")
        delta_time = time.time()
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

    if TIME:
        print("SED generated in %s seconds" % (time.time() - delta_time))
        print("Start generating Sky and Transmission")
        delta_time = time.time()

    # 2nd step: get sky radiance and transission
    # [-] ,[ph/s/cm2/A]

    if len(spectrograph.slices) > 1:
        sky_slit_y = acquisition.characteristics.slit_size_y / \
            len(spectrograph.slices)
    else:
        sky_slit_y = acquisition.characteristics.slit_size_y

    transimission, radiance = acquisition.sky.get_sky(
        acquisition.characteristics.slit_size_x, sky_slit_y)
    if DEBUG:
        plot.generic("Sky transmission", transimission.wavelength,
                     transimission.transmission, ["wavelength [$\AA$]", "[-]"])

        if radiance.flux is not None:
            plot.generic("Sky radiance", radiance.wavelength,
                         radiance.flux, ["wavelength [$\AA$]", "[ph/s/cm2/A]"])

        for i in range(0, spectrograph.len_n_orders):
            plot.generic("Spectrograph + Telesciope Efficiency",
                         spectrograph.wavematrix[i], spectrograph.telescope_spectrograph_efficiency_fdr[i], ["wavelength [$\AA$]", "[-]"])

    if TIME:
        print("Sky and Transmission generated in %s seconds" %
              (time.time() - delta_time))
        print("Start setting Sky and Obj Efficiency")
        delta_time = time.time()
    # 3th step: Set Sky and Obj Efficiency
    acquisition.sky.set_efficiency(spectrograph.wavematrix,
                                   spectrograph.telescope_spectrograph_efficiency_fdr)
    acquisition.sed.set_efficiency(acquisition.sky.transmission.transmission_matrix,
                                   spectrograph.wavematrix, spectrograph.telescope_spectrograph_efficiency_fdr)

    # 4th step: Get Slit Efficiency and Image Quality
    # slit_efficiency_matrix, fwhm_iq_matrix = efficiency.get_slit_efficiency(spectrograph.wavematrix, acquisition.sky.airmass,
    #                                                                        acquisition.characteristics.slit_size_x, acquisition.characteristics.slit_size_y,
    #                                                                        acquisition.sky.seeing, spectrograph.fwhm_instrument, (telescope.diameter/100), telescope.l_zero)

    if TIME:
        print("Sky and Obj Efficiency set in %s seconds" %
              (time.time() - delta_time))
        print("Start setting subpixels and orders/slices")
        delta_time = time.time()

    spectrograph.set_subpixels(
        parameter.pixel_oversampling, parameter.psf_map_pixel_number)

    order_y_subpix_min = np.zeros(
        (spectrograph.len_n_orders, spectrograph.len_n_slices))

    if parameter.orders_index != None:
        orders = np.array(parameter.orders_index)
    else:
        orders = np.arange(0, spectrograph.len_n_orders)

    # TBD: Selection of slices from the configuration file
    slices = np.arange(0, spectrograph.len_n_slices)

    orders_slices = np.array(np.meshgrid(orders, slices)).T.reshape(-1, 2)

    if TIME:
        print("Subpixels and orders/slices set in %s seconds" %
              (time.time() - delta_time))
        print("Start Calculations")

    # order_y_subpix_min = do_in_parallel(orders, configuration)
    if PARALLEL:
        print("Start Parallel Calculation after ",
              time.time() - start_time, " seconds")
        with futures.ProcessPoolExecutor() as executor:
            parallel_results = executor.map(
                slice_calculation, orders_slices, repeat(configuration))

            for i, result in enumerate(parallel_results):
                spectrograph.detector_subpixel[:, :,
                                               result[2], result[3]] = result[0]
                order_y_subpix_min[result[2], result[3]] = result[1]
    else:
        print("Start Calculations after ", time.time() - start_time, " seconds")
        for order_slice in (orders_slices):
            spectrograph.detector_subpixel[:, :, order_slice[0], order_slice[1]], order_y_subpix_min[order_slice[0],
                                                                                                     order_slice[1]], i, j = slice_calculation(order_slice, configuration)

        # To keep just if i need to check with older version
        # for order in orders:
        #    spectrograph.detector_subpixel[:, :, order, 0], order_y_subpix_min[order, 0], i, j = order_calculation(
        #        order, configuration)

 # -------------------------------------------------------------------------
    # Recombination of the orders
    # When if pixel_oversampling < 10, the detector_recombined is 310x:
    # To simulate the AFC and not get MEMORY ERROR in the recombination of the orders
    # TBD: Check if this is the best way to do it
    if parameter.pixel_oversampling < 10:
        detector_recombined = np.zeros(
            (spectrograph.n_pixels_subpixel, spectrograph.n_pixels_subpixel))
    else:
        detector_recombined = np.zeros(
            (310*parameter.pixel_oversampling, spectrograph.n_pixels_subpixel))

    unkown3 = 3 if spectrograph.name == "SOXS" else 5

    for t in orders:
        for s in slices:
            print("ORDER:")
            print(order_y_subpix_min[t][s])
            if parameter.pixel_oversampling < 10:
                yy_start = int(order_y_subpix_min[t][s] -
                               (unkown3*spectrograph.psf_map_pixel_number_subpixel)) - 1
                yy_end = int(yy_start + (310*parameter.pixel_oversampling))
            else:
                yy_start = 0
                yy_end = 310*parameter.pixel_oversampling

            print("Start: ", str(yy_start), "End: ", str(yy_end))

            detector_recombined[yy_start:yy_end, :] = detector_recombined[yy_start:yy_end, :] + \
                spectrograph.detector_subpixel[:, int(spectrograph.subpixel_edge/2):int(
                    spectrograph.subpixel_edge/2)+spectrograph.n_pixels_subpixel, t, s]
            '''
            # RIX
            detector_recombined[yy_start:yy_end, :] = detector_recombined[yy_start:yy_end, :] + \
                spectrograph.detector_subpixel[:,
                                               :spectrograph.n_pixels_subpixel, t, s]
            '''

    # Binnig to the final detector image with real size
    detector_recombined_binned = tools.rebin_image(
        detector_recombined, [parameter.pixel_oversampling, parameter.pixel_oversampling])

    # Display the final detector image [Binned]

    # Moltilpy for telescope are to obtain [/s]
    if calibration is None and acquisition.sed.calibration is False:
        M2_obs = 0.03 if spectrograph.name == "CUBES" else 0
        detector_final = detector_recombined_binned * \
            ((telescope.diameter/2) ** 2) * math.pi * (1-M2_obs)
    else:
        detector_final = detector_recombined_binned

    # Moltiply for time exposure
    detector_final = detector_final * \
        acquisition.characteristics.detector_integration_time

    if spectrograph.name == "CUBES":
        # SHIFT
        if parameter.pixel_oversampling < 10:
            Det_s = 704
            detector_shifted = np.zeros(
                (detector_final.shape[0], detector_final.shape[1]))
            detector_shifted[:-Det_s, :] = detector_final[Det_s:, :]
        else:
            detector_shifted = detector_final

        # Piston Counts w.r.t. to Bias-Dark
        piston = 0  # 50
        detector_final = detector_shifted + piston
        # Add pre/over scan region
        # detector_final = tools.add_pre_over_scan(detector_final)

    # --- END OF THE SIMULATION ------------------------------------------------

    print("End Calculations")
    print("--- %s seconds ---" % (time.time() - start_time))

    # plot.detector('Detector - Real Pixel Scale - PSF Incl' + ' - DIT = ' +
    #              str(acquisition.characteristics.detector_integration_time) + 's - [e-]', detector_final)

    # Get current time and use it to save the fits file
    if output.SAVE:
        current_time = time.strftime("%Y%m%d-%H%M%S")
        folder = output.output_folder
        filename = current_time + '.fits'
        fits.writeto(path.join(folder, filename),
                     detector_final, overwrite=True)


def slice_calculation(order_slice, configuration):
    i = order_slice[0]
    slice_index = order_slice[1]

    # Unpack configuration
    telescope = configuration.telescope
    spectrograph = configuration.spectrograph
    calibration = configuration.calibration
    acquisition = configuration.acquisition
    parameter = configuration.parameters

    # ---------------------------------------------------------------------
    # slice and single-order-table
    order_number = spectrograph.grating.n_orders[i]
    slice_number = spectrograph.slices[slice_index].slice_id

    print("Order Number: ", str(order_number))
    print("Slice Number: ", str(slice_number))
    slice_time = time.time()

    order = spectrograph.grating.order_table[spectrograph.grating.order_table.T[0] == order_number]
    slice = order[spectrograph.grating.order_table.T[1] ==
                  slice_number] if len(spectrograph.slices) > 1 else order

    if len(spectrograph.slices) > 1:
        # Slice Index Table (SIT), due to the fact that the order table without slice has one column less
        SIT = 1
    else:
        SIT = 0

    # wavlength vector in meters
    wavelength = unit_converter.wavelength(
        slice.T[2+SIT], "um", "m")  # in meters
    # ---------------------------------------------------------------------

    # ---------------------------------------------------------------------
    # ### Sub-Pixel scale: x, lambda, y, Sx
    # X
    order_x_subpix = slice.T[3+SIT] * parameter.pixel_oversampling
    order_x_subpix_start = math.ceil(
        slice.T[3+SIT][0] * parameter.pixel_oversampling)
    order_x_subpix_end = math.ceil(
        slice.T[3+SIT][-1] * parameter.pixel_oversampling)
    order_x_subpix_new = np.arange(
        order_x_subpix_start, order_x_subpix_end+1, 1)

    # Lambda
    order_wavelength_subpix = tools.interp(
        order_x_subpix, wavelength, order_x_subpix_new, "extrapolate", "cubic")
    delta_lambda_subpix = np.diff(order_wavelength_subpix)

    # Y
    order_y_subpix = slice.T[4+SIT] * parameter.pixel_oversampling
    order_y_subpix = tools.interp(
        order_x_subpix, order_y_subpix, order_x_subpix_new, "extrapolate", "cubic")

    # Sx
    order_sx = slice.T[5+SIT]
    order_sx = tools.interp(
        order_x_subpix, order_sx, order_x_subpix_new, "extrapolate", "cubic")

    # DELATE PSF_XX PSF_YY BEFORE CONTNUING
    # Tilt Order
    order_tilt = slice.T[11+SIT]
    order_tilt = tools.interp(
        order_x_subpix, order_tilt, order_x_subpix_new, "extrapolate", "cubic")

    # reference sub-pixel map
    order_x_subpix = order_x_subpix_new
    order_y_subpix = np.round(order_y_subpix)
    order_len_wavelength_subpix = len(order_wavelength_subpix)

    # ---------------------------------------------------------------------
    # Load wavematrix in meters x order
    spectrograph_wavematrix_order_i = unit_converter.wavelength(
        spectrograph.wavematrix[i], "A", "m")
    # Efficiency interpolation and application
    # SKY
    order_total_efficiency_sky = np.zeros((order_len_wavelength_subpix, 2))
    order_total_efficiency_sky.T[0] = order_wavelength_subpix
    order_total_efficiency_sky.T[1] = tools.interp(
        spectrograph_wavematrix_order_i, acquisition.sky.sky_efficiency_matrix[i], order_wavelength_subpix, "extrapolate", "cubic")
    # OBJ
    order_total_efficiency_object = np.zeros(
        (order_len_wavelength_subpix, 2))
    order_total_efficiency_object.T[0] = order_wavelength_subpix
    order_total_efficiency_object.T[1] = tools.interp(
        spectrograph_wavematrix_order_i, acquisition.sed.sed_total_efficincy[i], order_wavelength_subpix, "extrapolate", "cubic")

    # ---------------------------------------------------------------------
    # Initializing and interpolate PSF map
    psf_map_shape = ((parameter.psf_field_sampling,
                      parameter.psf_field_sampling, order_len_wavelength_subpix))

    # Extract psf map where spectrograph.psf_map[i][0] is equal to order number
    if spectrograph.psf_map[i][0] == order_number:
        if spectrograph.psf_map[i][1][slice_index][0] == slice_number:
            psf_map = tools.interpolate_psf_map(
                psf_map_shape, spectrograph.psf_map[i][1][slice_index][1], wavelength, order_wavelength_subpix)
        else:
            raise Exception("PSF MAP NOT FOUND:\nSimulation Slice: " + str(slice_number) + "does not match with the PSF map slice: " +
                            str(spectrograph.psf_map[i][1][slice_index][0]))
    else:
        raise Exception("PSF MAP NOT FOUND:\nSimulation Order: " + str(order_number) + "does not match with the PSF map order: " +
                        str(spectrograph.psf_map[i][0]))

    # ---------------------------------------------------------------------
    # Effective slit length / height and width/sampling_x

    slice_shape = slice.shape
    slice_center_index = int(np.round(slice_shape[0])/2)
    order_efficiency_subpix = np.round(
        slice[slice_center_index][6+SIT]) * parameter.pixel_oversampling  # Should be called as length

    uknown5 = 48 if spectrograph.name == "SOXS" else 100

    ps_y_fact = order_efficiency_subpix / \
        (uknown5 * parameter.pixel_oversampling)
    sx = order_sx * parameter.pixel_oversampling

    # Slit Image Simulation

    # ---------------------------------------------------------------------
    # Vectors and matrices initializing
    object_counts = np.zeros((order_len_wavelength_subpix, 1))
    sky_counts = np.zeros((order_len_wavelength_subpix, 1))
    object_efficiency = np.zeros((order_len_wavelength_subpix, 1))
    sky_efficiency = np.zeros((order_len_wavelength_subpix, 1))

    order_detector_subpixel = np.zeros(
        (310*parameter.pixel_oversampling, spectrograph.n_pixels_subpixel + spectrograph.subpixel_edge))

    order_y_subpix_min = min(order_y_subpix)
    unkown = 3 if spectrograph.name == "SOXS" else 5
    order_y_subpix = order_y_subpix - \
        (order_y_subpix_min) + 1 + \
        (unkown * spectrograph.psf_map_pixel_number_subpixel)

    # ----------------------------------------------------
    '''
    features_positions_sim = {
                          23: [[273, 1450], [1204, 1350], [1396, 1355]]
                        }
    '''
    # temp_pixel = 1000 if spectrograph.name == "SOXS" else 2955
    # rng = np.array([273, 1204, 1396])*parameter.pixel_oversampling
    # rng = range(2640 * parameter.pixel_oversampling,
    #            4453 * parameter.pixel_oversampling)
    rng = range(0, order_len_wavelength_subpix-1)
    print("SEEING: ", acquisition.sky.seeing)

    for j in tqdm(rng):
        subpix_time = time.time()
        time_save = time
        if TIME:
            print("Starting Slice: ", slice_index, " sub-pixel: ", j)
            print("Start init")
            subpix_time = time.time()
            subpix_time_save = subpix_time

        # Calculate plate scale
        if len(spectrograph.slices) > 1:
            pixel_slice = order_sx[j]  # pixel/slice
            slice_arcsec = acquisition.characteristics.slit_size_x / \
                len(spectrograph.slices)  # arcsec/slice
            pixel_platescale = pixel_slice/slice_arcsec  # pixel/arcsec
        else:
            pixel_slit = order_sx[j]  # pixel/slit
            slit_arcsec = acquisition.characteristics.slit_size_x  # arcsec/slit
            pixel_platescale = pixel_slit/slit_arcsec  # pixel/arcsec

        # Obj Counts & Efficiency
        object_counts[j] = tools.integration(
            order_wavelength_subpix[j], delta_lambda_subpix[j], acquisition.sed)

        object_efficiency[j] = object_counts[j] * \
            order_total_efficiency_object.T[1][j]

        # Sky Counts & Efficiency
        if acquisition.sky.radiance.flux is None:
            sky_counts[j] = 0
        else:
            sky_counts[j] = tools.integration(
                order_wavelength_subpix[j], delta_lambda_subpix[j], acquisition.sky.radiance)
            sky_efficiency[j] = sky_counts[j] * \
                order_total_efficiency_sky.T[1][j]

        if object_counts[j] == 0 and sky_counts[j] == 0:
            continue

        # --- Image size --------------------------------------------------
        uknown2 = 6 if spectrograph.name == "SOXS" else 10
        image_size = np.array([uknown2*spectrograph.psf_map_pixel_number_subpixel,
                              int(2*spectrograph.psf_map_pixel_number_subpixel)])
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
                    pixel_platescale,
                    spectrograph.dimension_pixel,
                    acquisition.sky.seeing,
                    parameter.pixel_oversampling,
                    image_size, ps_y_fact)
        else:
            calibration_slit = tools.calibration_slit(
                object_efficiency[j],
                image_size,
                calibration.area,
                ps_y_fact,
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
            print("Start shift")
            subpix_time = time.time()

            # --- SHIFT ---- to be applied -----------------------
        if len(spectrograph.slices) > 1:
            pixel_oversampling_shift = parameter.mask_oversampling
        else:
            pixel_oversampling_shift = 1

        if len(spectrograph.slices) > 1:
            subpixel_slice = sx[j]  # subpixel/slice
            slice_arcsec = acquisition.characteristics.slit_size_x / \
                len(spectrograph.slices)  # arcsec/slice
            subpixel_platescale = subpixel_slice/slice_arcsec  # subpixel/arcsec

            spectrograph.slices[slice_index].to_pix_arcsec(subpixel_platescale)
            spectrograph.slices[slice_index].to_subpix(
                pixel_oversampling_shift)  #  subpixel to subpixel_mask

        # --- SHIFT ---------------------------------------------
        if len(spectrograph.slices) > 1:
            shift = int(
                np.round(spectrograph.slices[slice_index].shift_subpix))

            vs0_1 = np.linspace(1, image_size[0], image_size[0])
            vs0_2 = np.linspace(1, image_size[1], image_size[1])

            vs1 = np.linspace(1, image_size[0],
                              image_size[0]*pixel_oversampling_shift)
            vs2 = np.linspace(1, image_size[1],
                              image_size[1]*pixel_oversampling_shift)
            '''
            Xs0, Ys0 = np.meshgrid(vs0_2, vs0_1)
            Xs1, Ys1 = np.meshgrid(vs2, vs1)

            psf_interp_0 = interpolate.griddata(
                (Xs0.flatten(), Ys0.flatten()), detector.flatten(), (Xs1, Ys1))
            psf_interp_0 = psf_interp_0 * (vs1[1] - vs1[0]) * (vs2[1] - vs2[0])
            '''
            # GRidata is SLOW and interp2d is DEPRECATED
            # It took 0.2 seconds in contrast to 0.0005 seconds with interpn
            # The difference in accuracy is negligible, the maximum difference is 0.00010 while the average is 2.3e-8
            # The difference in speed is huge, interp2d is 400 times faster
            # So I decided to use interpn instead of griddata

            Xs1, Ys1 = np.meshgrid(vs1, vs2, indexing="ij")

            # Interpolate with interpn
            detector = interpolate.interpn(
                (vs0_1, vs0_2), detector, (Xs1, Ys1), method="linear", bounds_error=False)

            detector = detector * (vs1[1] - vs1[0]) * (vs2[1] - vs2[0])

            if shift < 0:
                detector[:, (-shift):] = detector[:, 0:+shift]
            elif shift > 0:
                detector[:, 0:-shift] = detector[:, shift:]

        if TIME:
            print("Ended shifting in:")
            print("--- %s seconds ---" % (time.time() - subpix_time))
            print("Start mask")
            subpix_time = time.time()

        # --- MASK --------------------------------------------------------
        if calibration is None:
            sy_m = order_efficiency_subpix * pixel_oversampling_shift
            sx_m = sx[j] * pixel_oversampling_shift

            mask = tools.mask_ideal_slit(
                image_size, sy_m, sx_m, pixel_oversampling_shift)
        else:
            # Pinhole
            mask = calibration.get_mask(
                ps_y_fact, image_size, parameter.pixel_oversampling, parameter.mask_oversampling)

        # HERE
        detector = detector * mask

        if pixel_oversampling_shift > 1:
            detector = tools.rebin_image(
                detector, [pixel_oversampling_shift, pixel_oversampling_shift])

        if TIME:
            print("Ended mask in:")
            print("--- %s seconds ---" % (time.time() - subpix_time))
            print("Start tilting")
            subpix_time = time.time()

        # ----- FLIP SLIT(potrebbe essere che basta ruotare senza fare resampling).
        # detector_rotated = tools.flip_slit(detector, order_tilt[j])
        if spectrograph.name == "CUBES":
            tilt = order_tilt[j]
        else:
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

        if DEBUG:
            plot.detector("Order: " + str(i) + ": Slit Image",
                          detector_rotated)

        if TIME:
            print("Ended normalize and rebin convoluting in:")
            print("--- %s seconds ---" % (time.time() - subpix_time))
            print("Start convoluting")
            subpix_time = time.time()

        detector_conv = tools.convolve(detector_rotated, psf_bin)

        if TIME:
            print("Ended convoluting in:")
            print("--- %s seconds ---" % (time.time() - subpix_time))
            print("Start Positioning the image")
            subpix_time = time.time()

        # --- POSITION IMAGE SIMULATION ELEMENT ---------------------------
        ref_x = order_x_subpix[j] + (spectrograph.subpixel_edge/2)
        ref_y = order_y_subpix[j]
        ref_y_start = int(ref_y - (detector_conv.shape[0]/2)) - 1
        ref_y_end = int(ref_y_start + detector_conv.shape[0])

        ref_x_start = int(
            ref_x - (detector_conv.shape[1]/2)) + int(parameter.pixel_oversampling/2) - 1
        ref_x_end = int(ref_x_start + detector_conv.shape[1])  # Nessun +1

        order_detector_subpixel[ref_y_start:ref_y_end,
                                ref_x_start: ref_x_end] = order_detector_subpixel[ref_y_start: ref_y_end, ref_x_start: ref_x_end] + detector_conv

        if TIME:
            print("Ended Positioning the image in:")
            print("--- %s seconds ---" % (time.time() - subpix_time))
            print("Ended Subpixel: " + str(j))
            print("--- %s seconds ---" % (time.time() - subpix_time_save))

    # spectrograph.detector_subpixel[:, :, i] = order_detector_subpixel

    print("End Slice: " + str(i))
    print("--- %s seconds ---" % (time.time() - slice_time))

    # Display the final detector image
    if DEBUG:
        plot.detector("Order: " + str(i),
                      spectrograph.detector_subpixel[:, :, i])

    return order_detector_subpixel, order_y_subpix_min, i, slice_index
