from .simulation import Configuration
from .tool import tools
from .tool import efficiency
from .tool import unit_converter
from .tool import generic as plt
import numpy as np
import math
import time


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

    plt("SED", sed_wavelength, sed_flux, [
        "wavelength [$\AA$]", "flux " + r"[Ph/s/cm$^2$/$\AA$]"])

    sed_wavelength, sed_flux = acquisition.sed.normalize()

    plt("SED Normalized", sed_wavelength, sed_flux, [
        "wavelength [$\AA$]", "flux " + r"[Ph/s/cm$^2$/$\AA$]"])

    # 2nd step: get sky radiance and transission
    # [-] ,[ph/s/cm2/A]
    transimission, radiance = acquisition.sky.get_sky(
        acquisition.characteristics.slit_size_x, acquisition.characteristics.slit_size_y)

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

    spectrograph.set_subpixels(parameter.pixel_oversampling)

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
    refy_vect_min = np.zeros((spectrograph.len_n_orders))
    psf_bin_mat = np.zeros(
        (spectrograph.n_pixels_sub, spectrograph.n_pixels_sub, spectrograph.len_n_orders))

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
        wavelength = tools.interp(
            order_x_subpix, wavelength, order_x_subpix_new, "extrapolate", "cubic")
        delta_lambda = np.diff(wavelength)

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
        order_len_wavelength = len(wavelength)
        number_of_lambda[i] = order_len_wavelength

        # Load wavematrix in meters x order
        spectrograph_wavematrix_order_i = unit_converter.wavelength(
            spectrograph.wavematrix[i], "A", "m")

        # Efficiency interpolation and application
        order_total_efficiency_sky = np.zeros((order_len_wavelength, 2))
        order_total_efficiency_sky.T[0] = wavelength
        order_total_efficiency_sky.T[1] = tools.interp(
            spectrograph_wavematrix_order_i, acquisition.sky.sky_efficiency_matrix[i], wavelength, "extrapolate", "cubic")

        order_total_efficiency_object = np.zeros((order_len_wavelength, 2))
        order_total_efficiency_object.T[0] = wavelength
        order_total_efficiency_object.T[1] = tools.interp(
            spectrograph_wavematrix_order_i, acquisition.sed.sed_total_efficincy[i], wavelength, "extrapolate", "cubic")

        plt("Obj total efficiency", order_total_efficiency_object.T[0], order_total_efficiency_object.T[1], [
            "wavelength [$\AA$]", "[-]"])

        plt("Sky total efficiency", order_total_efficiency_sky.T[0], order_total_efficiency_sky.T[1], [
            "wavelength [$\AA$]", "[-]"])

        # CONROLLA TUTTI I PARAMETRI PRIMA

        '''
        jj=j-m_vect(end)+1;
        Tot_Eff_Obj_1ord(:,2)=interp1(Wave_mat(jj,:),Tot_Eff_Obj(jj,:),lambda, 'spline','extrap');
        Tot_Eff_1ord(:,2)=interp1(Wave_mat(jj,:),Tot_Eff(jj,:),lambda, 'spline','extrap');
        '''
