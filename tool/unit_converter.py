# coding=utf-8
from scipy import constants

FLUX_UNIT_OF_MEASURE = ["Jy", "photons/cm^2/s/KeV", "photons/cm^2/s/A",
                        "ergs/cm^2/s/A", "ergs/cm^2/s/Hz", "W/m^2/um", "J/(s * m^3)", "W/m^3"]

WAVELENGTH_UNIT_OF_MEASURE = ["um", "nm", "A"]


def flux(flux, frm, to, E_LL=[]):
    if(frm not in FLUX_UNIT_OF_MEASURE):
        raise Exception(
            "ERROR: FROM - Unit of measure NOT present -> choose from [photons/cm^2/s/KeV, photons/cm^2/s/A, ergs/cm^2/s/A, ergs/cm^2/s/Hz]")

    if(to not in FLUX_UNIT_OF_MEASURE):
        raise Exception(
            "ERROR: TO - Unit of measure NOT present -> choose from [photons/cm^2/s/KeV, photons/cm^2/s/A, ergs/cm^2/s/A, ergs/cm^2/s/Hz]")

    if(frm == to):
        return flux

    # FROM Jy
    if(frm == FLUX_UNIT_OF_MEASURE[0]):
        # TO photons/cm^2/s/KeV
        if(to == FLUX_UNIT_OF_MEASURE[1]):
            if(E_LL != []):
                return 1.51 * 10**3 * flux / E_LL
            else:
                raise Exception(
                    "ERROR: FROM Jy TO photons/cm^2/s/keV -> Energy in keV needed, pass as -> flux_converter(flux, from, to, E)")
        # TO photons/cm^2/s/A
        if(to == FLUX_UNIT_OF_MEASURE[2]):
            if(E_LL != []):
                return 1.51 * 10**3 * flux / E_LL
            else:
                raise Exception(
                    "ERROR: FROM Jy TO photons/cm^2/s/A -> Wavelength in A needed, pass as -> flux_converter(flux, from, to, LL)")
        # TO ergs/cm^2/s/A
        if(to == FLUX_UNIT_OF_MEASURE[3]):
            if(E_LL != []):
                return 3.00 * 10**-5 * flux / E_LL**2
            else:
                raise Exception(
                    "ERROR: FROM Jy TO ergs/cm^2/s/A -> Wavelength in A needed, pass as -> flux_converter(flux, from, to, LL)")
        # TO ergs/cm^2/s/Hz
        if(to == FLUX_UNIT_OF_MEASURE[4]):
            return 10**-23 * flux

    # FROM photons/cm^2/s/KeV
    if(frm == FLUX_UNIT_OF_MEASURE[1]):
        # TO Jy
        if(to == FLUX_UNIT_OF_MEASURE[0]):
            if(E_LL != []):
                return 6.63 * 10**-4 * E_LL * flux
            else:
                raise Exception(
                    "ERROR: FROM photons/cm^2/s/KeV TO Jy -> Energy in keV needed, pass as -> flux_converter(flux, from, to, E)")
        # TO photons/cm^2/s/A
        if(to == FLUX_UNIT_OF_MEASURE[2]):
            if(E_LL != []):
                return 8.07 * 10**-2 * E_LL**2 * flux
            else:
                raise Exception(
                    "ERROR: FROM photons/cm^2/s/KeV TO photons/cm^2/s/A -> Energy in keV needed, pass as -> flux_converter(flux, from, to, E)")
        # TO ergs/cm^2/s/A
        if(to == FLUX_UNIT_OF_MEASURE[3]):
            if(E_LL != []):
                return 1.29 * 10**-10 * E_LL**3 * flux
            else:
                raise Exception(
                    "ERROR: FROM photons/cm^2/s/KeV TO ergs/cm^2/s/A -> Energy in keV needed, pass as -> flux_converter(flux, from, to, E)")
        # TO ergs/cm^2/s/Hz
        if(to == FLUX_UNIT_OF_MEASURE[4]):
            if(E_LL != []):
                return 6.63 * 10**-27 * E_LL * flux
            else:
                raise Exception(
                    "ERROR: FROM photons/cm^2/s/KeV TO ergs/cm^2/s/Hz -> Energy in keV needed, pass as -> flux_converter(flux, from, to, E)")

    # FROM photons/cm^2/s/A
    if(frm == FLUX_UNIT_OF_MEASURE[2]):
        # TO Jy
        if(to == FLUX_UNIT_OF_MEASURE[0]):
            if(E_LL != []):
                return 6.63 * 10**-4 * E_LL * flux
            else:
                raise Exception(
                    "ERROR: FROM photons/cm^2/s/A TO Jy -> Wavelength in A needed, pass as -> flux_converter(flux, from, to, LL)")
        # TO photons/cm^2/s/keV
        if(to == FLUX_UNIT_OF_MEASURE[1]):
            if(E_LL != []):
                return 8.07 * 10**-2 * E_LL**2 * flux
            else:
                raise Exception(
                    "ERROR: FROM photons/cm^2/s/A TO photons/cm^2/s/keV -> Wavelength in A needed, pass as -> flux_converter(flux, from, to, LL)")
        # TO ergs/cm^2/s/A
        if(to == FLUX_UNIT_OF_MEASURE[3]):
            if(E_LL != []):
                return 1.99 * 10**-8 * flux / E_LL
            else:
                raise Exception(
                    "ERROR: FROM photons/cm^2/s/A TO ergs/cm^2/s/A -> Wavelength in A needed, pass as -> flux_converter(flux, from, to, LL)")
        # TO ergs/cm^2/s/Hz
        if(to == FLUX_UNIT_OF_MEASURE[4]):
            if(E_LL != []):
                return 6.63 * 10**-27 * E_LL * flux
            else:
                raise Exception(
                    "ERROR: FROM photons/cm^2/s/A TO ergs/cm^2/s/Hz -> Wavelength in A needed, pass as -> flux_converter(flux, from, to, LL)")

    # FROM ergs/cm^2/s/A
    if(frm == FLUX_UNIT_OF_MEASURE[3]):
        # TO Jy
        if(to == FLUX_UNIT_OF_MEASURE[0]):
            if(E_LL != []):
                return 3.34 * 10**4 * E_LL**2 * flux
            else:
                raise Exception(
                    "ERROR: FROM ergs/cm^2/s/A TO Jy -> Wavelength in A needed, pass as -> flux_converter(flux, from, to, LL)")
        # TO photons/cm^2/s/keV
        if(to == FLUX_UNIT_OF_MEASURE[1]):
            if(E_LL != []):
                return 4.06 * 10**6 * E_LL**3 * flux
            else:
                raise Exception(
                    "ERROR: FROM ergs/cm^2/s/A TO photons/cm^2/s/keV -> Wavelength in A needed, pass as -> flux_converter(flux, from, to, LL)")
        # TO photons/cm^2/s/A
        if(to == FLUX_UNIT_OF_MEASURE[2]):
            if(E_LL != []):
                return 5.03411250 * 10**7 * E_LL * flux
                # ph/s/cm/A = energia_osservata * wave_strumento_tot * 10 ** (-10) * 10 ** (-7) / constants.h / constants.c EQUALS
            else:
                raise Exception(
                    "ERROR: FROM ergs/cm^2/s/A TO photons/cm^2/s/A -> Wavelength in A needed, pass as -> flux_converter(flux, from, to, LL)")
        # TO ergs/cm^2/s/Hz
        if(to == FLUX_UNIT_OF_MEASURE[4]):
            if(E_LL != []):
                return 3.34 * 10**-19 * E_LL**2 * flux
            else:
                raise Exception(
                    "ERROR: FROM ergs/cm^2/s/A TO ergs/cm^2/s/Hz -> Wavelength in A needed, pass as -> flux_converter(flux, from, to, LL)")

    # FROM ergs/cm^2/s/Hz
    if(frm == FLUX_UNIT_OF_MEASURE[4]):
        # TO Jy
        if(to == FLUX_UNIT_OF_MEASURE[0]):
            return 10**23 * flux
        # TO photons/cm^2/s/keV
        if(to == FLUX_UNIT_OF_MEASURE[1]):
            if(E_LL != []):
                return 1.51 * 10**26 * flux / E_LL
            else:
                raise Exception(
                    "ERROR: FROM ergs/cm^2/s/Hz TO photons/cm^2/s/keV -> Energy in keV needed, pass as -> flux_converter(flux, from, to, LL)")
        # TO photons/cm^2/s/A
        if(to == FLUX_UNIT_OF_MEASURE[2]):
            if(E_LL != []):
                return 1.51 * 10**26 * flux / E_LL
            else:
                raise Exception(
                    "ERROR: FROM ergs/cm^2/s/Hz TO photons/cm^2/s/A -> Wavelength in A needed, pass as -> flux_converter(flux, from, to, LL)")
        # TO ergs/cm^2/s/A
        if(to == FLUX_UNIT_OF_MEASURE[3]):
            if(E_LL != []):
                return 3.00 * 10**18 * flux / E_LL**2
            else:
                raise Exception(
                    "ERROR: FROM ergs/cm^2/s/Hz TO ergs/cm^2/s/A -> Wavelength in A needed, pass as -> flux_converter(flux, from, to, LL)")

    # FROM W/m^2/um
    if(frm == FLUX_UNIT_OF_MEASURE[5]):
        # TO photons/cm^2/s/A
        if(to == FLUX_UNIT_OF_MEASURE[2]):
            if(E_LL != []):  # to A
                return ((5.03411250 * 10**14 * flux * (E_LL)) / 10000)
            else:
                raise Exception(
                    "ERROR: FROM W/m^2/um TO photons/cm^2/s/A -> Wavelength in um needed, pass as -> flux_converter(flux, from, to, LL) / LL in um!!")

    # FROM J/(s * m^3)
    if(frm == FLUX_UNIT_OF_MEASURE[6]):
        # TO photons/cm^2/s/A
        if(to == FLUX_UNIT_OF_MEASURE[2]):
            if(E_LL != []):  # to A
                # [phot/s/cm2/A]
                return flux * E_LL * 10 ** (-10) * 10 ** (-7) / constants.h / constants.c
            else:
                raise Exception(
                    "ERROR: FROM W/m^2/um TO photons/cm^2/s/A -> Wavelength in um needed, pass as -> flux_converter(flux, from, to, LL) / LL in um!!")


def wavelength(wavelength, frm, to):
    if(frm not in WAVELENGTH_UNIT_OF_MEASURE):
        raise Exception(
            "ERROR: FROM - Unit of measure NOT present -> choose from " + str(WAVELENGTH_UNIT_OF_MEASURE))

    if(to not in WAVELENGTH_UNIT_OF_MEASURE):
        raise Exception(
            "ERROR: TO - Unit of measure NOT present -> choose from " + str(WAVELENGTH_UNIT_OF_MEASURE))

    if(frm == to):
        return flux

    # FROM um
    if(frm == WAVELENGTH_UNIT_OF_MEASURE[0]):
        # TO nm
        if(to == WAVELENGTH_UNIT_OF_MEASURE[1]):
            return wavelength * 1000
        # TO A
        if(to == WAVELENGTH_UNIT_OF_MEASURE[2]):
            return wavelength * 10000

    # FROM nm
    if(frm == WAVELENGTH_UNIT_OF_MEASURE[1]):
        # TO um
        if(to == WAVELENGTH_UNIT_OF_MEASURE[0]):
            return wavelength / 1000
        # TO A
        if(to == WAVELENGTH_UNIT_OF_MEASURE[2]):
            return wavelength * 10

    # FROM A
    if(frm == WAVELENGTH_UNIT_OF_MEASURE[2]):
        # TO um
        if(to == WAVELENGTH_UNIT_OF_MEASURE[0]):
            return wavelength / 10000
        # TO nm
        if(to == WAVELENGTH_UNIT_OF_MEASURE[1]):
            return wavelength / 10


# TEST
'''
print(flux_converter(3.6*10**-17, "ergs/cm^2/s/A", "photons/cm^2/s/A", 5550))  # 1.004 x 10^-5

print(flux_converter(3.631*10**-5, "Jy", "photons/cm^2/s/A", 5500)) #9.9 x 10^-6

print(flux_converter(3.631*10**-5, "Jy", "ergs/cm^2/s/Hz")) #3.6 x 10^-28

print(flux_converter(5.27451*10**-16, "W/m^2/um", "photons/cm^2/s/A", 4786/10000)) #9.9 x 10^-6

'''
