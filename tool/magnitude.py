# coding=utf-8

MAG_SYSTEMS = ["Vega", "AB"]

filters_vega = ["U", "B", "V", "R", "I", "J", "H", "K"]
ab_vega_diff_vf = [0.79, -.09, .02, .21, .45, .91, 1.39, 1.85]

filters_AB = ["u", "g", "r", "i", "z"]
ab_vega_diff_abf = [0.91, -.08, .16, .37, .54]


# SOURCE http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=Paranal/FORS1.ESO1033&&mode=browse&gname=Paranal&gname2=FORS1#filter
#                       BP: [0: ll in A, 1: Fv in [erg/cm^2/s/Hz], 2: Fll in [erg/cm^2/s/A], 3: PHll in [Photons/cm^2/s/A] ]
VEGA_flux_zeropoints_svo = {"U": [3596.21, None, 3.50079 * 10**-9, 633.8],  # NON HA SENSO è IL RISULTATO DELL'INTEGRALE
                            "B": [4242.96, None, 6.42193 * 10**-9, 1371.7],
                            "V": [5457.39, None, 3.58159 * 10**-9, 983.9],
                            "R": [6474.93, None, 2.10171 * 10**-9, 685.],
                            "I": [7864.59, None, 1.18093 * 10**-9, 467.5],
                            }


# SOURCE https://www.eso.org/observing/etc/doc/skycalc/helpskycalc.html#mags
#                       BP: [0: ll in A, 1: Fv in [erg/cm^2/s/Hz], 2: Fll in [erg/cm^2/s/A], 3: PHll in [Photons/cm^2/s/A] ]
VEGA_flux_zeropoints = {"U": [3600., None, 4.18023 * 10**-9, 757.5],
                        "B": [4380., None, 6.60085 * 10**-9, 1455.4],
                        "V": [5450., None, 3.60994 * 10**-9, 990.4],
                        "R": [6410., None, 2.28665 * 10**-9, 737.9],
                        "I": [7980., None, 1.22603 * 10**-9, 492.5],
                        "J": [12200., None, 3.12 * 10**-10, 191.6],
                        "H": [16300., None, 1.14 * 10**-10, 93.5],
                        "K": [21900., None, 3.94 * 10**-11, 43.4]
                        }


# SOURCE: http://www.astronomy.ohio-state.edu/~martini/usefuldata.html
#                       BP: [0: ll in A, 1: Fv in [erg/cm^2/s/Hz], 2: Fll in [erg/cm^2/s/A], 3: PHll in [Photons/cm^2/s/A] ]
VEGA_flux_zeropoints_ua = {"U": [3600., 1.79 * 10**-20, 417.5 * 10**-11, 756.1],
                           "B": [4380., 4.063 * 10**-20, 632 * 10**-11, 1392.6],
                           "V": [5450., 3.636 * 10**-20, 363.1 * 10**-11, 995.5],
                           "R": [6410., 3.064 * 10**-20, 217.7 * 10**-11, 702.],
                           "I": [7980., 2.416 * 10**-20, 112.6 * 10**-11, 452.],
                           "J": [12200., 1.589 * 10**-20, 31.47 * 10**-11, 193.1],
                           "H": [16300., 1.021 * 10**-20, 11.38 * 10**-11, 93.3],
                           "K": [21900., .64 * 10**-20, 3.961 * 10**-11, 43.6]
                           }

#                       BP: [0: ll in A, 1: Fv in Jy, 2: Fll in [erg/cm^2/s/A], 3: PHll in [Photons/cm^2/s/A] ]
AB_flux_zeropoints = {"u": [3560., 3631, 859.5 * 10**-11, 1539.3],
                      "g": [4830., 3631, 466.9 * 10**-11, 1134.6],
                      "r": [6260., 3631, 278. * 10**-11, 875.4],
                      "i": [7670., 3631, 185.2 * 10**-11, 714.5],
                      "z": [9100., 3631, 131.5 * 10**-11, 602.2]
                      }


def get_vega_flux_zeropoints(bandm, quantity):
    if(bandm not in filters_vega):
        raise Exception(
            "ERROR: Filter not present -> choose from [U,B,V,R,I,J,H,K]")

    if(quantity not in ["Fv", "Fll", "PHll"]):
        raise Exception(
            "ERROR: Filter not present -> choose from [Fv, Fll, PHll]")

    if(quantity == "Fv"):
        sel = 1
    elif(quantity == "Fll"):
        sel = 2
    elif(quantity == "PHll"):
        sel = 3

    #       lambda, zeropoint
    return VEGA_flux_zeropoints[bandm][0], VEGA_flux_zeropoints[bandm][sel]


def get_AB_flux_zeropoints(bandm, quantity):
    if(bandm not in filters_AB):
        raise Exception("ERROR: Filter not present -> choose from [u,g,r,i,z]")

    if(quantity not in ["Fv", "Fll", "PHll"]):
        raise Exception(
            "ERROR: Filter not present -> choose from [Fv, Fll, PHll]")

    if(quantity == "Fv"):
        sel = 1
    elif(quantity == "Fll"):
        sel = 2
    elif(quantity == "PHll"):
        sel = 3

    #       lambda, zeropoint
    return AB_flux_zeropoints[bandm][0], AB_flux_zeropoints[bandm][sel]

# print(get_vega_flux_zeropoints("U","PHll"))
# print(get_AB_flux_zeropoints("z","PHll"))

# ***********************************************************************************


# MAG CONVERTER

def vega_to_ab_converter(mag, bandm):

    if(bandm not in filters_AB):
        raise Exception(
            "ERROR: Bandpass NOT present -> choose from: " + ', '.join(filters_AB))

    return mag + ab_vega_diff_abf[filters_AB.index(bandm)]  # Mab = Mvega + MAG
    ####


def ab_to_vega_converter(mag, bandm):

    if(bandm not in filters_vega):
        raise Exception(
            "ERROR: Bandpass NOT present -> choose from: " + ', '.join(filters_vega))

    # Mvega = Mab - MAG
    return mag - ab_vega_diff_vf[filters_vega.index(bandm)]
