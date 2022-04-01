# coding=utf-8

import matplotlib.pyplot as plt


def generic(title, wavelength, flux, labels, legend=False):

    plt.title(title)

    plt.plot(wavelength, flux)

    plt.xlabel(labels[0])
    plt.ylabel(labels[1])

    if legend:
        plt.legend()

    plt.show()

# Plot the image of the detector


def detector(title, image):
    plt.figure()
    plt.title(title)
    plt.imshow(image)
    plt.colorbar()
    plt.show()
