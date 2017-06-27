from scipy.io.idl import readsav
from netCDF4 import Dataset, num2date
from mpl_toolkits.basemap import Basemap
from astropy.io import fits
import numpy as np
import healpy as hp
import scipy as sp
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os
import plothealpix_map_coords
import plothealpix_map


def plotstokesQ_U(filename_Q, filename_U, plotfile_base=None, plotdirectory=None, save_show=None):
    if not filename_U.endswith('.fits'):
        raise ValueError("Stokes U file is not a fits file. Files must be .fits files")
    if not filename_Q.endswith('.fits'):
        raise ValueError("Stokes Q file is not a fits file. Files must be .fits files")
    # get header and data info from the fits file
    contents_Q = fits.open(filename_Q)
    pixelnum_Q = contents_Q[1].header['naxis1'] * contents_Q[1].header['naxis2']
    data_Q = contents_Q[1].data
    nside_Q = contents_Q[1].header['nside']

    contents_U = fits.open(filename_U)
    pixelnum_U = contents_U[1].header['naxis1'] * contents_U[1].header['naxis2']
    data_U = contents_U[1].data
    nside_U = contents_U[1].header['nside']

    if pixelnum_U == pixelnum_Q and nside_U == nside_Q:
        pixelnum_Q = pixelnum_U
        nside_Q = nside_U
    if not pixelnum_U == pixelnum_Q and nside_U == nside_Q:
        raise ValueError("files do not have same indices.")

    # extract data from specified files
    pixels_Q = data_Q.field('PIXEL')
    signal_Q = data_Q.field('SIGNAL')
    pixels_U = data_U.field('PIXEL')
    signal_U = data_U.field('SIGNAL')

    # Finding x and y from Stokes parameters U and Q
    Q = signal_Q
    U = signal_U
    K = np.sqrt(U**2 + Q**2)
    U_pos = U[np.where(U >= 0)]
    U_neg = U[np.where(U < 0)]
    # theta is in radians
    theta = (.5 * np.arccos(((K + Q) / K) - 1))
    theta[np.where(U >= 0)] = theta
    theta[np.where(U < 0)] = theta + np.pi / 2
    # the x and y components of the theta-mag points.
    x_stokes = K * np.cos(theta)
    y_stokes = K * np.sin(theta)

    if plotfile_base is None:
        plotfile_base = os.path.splitext(filename_Q)[0] + '_' + os.path.splitext(filename_U)[0]

    if plotdirectory is None:
        plotdirectory = os.getcwd()
    if not os.path.isdir(plotdirectory):
        os.makedirs(plotdirectory)
    plotfile_base = os.path.join(plotdirectory, plotfile_base)
    plotfile = plotfile_base + '.png'
    # manipulate histogram for number of bins, color scheme.
    plt.hist2d(x_stokes, y_stokes, bins=150, norm=LogNorm())
    plt.colorbar()
    # either show the graph, or save it to a location.
    if save_show is None or save_show == 'show':
        plt.show()
    elif save_show == 'save':
        plt.savefig(plotfile)
        print "saved plot to " + plotfile
    else:
        raise ValueError("Do you want to save or show the image?")
    # using the function defined in plothealpix_map to graph the data on a globe.
    plothealpix_map.mapping(nside_Q, pixels_Q, plotfile, theta, 'nest')
    plt.savefig()
    return x_stokes, y_stokes
