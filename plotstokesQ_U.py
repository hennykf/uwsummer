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


def plotstokesQ_U(filename_Q, filename_U, plotfile_base=None, plotdirectory=None, save_show=None):
    if not filename_U.endswith('.fits'):
        raise ValueError("files must be .fits")
    if not filename_Q.endswith('.fits'):
        raise ValueError("files must be .fits")
    # get header and data info from the fits file
    contents_Q = fits.open(filename_Q)
    pixelnum_Q = contents_Q[1].header['naxis1'] * contents_Q[1].header['naxis2']
    data_Q = contents_Q[1].data
    nside_Q = contents_Q[1].header['nside']

    contents_U = fits.open(filename_U)
    pixelnum_U = contents_U[1].header['naxis1'] * contents_U[1].header['naxis2']
    data_U = contents_U[1].data
    nside_U = contents_U[1].header['nside']

    # extract data from specified files
    pixels_Q = data_Q.field('PIXEL')
    signal_Q = data_Q.field('SIGNAL')
    pixels_U = data_U.field('PIXEL')
    signal_U = data_U.field('SIGNAL')

    ra_Q, dec_Q = hp.pixelfunc.pix2ang(int(nside_Q), pixels_Q, nest=False, lonlat=True)
    ra_Q[np.where(ra_Q > 180)] -= 360
    ra_U, dec_U = hp.pixelfunc.pix2ang(int(nside_U), pixels_Q, nest=False, lonlat=True)
    ra_U[np.where(ra_U > 180)] -= 360

    # Finding x and y from Stokes parameters U and Q
    Q = signal_Q
    U = signal_U
    # print U
    K = np.sqrt(U**2 + Q**2)
    U_pos = U[np.where(U >= 0)]
    U_neg = U[np.where(U < 0)]

    # theta is in radians
    theta = (.5 * np.arccos(((K + Q) / K) - 1))
    theta[np.where(U >= 0)] = theta
    theta[np.where(U < 0)] = theta + np.pi / 2
    theta.shape

    x = K * np.cos(theta)
    y = K * np.sin(theta)

    if plotfile_base is None:
        plotfile_base = os.path.splitext(filename_Q)[0] + '_' + os.path.splitext(filename_U)[0]

    if plotdirectory is None:
        plotdirectory = os.getcwd()
    if not os.path.isdir(plotdirectory):
        os.makedirs(plotdirectory)
    plotfile_base = os.path.join(plotdirectory, plotfile_base)
    plotfile = plotfile_base + '.png'

    plt.hist2d(x, y, bins=150, norm=LogNorm())
    plt.colorbar()

    if save_show is None or save_show == 'show':
        plt.show()
    elif save_show == 'save':
        plt.savefig(plotfile)
        print "saved plot to " + plotfile
    else:
        raise ValueError("Do you want to save or show the image?")
