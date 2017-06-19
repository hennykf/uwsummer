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


def mapping(nside, pixelnum, plotfile, data, nest_or_ring='ring'):
    # file pixels may be organized as ring or nest
    if nest_or_ring is 'ring':
        ra, dec = hp.pixelfunc.pix2ang(int(nside), pixelnum, nest=False, lonlat=True)
    if nest_or_ring is 'nest':
        ra, dec = hp.pixelfunc.pix2ang(int(nside), pixelnum, nest=True, lonlat=True)
    ra[np.where(ra > 180)] -= 360
    # plot of Galactic gas with coordinate projection
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    # llc and urc are x & y ranges, and are specific to a location.
    # The latitude and longitude settings are part of basemap.
    m = Basemap(projection='hammer', llcrnrlon=-11, llcrnrlat=-15, urcrnrlon=13.5, urcrnrlat=-37, resolution='h', epsg=5520)
    x, y = m(ra, dec)
    # draw parallels and meridians. Labels are 1/0 as [Top,bottom,right,left]
    m.drawparallels(np.arange(-90., 120., 2.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(0., 420., 2.), labels=[0, 0, 0, 1])
    # creates a scatter plot of the selected data on a globe.
    m.scatter(x, y, 3, marker='o', linewidths=.1, c=data, cmap=plt.cm.coolwarm)
    m.colorbar()
    plt.show(plotfile)
