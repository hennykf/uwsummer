from scipy.io.idl import readsav
from netCDF4 import Dataset, num2date
from mpl_toolkits.basemap import Basemap
import numpy as np
import healpy as hp
import matplotlib
import matplotlib.pyplot as plt
import os
import warnings


def plothealpix(filename, frequency, plottype='residual_data', plotfile_base=None, plotdirectory=None):
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    m = Basemap(projection='moll', lon_0=90, resolution='c')

    # draw parallels and meridians. Degrees&direction intrinsic to basemap
    m.drawparallels(np.arange(-90., 120., 30.))
    m.drawmeridians(np.arange(0., 420., 60.))

    contents = readsav(filename, python_dict=True)
    nside = int(contents['nside'])
    pixels = contents['hpx_inds']

    ra, dec = hp.pixelfunc.pix2ang(nside, pixels, nest=False, lonlat=True)
    ra[np.where(ra > 180)] -= 360
    contents.keys()

    # plot types available
    residual_data = contents['dirty_cube'] - contents['model_cube']
    dirty_data = contents['dirty_cube']
    model_data = contents['model_cube']
    if plotfile_base is None:
        plotfile_base = os.path.splitext(filename)[0]
    # For 'all' frequencies, creates folder containing individual plots
    # Includes both filename and specified plot name
    if frequency == 'all':
        if plotdirectory is None:
            plotdirectory = os.getcwd() + '/' + os.path.splitext(plotfile_base)[0] + '_all'
    # makes directory if specified, adds file f and data type info to filename
    if plotdirectory is None:
        plotdirectory = os.getcwd()
    if not os.path.isdir(plotdirectory):
        os.makedirs(plotdirectory)
    plotfile_base = os.path.join(plotdirectory, plotfile_base)
    plotfile = plotfile_base + '_f' + str(frequency) + '_' + str(plottype) + '.png'

    # average of all frequencies in cube
    if frequency is 'average':
        dirty_data = np.mean(dirty_data, axis=0)
        model_data = np.mean(model_data, axis=0)
        residual_data = np.mean(residual_data, axis=0)

    # all frequencies in cube to individual files:
    elif frequency is 'all':
        for i in range(0, len(residual_data[:, 0])):
            plothealpix(filename, i, plottype=plottype, plotfile_base=plotfile_base)

    elif isinstance(frequency, (int, long)) and frequency >= 0 and frequency <= 192:
        residual_data = residual_data[frequency, :]
        dirty_data = dirty_data[frequency, :]
        model_data = model_data[frequency, :]

    else:
        raise ValueError("frequency must be 'average', 'all', or a positive integer")

    # plot of Galactic gas with coordinate projection
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    # llc and urc are x & y ranges.
    m = Basemap(projection='hammer', llcrnrlon=-11, llcrnrlat=-15, urcrnrlon=13.5, urcrnrlat=-37, resolution='h', epsg=5520)
    x, y = m(ra, dec)
    # draw parallels and meridians. Labels are 1/0 as [Top,bottom,right,left]
    m.drawparallels(np.arange(-90., 120., 2.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(0., 420., 2.), labels=[0, 0, 0, 1])

    if plottype is 'dirty_data':
        data = dirty_data
    elif plottype is 'model_data':
        data = model_data
    elif plottype is 'residual_data':
        data = residual_data
    else:
        raise ValueError("plottype must be dirty_data, model_data, or residual_data")

    # plot labeling
    m.scatter(x, y, 3, marker='o', linewidths=.1, c=data, cmap=plt.cm.coolwarm)
    m.colorbar()
    plt.title('Galactic gas f=' + str(frequency))

    # either plt.show or plt.savefig
    # plt.show()
    plt.savefig(plotfile)
    print "saved plot to " + plotfile
