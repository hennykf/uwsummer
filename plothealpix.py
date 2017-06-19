from scipy.io.idl import readsav
from netCDF4 import Dataset, num2date
from mpl_toolkits.basemap import Basemap
from astropy.io import fits
import numpy as np
import healpy as hp
import scipy as sp
import matplotlib
import matplotlib.pyplot as plt
import os
import warnings
import plothealpix_map


def plothealpix(filename, frequency, plottype='residual_data', plotfile_base=None, plotdirectory=None, save_show=None):
    warnings.filterwarnings("ignore", category=DeprecationWarning)

    # deals with .sav file types
    if filename.endswith('.sav'):
        contents = readsav(filename, python_dict=True)
        pixels = contents['hpx_inds']
        nside = int(contents['nside'])
        contents.keys()

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

        # all frequencies in cube to individual files
        elif frequency is 'all':
            for i in range(0, len(residual_data[:, 0])):
                plothealpix(filename, i, plottype=plottype, plotfile_base=plotfile_base)
        # if a specific frequency is given
        elif isinstance(frequency, (int, long)) and frequency >= 0 and frequency <= 192:
            residual_data = residual_data[frequency, :]
            dirty_data = dirty_data[frequency, :]
            model_data = model_data[frequency, :]

        else:
            raise ValueError("frequency must be 'average', 'all', or a positive integer")
        if plottype is 'dirty_data':
            data = dirty_data
        elif plottype is 'model_data':
            data = model_data
        elif plottype is 'residual_data':
            data = residual_data
        else:
            raise ValueError("plottype must be dirty_data, model_data, or residual_data")
        print data.shape
        plothealpix_map.mapping(nside, pixels, plotfile, data, 'ring')

        # plot labeling, specific to .sav files
        plt.title('Galactic gas f=' + str(frequency))

    plt.savefig(plotfile)
    print "saved plot to " + plotfile
