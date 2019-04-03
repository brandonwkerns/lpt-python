import matplotlib; matplotlib.use('agg')
import numpy as np
from context import lpt
import matplotlib.pylab as plt
import datetime as dt
import sys
import os
import os.path
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap


"""
This module file contains function for making LPT related plots.
Most of the plotting function take a Matplotlib axis object as input
So create Matlplotlib figures and axes before calling the functions.
"""

def cmap_map(function, cmap):
    """ Applies function (which should operate on vectors of shape 3: [r, g, b]), on colormap cmap.
    This routine will break any discontinuous points in a colormap.
    """
    cdict = cmap._segmentdata
    step_dict = {}
    # Firt get the list of points where the segments start or end
    for key in ('red', 'green', 'blue'):
        step_dict[key] = list(map(lambda x: x[0], cdict[key]))
    step_list = sum(step_dict.values(), [])
    step_list = np.array(list(set(step_list)))
    # Then compute the LUT, and apply the function to the LUT
    reduced_cmap = lambda step : np.array(cmap(step)[0:3])
    old_LUT = np.array(list(map(reduced_cmap, step_list)))
    new_LUT = np.array(list(map(function, old_LUT)))
    # Now try to make a minimal segment definition of the new LUT
    cdict = {}
    for i, key in enumerate(['red','green','blue']):
        this_cdict = {}
        for j, step in enumerate(step_list):
            if step in step_dict[key]:
                this_cdict[step] = new_LUT[j, i]
            elif new_LUT[j,i] != old_LUT[j, i]:
                this_cdict[step] = new_LUT[j, i]
        colorvector = list(map(lambda x: x + (x[1], ), this_cdict.items()))
        colorvector.sort()
        cdict[key] = colorvector

    return colors.LinearSegmentedColormap('colormap',cdict,1024)


def plot_map_background(plotArea=[0,360,-60,60], lon_labels = [1,0,0,0], lat_labels = [0,0,0,1], res='c', anchor='C'
                        , coast_color = 'k', fontsize=10):
    map = Basemap(projection='cyl',resolution=res, anchor=anchor,
                  llcrnrlat = plotArea[2],
                  urcrnrlat = plotArea[3],
                  llcrnrlon = plotArea[0],
                  urcrnrlon = plotArea[1]);

    #map.fillcontinents(color='lightgray', lake_color='lightgray')
    map.drawcoastlines(linewidth=0.5, color=coast_color) ;
    map.drawparallels(np.arange(-80,81,20), linewidth=0.5, labels = lon_labels, fontsize = fontsize);
    map.drawmeridians(np.arange(0,361,20), linewidth=0.5, labels = lat_labels, fontsize = fontsize);

    return map ;


def print_and_save(file_out_base):
    os.makedirs(os.path.dirname(file_out_base), exist_ok=True) # Make directory if needed.
    print(file_out_base + '.png')
    plt.savefig(file_out_base + '.png' ,bbox_inches='tight', dpi=150)


"""
####################################################################
## High level plotting functions follow below.
####################################################################
"""

def plot_rain_map_with_filtered_contour(ax, DATA_ACCUM, OBJ, plot_area=[50, 200, -30, 30]):
    lon = OBJ['grid']['lon']
    lat = OBJ['grid']['lat']

    map1 = plot_map_background(plot_area)
    cmap = cmap_map(lambda x: x/2 + 0.5, plt.cm.jet)
    cmap.set_under(color='white')
    H1 = map1.pcolormesh(lon, lat, DATA_ACCUM, cmap=cmap, vmin=1, vmax=50)

    label_im = np.array(OBJ['label_im'])
    label_im[label_im > 0.5] = 1
    Hobj = plt.contour(lon, lat, label_im, [0.5,], colors='k', linewidths=1.0)

    map1.plot(OBJ['lon'], OBJ['lat'], 'kx', markersize=7)
    CB = plt.colorbar(H1)
    return (map1, H1, Hobj, CB)


def plot_rain_map_with_lpt_history_and_contour(ax):
    pass

def plot_rain_time_longitude(ax):
    pass

def plot_lpt_systems_time_longitude(ax):
    pass
