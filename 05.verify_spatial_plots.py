#!/usr/bin/env python
##!/data/aqf2/barryb/anaconda2/envs/website/bin/python

###############################################################
# < next few lines under version control, D O  N O T  E D I T >
# $Date: 2018-03-29 10:12:00 -0400 (Thu, 29 Mar 2018) $
# $Revision: 100014 $
# $Author: Barry.Baker@noaa.gov $
# $Id: nemsio2nc4.py 100014 2018-03-29 14:12:00Z Barry.Baker@noaa.gov $
###############################################################

__author__ = "Patrick Campbell"
__email__ = "Patrick.C.Campbell@noaa.gov"
__license__ = "GPL"

import os
import subprocess
import sys

# sys.path.append('/data/aqf/patrickc/MONET/')
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

import dask
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import cartopy.crs as ccrs

mpl.use("agg")

import pandas as pd
import seaborn as sns
import numpy as np
import monet
from monet.util.tools import calc_8hr_rolling_max, calc_24hr_ave, get_relhum

sns.set_context("notebook")

plt.ioff()

"""
Simple utility to make spatial plots from the NAQFC forecast and overlay observations
"""

initial_datetime = None


def make_24hr_regulatory(df, col=None):
    """ Make 24-hour averages """
    return calc_24hr_ave(df, col)


def make_8hr_regulatory(df, col=None):
    """ Make 8-hour rolling average daily """
    return calc_8hr_rolling_max(df, col, window=8)


def map_projection(f):
    import cartopy.crs as ccrs

    proj = ccrs.LambertConformal(central_longitude=f.XCENT, central_latitude=f.YCENT)
    return proj


def chdir(fname):
    dir_path = os.path.dirname(os.path.realpath(fname))
    os.chdir(dir_path)
    return os.path.basename(fname)


def open_cmaq(finput):
    from monet.models import cmaq

    f = cmaq.open_mfdataset(finput)
    return f


def load_paired_data(fname):
    return pd.read_hdf(fname)


def make_spatial_plot(da, df, outname, proj, startdate, enddate, region="domain"):
    cbar_kwargs = dict(aspect=30, shrink=0.8)  # dict(aspect=30)

    if region == "domain":
        latmin = 25.0
        lonmin = -130.0
        latmax = 55.0
        lonmax = -55.0
    else:
        from monet.util.tools import get_epa_region_bounds as get_epa_bounds

        latmin, lonmin, latmax, lonmax, acro = get_epa_bounds(
            index=None, acronym=region
        )

    extent = [lonmin, lonmax, latmin, latmax]
    ax = da.monet.quick_map(
        cbar_kwargs=cbar_kwargs,
        figsize=(15, 8),
        map_kwarg={"states": True, "crs": proj, "extent": extent},
        robust=True,
    )
    plt.gcf().canvas.draw()
    plt.tight_layout(pad=-1)
    if startdate == None and enddate == None:
     date = pd.Timestamp(da.time.values) 
     dt = date - initial_datetime
     dtstr = str(dt.days * 24 + dt.seconds // 3600).zfill(2)
     plt.title(date.strftime('time=%Y/%m/%d %H:00 | CMAQ - AIRNOW '))
    else:
        plt.title("average time period | CMAQ - AIRNOW ")

    cbar = ax.figure.get_axes()[1]
    vmin, vmax = cbar.get_ybound()
    vars = df.keys()
    varname = [x for x in vars if x not in ["latitude", "longitude"]][0]
    ax.scatter(
        df.longitude.values,
        df.latitude.values,
        s=25,
        c=df[varname],
        transform=ccrs.PlateCarree(),
        edgecolor="w",
        linewidth=0.08,
        vmin=vmin,
        vmax=vmax,
    )
    ax.set_extent(extent, crs=ccrs.PlateCarree())

    if startdate == None and enddate == None:

     savename = "{}.{}.{}.jpg".format(outname,
                                     initial_datetime.strftime('sp'),
                                     dtstr)
    else:
        savename = "{}.{}.jpg".format(outname, "sp")
    print(savename)
    monet.plots.savefig(savename, bbox_inches="tight", dpi=100, decorate=True)
    plt.close()


def make_plots(
    finput,
    paired_data,
    variable,
    obs_variable,
    verbose,
    startdate,
    enddate,
    region,
    outname,
):
    if startdate == None and enddate == None:
        # open the files
        f = open_cmaq(finput)
        # get map projection
        proj = map_projection(f)
        if paired_data is not None:
            df = paired_data
        # loop over varaible list
        plots = []
        #    for index, var in enumerate(variable):
        obj = f[variable]
        # loop over time
        for t in obj.time:
            date = pd.Timestamp(t.values)
            print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            print("Creating Plot:", obs_variable, "at time:", date)
            print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            odf = df.loc[
                df.time == pd.Timestamp(t.values),
                ["latitude", "longitude", obs_variable],
            ]
            make_spatial_plot(
                obj.sel(time=t), odf, outname, proj, startdate, enddate, region=region
            )
    #            plots.append(dask.delayed(make_spatial_plot)
    #                         (obj.sel(time=t), odf, proj))
    #        plots dask.delayed(make_spatial_plot)(
    #            obj.sel(time=t), proj) for t in obj.time]
    #    dask.delayed(plots).compute()
    # if paired_data is not None:
    #     ov = obs_variable[[index, 'latitude', 'longitude'].loc[obs_variable.time == t]
    # ax.scatter()
    else:
        # open the files
        f = open_cmaq(finput)
        # get map projection
        proj = map_projection(f)
        if paired_data is not None:
            df = paired_data
        # loop over varaible list
        plots = []
        #    for index, var in enumerate(variable):
        obj = f[variable]
        sdate = pd.Timestamp(startdate)
        edate = pd.Timestamp(enddate)
        df_mean = df.groupby(["siteid"], as_index=False).mean()
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        print("Creating Plot:", obs_variable, "for period:", startdate, "to ", enddate)
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        odf_mean = df_mean[["latitude", "longitude", obs_variable]]
        mod_slice = obj.sel(time=slice(startdate, enddate))
        mod_mean = mod_slice.mean(dim="time")
        make_spatial_plot(
            mod_mean, odf_mean, outname, proj, startdate, enddate, region=region
        )


if __name__ == "__main__":

    parser = ArgumentParser(
        description="Make Spatial Plots for each time step in files",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-f",
        "--files",
        help="input model file names",
        nargs="+",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-p",
        "--paired_data",
        help="associated paired data input file name",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-s",
        "--species",
        nargs="+",
        help="species to plot",
        type=str,
        required=False,
        default={"OZONE"},
    )
    parser.add_argument(
        "-v",
        "--verbose",
        help="print debugging information",
        action="store_true",
        required=False,
    )
    parser.add_argument(
        "-b",
        "--subset_epa",
        help="EPA Region Subset true/false",
        type=bool,
        required=False,
        default=False,
    )
    parser.add_argument(
        "-e",
        "--epa_region",
        help="EPA Region Acronyms",
        type=str,
        required=False,
        default="domain",
    )
    parser.add_argument(
        "-n",
        "--output_name",
        help="Output base name",
        type=str,
        required=False,
        default="CMAQ_AIRNOW",
    )
    parser.add_argument(
        "-sup",
        "--suppress_xwindow",
        help="Suppress X Window",
        action="store_true",
        required=False,
    )
    parser.add_argument(
        "-r",
        "--regulatory",
        help="boolean set to True fore 8-hrmax  or 24-ave NAAQS regulatory calcs",
        type=bool,
        required=False,
        default=False,
    )
    parser.add_argument(
        "-sd",
        "--startdate",
        help="Startdate for bias plot statistics over a period YYYY-MM-DD HH:MM:SS",
        type=str,
        required=False,
        default=None,
    )
    parser.add_argument(
        "-ed",
        "--enddate",
        help="Enddate for bias plot statisics over a period YYYY-MM-DD HH:MM:SS",
        type=str,
        required=False,
        default=None,
    )

    args = parser.parse_args()

    finput = args.files
    verbose = args.verbose
    paired_data = args.paired_data
    species = args.species
    out_name = args.output_name
    subset = args.subset_epa
    region = args.epa_region
    startdate = args.startdate
    enddate = args.enddate
    reg = args.regulatory

    # load the paired dataframe

    df = load_paired_data(paired_data)
    mapping_table = {
        "OZONE": "O3",
        "PM2.5": "PM25",
        "PM10": "PM10",
        "CO": "CO",
        "NO": "NO",
        "NO2": "NO2",
        "SO2": "SO2",
        "NOX": "NOX",
        "NO2Y": "NOY",
        "TEMP": "TEMP2",
        "WS": "WSPD10",
        "WD": "WDIR10",
        "SRAD": "GSW",
        "BARPR": "PRSFC",
        "PRECIP": "RT",
        "RHUM": "Q2",
    }
    sub_map = {i: mapping_table[i] for i in species if i in mapping_table}
    if region is "domain":
        subset = False
    # subset  only the correct region
    if subset is True:
        df.query("epa_region == " + '"' + region + '"', inplace=True)

    # Loop through species
    for jj in species:

