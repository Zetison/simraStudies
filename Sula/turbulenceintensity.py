#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  2 15:20:44 2019

@author: midjiyawaz

Modified by Jon Vegard Venås at SINTEF Digital
"""
from __future__ import division
import netCDF4
import numpy as np
import pandas as pd
from datetime import datetime
from dateutil.relativedelta import relativedelta
import requests
import time
import click
from os.path import expanduser

############## Perdelta function ############################################## 
def perdelta(start, end, delta):
    curr = start
    while curr < end:
        yield curr
        curr += delta
############## Mean dir calculation ###########################################
def calc_mean_wind_direction(wdir):
    """
    Routine calculates the mean wind direction for the dataset
    wdir.
    Is this the best way of doing this?
    """
    # Calculate mean wind direction
    d_rads = wdir*np.pi/180.
    mean_rads = np.arctan2(np.sin(d_rads).sum(), np.cos(d_rads).sum())
    mean_wdir = mean_rads*180./np.pi
    mean_wdir = mean_wdir % 360.
    # Stand and deliver
    return mean_wdir

def binarySearch(nc,T,n,units):
    L = 0
    R = n-1
    while L <= R:
        m = (L+R)//2
        A = netCDF4.num2date(nc.variables['time'][m],units,only_use_cftime_datetimes=False)
        if A < T:
            L = m + 1
        elif A > T:
            R = m - 1
        else:
            return m

def extractData(year=2020,month=11,day=1,hour=0,frequency='hz',time_interval=60,totPeriod=60,midSampled=False,location='Kvitneset'):
    if midSampled:
        start_dt = datetime(year, month, day, hour) + relativedelta(minutes=-30)
        dataType = 'rawMid'
    else:
        start_dt = datetime(year, month, day, hour)
        dataType = 'raw'

    end_dt = start_dt + relativedelta(minutes=+totPeriod) #40,320
    ############### Read data and clean ###########################################
    t0 = time.perf_counter()
    base_url = 'https://thredds.met.no/thredds/dodsC/obs/mast-svv-e39/%d/%d/%d%02d_%s_10%s.nc' % (year,month,year,month,location,frequency)
    nc_wind = netCDF4.Dataset(base_url)

    units = nc_wind.variables['time'].units
    sz = nc_wind.dimensions['time'].size
    t_coverage_start = datetime.strptime(nc_wind.time_coverage_start, '%Y-%m-%dT%H:%M:%S')
    index_start = binarySearch(nc_wind, start_dt, sz, units) 
    index_end   = binarySearch(nc_wind, end_dt, sz, units) 
    indices = slice(index_start,index_end+1)

    tAll = netCDF4.num2date(nc_wind.variables['time'][indices],units,only_use_cftime_datetimes=False)
    uAll = nc_wind.variables['windspeed'][indices,:]
    wAll = nc_wind.variables['vws'][indices,:]
    winddirAll = nc_wind.variables['winddirection'][indices,:]
    z = nc_wind.variables['alt'][:]

    print(time.perf_counter() - t0)
    t0 = time.perf_counter()
    for height_level in range(len(z)):
        t0 = time.perf_counter()
        df0 = pd.DataFrame({'Time':tAll,
                            'u':uAll[:,height_level],
                            'w':wAll[:,height_level],
                            'dir':winddirAll[:,height_level]})
        df0.dropna()
        df0['Time'] = pd.to_datetime(df0.Time)

        print(time.perf_counter() - t0)
        t0 = time.perf_counter()
        #################################### Program Start ############################
        print (location,z[height_level],start_dt)

        meanDir = []
        u_mag   = []
        mean_U  = []
        mean_u  = []
        mean_v  = []
        mean_w  = []
        alpha   = []
        dates   = []
        ###################### Gust wind speed ########################################
        for dt in perdelta(start_dt, end_dt,  relativedelta(minutes=+time_interval)):
            try:
                start = dt.strftime('%Y-%m-%d %H:%M')
                end = (dt +  relativedelta(minutes=+time_interval)).strftime('%Y-%m-%d %H:%M')
                print (start, end)
                date_rng10m = pd.date_range(start, end, freq='100ms')
                df10 = df0[np.in1d(df0['Time'], date_rng10m)]
                ########## Wind speed decomposition #######################################
                meandir = calc_mean_wind_direction(df10['dir'].dropna())
                theta = df10['dir'].dropna() - meandir
                Cosines = np.cos(np.deg2rad(theta))
                Sines = np.sin(np.deg2rad(theta))
                u = [a*b for a,b in zip(df10['u'].dropna(), Cosines)]
                v = [a*b for a,b in zip(df10['u'].dropna(), Sines)]
                meanu = np.mean(u)
                meanv = np.mean(v)
                ########## Statistical quatities ##########################################
                meanU = np.mean(df10['u'].dropna())
                meanw = np.mean(df10['w'].dropna())
                ########## Save everything in a array #####################################
                meanDir.append(meandir)
                u_mag.append(np.sqrt(meanu**2+meanv**2+meanw**2))
                mean_U.append(meanU)
                mean_u.append(meanu)
                mean_v.append(meanv)
                mean_w.append(meanw)
                alpha.append(np.degrees(np.arctan2(meanw,meanU)))
                if midSampled:
                    dates.append(np.array(np.mean(df10['Time'].dropna()), dtype='datetime64[m]'))
                else:
                    dates.append(np.array(df10['Time'].to_numpy()[0], dtype='datetime64[m]'))
            except Exception:
                   print("Warning - Bad data")

        ################### Write all result into a panda dataframe ####################
        df_turbstat = pd.DataFrame()
        df_turbstat['meandir'] = meanDir
        df_turbstat['date'] = dates

        df_turbstat['u_mag'] = u_mag
        df_turbstat['meanU'] = mean_U
        df_turbstat['meanu'] = mean_u
        df_turbstat['meanv'] = mean_v
        df_turbstat['meanw'] = mean_w

        df_turbstat['alpha'] = alpha

        home = expanduser("~")
        df_turbstat.to_csv(home+'/results/simra/Sula/measurements/'+dataType+'/10%s_%s_60mnturbulence_statistics_%d_%d%02d.csv' % (frequency,location, z[height_level], year, month), index=False)
        print(time.perf_counter() - t0)

@click.command()
def main():
    midSampled = False
    totPeriod = 30*60*24 - 60  # in minutes
    #totPeriod = 60  # in minutes

    for location in ['Kvitneset','Traelboneset','Langeneset','Kaarsteinen']:
        for midSampled in [False,True]:
            extractData(totPeriod=totPeriod,midSampled=midSampled,location=location)


if __name__ == '__main__':
    main()
