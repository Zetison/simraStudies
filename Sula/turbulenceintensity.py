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
import click
from os.path import expanduser

############## Perdelta function ############################################## 
def nearest_ind(items, pivot):
    time_diff = np.abs([date - pivot for date in items])
    boolArr = np.zeros((len(items)), dtype=bool)
    boolArr[time_diff.argmin(0)] = True
    return boolArr

def perdelta(start, end, delta):
    curr = start
    while curr < end:
        yield curr
        curr += delta

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
        dataType = 'rawMidNew'
    else:
        start_dt = datetime(year, month, day, hour)
        dataType = 'rawNew'

    end_dt = start_dt + relativedelta(minutes=+totPeriod) #40,320
    ############### Read data and clean ###########################################
    if location == 'Bridgecenter':
        base_url = '/home/zetison/results/simra/Sula/measurements/%d%02d_%s_1%s.nc' % (year,month,location,frequency)
    else:
        base_url = 'https://thredds.met.no/thredds/dodsC/obs/mast-svv-e39/%d/%d/%d%02d_%s_10%s.nc' % (year,month,year,month,location,frequency)

    nc_wind = netCDF4.Dataset(base_url)
    units = nc_wind.variables['time'].units
    if location == 'Bridgecenter':
        indices = slice(len(nc_wind.variables['time']))
    else:
        sz = nc_wind.dimensions['time'].size
        t_coverage_start = datetime.strptime(nc_wind.time_coverage_start, '%Y-%m-%dT%H:%M:%S')
        index_start = binarySearch(nc_wind, start_dt, sz, units) 
        index_end   = binarySearch(nc_wind, end_dt, sz, units) 
        indices = slice(index_start,index_end+1)

    tAll = netCDF4.num2date(nc_wind.variables['time'][indices],units,only_use_cftime_datetimes=False)
    uAll = nc_wind.variables['windspeed'][indices,:]
    if location == 'Bridgecenter':
        wAll = uAll * 0 
    else:
        wAll = nc_wind.variables['vws'][indices,:]
    winddirAll = nc_wind.variables['winddirection'][indices,:]
    z = nc_wind.variables['alt'][:]

    for height_level in range(len(z)):
        df0 = pd.DataFrame({'Time':tAll,
                            'windspeed':uAll[:,height_level],
                            'w':wAll[:,height_level],
                            'winddirection':winddirAll[:,height_level]})
        df0.dropna()
        df0['Time'] = pd.to_datetime(df0.Time)
        df0['u'] = -df0['windspeed']*np.sin(np.radians(df0['winddirection']))
        df0['v'] = -df0['windspeed']*np.cos(np.radians(df0['winddirection']))

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
            start = dt.strftime('%Y-%m-%d %H:%M')
            end = (dt + relativedelta(minutes=+time_interval)).strftime('%Y-%m-%d %H:%M')
            print (start, end)

            indices = np.logical_and(start <= df0['Time'], df0['Time'] <= end)
            if np.any(indices):
                dfUniform = pd.DataFrame()
                uniformTime = pd.date_range(start, end, freq='100ms')
                dfUniform['Time'] = uniformTime
                for field in ['windspeed','u','v','w']:
                    dfUniform[field] = np.interp(uniformTime,df0['Time'][indices],df0[field][indices])
            else:
                print('Warning - Bad data (no overlap in time)')
                continue

            ########## Wind speed decomposition #######################################
            meanu = np.mean(dfUniform['u'])
            meanv = np.mean(dfUniform['v'])
            meandir = (180+90-np.degrees(np.arctan2(meanv,meanu))) % 360.
            meanU = np.mean(dfUniform['windspeed'])
            meanw = np.mean(dfUniform['w'])

            ########## Save everything in a array #####################################
            meanDir.append(meandir)
            u_mag.append(np.sqrt(meanu**2+meanv**2+meanw**2))
            mean_U.append(meanU)
            mean_u.append(meanu)
            mean_v.append(meanv)
            mean_w.append(meanw)
            alpha.append(np.degrees(np.arctan2(meanw,meanU)))
            if midSampled:
                dates.append(np.array(np.mean(dfUniform['Time']).round(str(time_interval)+'min'), dtype='datetime64[m]'))
            else:
                dates.append(np.array(dfUniform['Time'].round(str(time_interval)+'min').to_numpy()[0], dtype='datetime64[m]'))

        ################### Write all result into a panda dataframe ####################
        df_turbstat = pd.DataFrame()
        df_turbstat['date'] = dates
        df_turbstat['meandir'] = meanDir

        
        if location == 'Bridgecenter':
            df_turbstat['u_mag'] = np.nan
            df_turbstat['meanw'] = np.nan
            df_turbstat['alpha'] = np.nan
        else:
            df_turbstat['u_mag'] = u_mag
            df_turbstat['meanw'] = mean_w
            df_turbstat['alpha'] = alpha

        df_turbstat['meanU'] = mean_U
        df_turbstat['meanu'] = mean_u
        df_turbstat['meanv'] = mean_v

        home = expanduser("~")
        filename = home+'/results/simra/Sula/measurements/'+dataType+'/10%s_%s_60mnturbulence_statistics_%d_%d%02d.csv' % (frequency,location, z[height_level], year, month)
        df_turbstat.to_csv(filename, index=False)
        print("File %s written." % filename)

@click.command()
def main():
    totPeriod = 30*60*24 - 60  # in minutes
    #totPeriod = 365*60*24 - 60  # in minutes
    #totPeriod = 60  # in minutes
    locations = ['Kvitneset','Traelboneset','Langeneset','Kaarsteinen','Bridgecenter']
    locations = ['Bridgecenter']
    locations = ['Kvitneset']
    for location in locations:
        for midSampled in [True]:
            #extractData(totPeriod=totPeriod,midSampled=midSampled,location=location,day=19,hour=15)
            extractData(totPeriod=totPeriod,midSampled=midSampled,location=location,month=11)


if __name__ == '__main__':
    main()
