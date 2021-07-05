#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 08:33:08 2020

@author: midjiyawaz

Modified by Jon Vegard Venås at SINTEF Digital
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from convertCoords import computeSensorLoc 
from os.path import expanduser

home = expanduser("~")
year = 2020
month = 11

simraResultsFolder = home+'/results/simra/Sula/'
bwr=plt.get_cmap('bwr')
twilight=plt.get_cmap('twilight')

viridis=plt.get_cmap('viridis')

########################### Give parameters and file names ######################################
dataTypes = ['raw','rawMid','rawNew','rawMidNew']
mastNames = ['Kvitneset','Traelboneset','Langeneset','Kaarsteinen']
layoutNames = ['VelocityProfiles', 'WindDirProfiles', 'alphaProfiles']
xArrayNames = ['u_mag', 'meandir', 'alpha']
noPlots = len(layoutNames)
noDataTypes = len(dataTypes)
sensorLoc,CoordUTM32,mastb,masth,Sensorh = computeSensorLoc(originx=-200,originy=6899800)
resultFileName = simraResultsFolder+'sampledResults.csv'
dfSim = pd.read_csv(resultFileName)
for i in range(len(mastNames)):
    noSensors = len(Sensorh[i])
    for k in range(noSensors):
        z = np.floor(Sensorh[i][k]).astype(int)
        dfSim_i = dfSim[(dfSim.name == 'M1') & (dfSim.location == mastNames[i]) & (dfSim.z == sensorLoc[i][k,-1])]
        for j in range(noDataTypes):
            filename = simraResultsFolder+'measurements/'+dataTypes[j]+'/10hz_'+mastNames[i]+'_60mnturbulence_statistics_'+str(z)+'_202011.csv'
            df_all = pd.read_csv(filename)
            dfObs = df_all[np.isin(np.array(df_all.date,dtype='datetime64[m]'),np.array(dfSim_i.date,dtype='datetime64[m]'))]
            for i_l in range(noPlots):
                QoI_sim = dfSim_i[xArrayNames[i_l]]
                QoI_obs = dfObs[xArrayNames[i_l]]
                corr_uu = stats.pearsonr(QoI_obs, QoI_sim)
                sc = plt.scatter(QoI_obs, QoI_sim,c=dfObs['meandir'],edgecolors='black',s=None, cmap=twilight)
                plt.colorbar(sc,label="Observation"+" "+"wind"+" "+"direction"+" "+"($^{\circ}$)")
                plt.title(mastNames[i] +", z = "+ str(z) +": " + "Corr. "+"{0:.2f}".format(corr_uu[0]))
                plt.plot([0, 25], [0, 25],'k-',linewidth=2)
                plt.xlabel('Obs. Mean wind speed  $(ms^{-1})$',fontsize=14)
                plt.ylabel('Sim. Mean wind speed  $(ms^{-1})$',fontsize=14)
                plt.ylim((0,25))
                plt.xlim((0,25))
                plt.clim((0,360))
                plt.grid(True)
                plt.tight_layout()
                plt.savefig(simraResultsFolder+'scatterPlots/'+"Mean_wind"+"_"+mastNames[i]+"_"+str(z)+"_"+dataTypes[j]+"_%d%02d.pdf"  % (year,month))
                plt.close()
                #plt.show()
