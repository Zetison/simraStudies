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
#dataTypes = ['raw']
#mastNames = ['Kvitneset']
#layoutNames = ['VelocityProfiles']
#xArrayNames = ['u_mag']
#xArrayLims = [[0,25]]
dataTypes = ['raw','rawMid','rawNew','rawMidNew']
mastNames = ['Kvitneset','Traelboneset','Langeneset','Kaarsteinen']
layoutNames = ['VelocityProfiles', 'WindDirProfiles', 'alphaProfiles']
xArrayNames = ['u_mag', 'meandir', 'alpha']
xArrayLims = [[0,25],[285,330],[-25,25]]
noPlots = len(layoutNames)
noDataTypes = len(dataTypes)
sensorLoc,CoordUTM32,mastb,masth,Sensorh = computeSensorLoc(originx=-200,originy=6899800)
resultFileName = simraResultsFolder+'sampledResults.csv'
dfSim = pd.read_csv(resultFileName)

# Find dates at which the wind speed at Kvitneset (92m) is in the interval [285,330]
filenameKvitneset = simraResultsFolder+'measurements/rawMidNew/10hz_Kvitneset_60mnturbulence_statistics_92_202011.csv'
df_Kvitneset = pd.read_csv(filenameKvitneset)
dfObs_nw = df_Kvitneset[(285 < df_Kvitneset.meandir) & (df_Kvitneset.meandir < 330) & (df_Kvitneset.u_mag > 12)].copy()
noMasts = len(mastNames)
totNoSensors = 0
for i in range(noMasts):
    totNoSensors += len(Sensorh[i])

corrArr = np.zeros((totNoSensors,noDataTypes,noPlots))
i_s = 0
for i in range(noMasts):
    noSensors = len(Sensorh[i])
    for k in range(noSensors):
        z = np.floor(Sensorh[i][k]).astype(int)
        dfSim_i = dfSim[(dfSim.name == 'M1') & (dfSim.location == mastNames[i]) & (dfSim.z == sensorLoc[i][k,-1]) & (dfSim.addtime < 6)]
        for j in range(noDataTypes):
            # Filter data based on the criteria above
            filename = simraResultsFolder+'measurements/'+dataTypes[j]+'/10hz_'+mastNames[i]+'_60mnturbulence_statistics_'+str(z)+'_202011.csv'
            df_all = pd.read_csv(filename)
            indices_obs = np.isin(np.array(df_all.date,dtype='datetime64[m]'),np.intersect1d(np.array(dfObs_nw.date,dtype='datetime64[m]'),np.array(dfSim_i.date,dtype='datetime64[m]')))
            dfObs = df_all[indices_obs].copy()
            dfObs = dfObs.reset_index()
            indices_sim = np.isin(np.array(dfSim_i.date,dtype='datetime64[m]'),np.array(dfObs.date,dtype='datetime64[m]'))
            dfSim_j = dfSim_i[indices_sim].copy()
            dfSim_j = dfSim_j.reset_index()
            for i_l in range(noPlots):
                QoI_sim = dfSim_j[xArrayNames[i_l]]
                QoI_obs = dfObs[xArrayNames[i_l]]
                corr = stats.pearsonr(QoI_obs, QoI_sim)
                corrArr[i_s][j][i_l] = corr[0]
                sc = plt.scatter(QoI_obs, QoI_sim,c=dfObs['meandir'],edgecolors='black',s=None, cmap=twilight)
                plt.colorbar(sc,label="Observation"+" "+"wind"+" "+"direction"+" "+"($^{\circ}$)")
                plt.title(mastNames[i] +", z = "+ str(z) +": " + "Corr. "+"{0:.2f}".format(corr[0]))
                plt.plot(xArrayLims[i_l], xArrayLims[i_l],'k-',linewidth=2)
                plt.xlabel('Obs. Mean wind speed  $(ms^{-1})$',fontsize=14)
                plt.ylabel('Sim. Mean wind speed  $(ms^{-1})$',fontsize=14)
                plt.xlim(xArrayLims[i_l])
                plt.ylim(xArrayLims[i_l])
                plt.clim((0,360))
                plt.grid(True)
                plt.tight_layout()
                caseName = layoutNames[i_l]+'_'+mastNames[i]+"_"+str(z)+"_"+dataTypes[j]+"_%d%02d"  % (year,month)
                print(caseName)
                plt.savefig(simraResultsFolder+'scatterPlots/northwest/'+caseName+'.pdf')
                plt.close()
                dfResult = dfSim_j[['date', xArrayNames[i_l]]].copy()
                dfResult[xArrayNames[i_l]+'_obs'] = dfObs[xArrayNames[i_l]]
                dfResult.to_csv(simraResultsFolder+'scatterPlots/northwest/'+caseName+'.csv',index=False)
                #plt.show()
        i_s += 1


with open(simraResultsFolder+'/scatterPlots/northwest/metadata.txt', 'w') as f:
    for j in range(noDataTypes):
        meanCorrArr = np.mean(corrArr[:,j,:])
        print('Average corr for '+dataTypes[j]+' is '+str(meanCorrArr), file=f)

    print('\nTable of mean results', file=f)
    print(('\n%20s'+'%20s'*noPlots) % tuple(('',)+tuple(layoutNames)), file=f)
    for j in range(noDataTypes):
        meanCorrArr = np.mean(corrArr[:,j,:],axis=0)
        print(('%20s' + '%20f'*len(meanCorrArr)) % tuple(tuple((dataTypes[j],))+tuple(meanCorrArr)), file=f)

    for j in range(noDataTypes):
        print(('\nResults for '+dataTypes[j]+'\n%20s'+'%20s'*noPlots) % tuple(('',)+tuple(layoutNames)), file=f)
        i_s = 0
        for i in range(noMasts):
            noSensors = len(Sensorh[i])
            for k in range(noSensors):
                z = np.floor(Sensorh[i][k]).astype(int)
                meanCorrArr = corrArr[i_s,j,:]
                location = mastNames[i]
                print(('%20s' + '%20f'*len(meanCorrArr)) % tuple(tuple((location+' ['+str(z)+'m]',))+tuple(meanCorrArr)), file=f)
                i_s += 1

    print('\nDates satisfying the criterias meandir in [285,330] and u_mag > 21', file=f)
    for _, row in dfObs_nw.iterrows():
        datestr = pd.to_datetime(row['date']).strftime('%Y-%m-%d %H:%M')
        print(datestr, file=f)

