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
xLabels = ['Mean wind speed $[m/s]$', 'Mean wind direction $[^\circ]$', 'Mean angle of attack $[^\circ]$']
xArrayNames = ['u_mag', 'meandir', 'alpha']
#winddirWindow = np.array([0.,360.])
#winddirWindow = np.array([285., 345])
winddirWindow = np.array([292.5, 337.5])
xArrayLims = [[0,30],[0,360],[-30,30]]
noPlots = len(layoutNames)
noDataTypes = len(dataTypes)
sensorLoc,CoordUTM32,mastb,masth,Sensorh = computeSensorLoc(originx=-200,originy=6899800)
resultFileName = simraResultsFolder+'sampledResults.csv'
dfSim = pd.read_csv(resultFileName)

# Find dates at which the wind speed at Kvitneset (92m) is in the interval [285,330]
filenameKvitneset = simraResultsFolder+'measurements/rawMidNew/10hz_Kvitneset_60mnturbulence_statistics_92_202011.csv'
df_Kvitneset = pd.read_csv(filenameKvitneset)
dfObs_nw = df_Kvitneset[(winddirWindow[0] < df_Kvitneset.meandir) & (df_Kvitneset.meandir < winddirWindow[1]) & (df_Kvitneset.u_mag > 12)].copy()
noMasts = len(mastNames)
totNoSensors = 0
for i in range(noMasts):
    totNoSensors += len(Sensorh[i])

corrArr = np.zeros((totNoSensors,noDataTypes,noPlots))
errArr = np.zeros((totNoSensors,noDataTypes,noPlots))
QoImax_sim = -np.Inf*np.ones(3)
QoImax_obs = -np.Inf*np.ones(3)
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

                maxQoI_sim = np.abs(np.max(QoI_sim))

                if QoImax_sim[i_l] < maxQoI_sim:
                    QoImax_sim[i_l] = maxQoI_sim

                maxQoI_obs = np.abs(np.max(QoI_obs))
                if QoImax_obs[i_l] < maxQoI_obs:
                    QoImax_obs[i_l] = maxQoI_obs

                if xArrayNames[i_l] == 'meandir':
                    shift = np.mean(winddirWindow) - 180
                    QoI_sim = ((QoI_sim - shift) % 360) + shift
                    QoI_obs = ((QoI_obs - shift) % 360) + shift
                else:
                    shift = 0


                corr = stats.pearsonr(QoI_obs, QoI_sim)
                #if xArrayNames[i_l] == 'meandir' or xArrayNames[i_l] == 'alpha':
                #    err = np.linalg.norm(QoI_obs - QoI_sim)
                #else:
                #    err = 100*np.linalg.norm(QoI_obs - QoI_sim)/np.linalg.norm(QoI_obs)
                err = 100*np.linalg.norm(QoI_obs - QoI_sim)/np.linalg.norm(QoI_obs)
                corrArr[i_s][j][i_l] = corr[0]
                errArr[i_s][j][i_l] = err
                sc = plt.scatter(QoI_obs, QoI_sim,c=dfObs['meandir'],edgecolors='black',s=None, cmap=twilight)
                plt.colorbar(sc,label="Observation"+" "+"wind"+" "+"direction"+" "+"($^{\circ}$)")
                plt.title(mastNames[i] +", z = "+ str(z) +"m: Corr. "+"{0:.2f}".format(corr[0]) +", relative error "+"{0:.2f}%".format(err))
                xlim = np.array(xArrayLims[i_l]) + shift
                plt.plot(xlim, xlim,'k-',linewidth=2)
                plt.xlabel('Obs. '+xLabels[i_l],fontsize=14)
                plt.ylabel('Sim. '+xLabels[i_l],fontsize=14)
                xticks = np.linspace(xlim[0],xlim[1],9)
                if xArrayNames[i_l] == 'meandir':
                    plt.xticks(xticks, (xticks % 360).astype(int))
                    plt.yticks(xticks, (xticks % 360).astype(int))
                plt.xlim(xlim)
                plt.ylim(xlim)
                plt.clim((0,360))
                plt.grid(True)
                plt.tight_layout()
                caseName = layoutNames[i_l]+'_'+mastNames[i]+"_"+str(z)+"_"+dataTypes[j]+"_%d%02d"  % (year,month)
                print(caseName)
                plt.savefig(simraResultsFolder+'scatterPlots/'+caseName+'.pdf')
                plt.close()
                dfResult = dfSim_j[['date', xArrayNames[i_l]]].copy()
                dfResult[xArrayNames[i_l]+'_obs'] = dfObs[xArrayNames[i_l]]
                dfResult.to_csv(simraResultsFolder+'scatterPlots/'+caseName+'.csv',index=False)
                #plt.show()
        i_s += 1


with open(simraResultsFolder+'/scatterPlots/metadata.txt', 'w') as f:
    for i_l in range(noPlots):
        print('Max '+layoutNames[i_l]+' for simulations is '+str(QoImax_sim[i_l]), file=f)
        print('Max '+layoutNames[i_l]+' for observations is '+str(QoImax_obs[i_l]), file=f)

    i_arr = 0
    names = ['error','correlation']
    print('\n\n',file=f)
    for arr in [errArr, corrArr]:
        for j in range(noDataTypes):
            meanArr = np.mean(arr[:,j,:])
            print('Average '+names[i_arr]+' for '+dataTypes[j]+' is '+str(meanArr), file=f)
        i_arr += 1

    i_arr = 0
    for arr in [errArr, corrArr]:
        print('\nTable of mean '+names[i_arr]+' results', file=f)
        print(('\n%20s'+'%20s'*noPlots) % tuple(('',)+tuple(layoutNames)), file=f)
        for j in range(noDataTypes):
            meanArr = np.mean(arr[:,j,:],axis=0)
            print(('%20s' + '%20f'*len(meanArr)) % tuple(tuple((dataTypes[j],))+tuple(meanArr)), file=f)
        i_arr += 1

    i_arr = 0
    for arr in [errArr, corrArr]:
        print('\nTable of mean '+names[i_arr]+' results', file=f)
        for j in range(noDataTypes):
            print(('\nResults for '+dataTypes[j]+'\n%20s'+'%20s'*noPlots) % tuple(('',)+tuple(layoutNames)), file=f)
            i_s = 0
            for i in range(noMasts):
                noSensors = len(Sensorh[i])
                for k in range(noSensors):
                    z = np.floor(Sensorh[i][k]).astype(int)
                    meanArr = arr[i_s,j,:]
                    location = mastNames[i]
                    print(('%20s' + '%20f'*len(meanArr)) % tuple(tuple((location+' ['+str(z)+'m]',))+tuple(meanArr)), file=f)
                    i_s += 1
        i_arr += 1

    print('\nDates satisfying the criterias meandir in [%f,%f] and u_mag > 21' % tuple(winddirWindow), file=f)
    for _, row in dfObs_nw.iterrows():
        datestr = pd.to_datetime(row['date']).strftime('%Y-%m-%d %H:%M')
        print(datestr, file=f)

