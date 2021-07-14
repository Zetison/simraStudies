import numpy as np
import matplotlib.pylab as plt
from scipy import interpolate
import pandas as pd
import click
from vtk import vtkXMLUnstructuredGridReader, vtkXMLStructuredGridReader, vtkPoints, vtkStructuredGrid, vtkProbeFilter
from vtk import vtkXMLStructuredGridWriter
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy
from glob import glob
from os.path import expanduser, isfile
from convertCoords import computeSensorLoc 
from datetime import datetime, timedelta
from matplotlib import cm
from siso.coords import graph, Coords

def getQoI(name,utmPoints,uUTM,u_mag,useDeg=False):
    src = Coords.find('utm:33n')
    tgt = Coords.find('geodetic')
    converter = graph.path(src, tgt)
    out_points = converter.points(src, tgt, utmPoints, None)
    u = converter.vectors(src, tgt, uUTM, None)
    if name == 'VelocityProfiles':
        QoI = u_mag
    elif name == 'WindDirProfiles':
        if useDeg:
            QoI = (180+90-np.degrees(np.arctan2(u[:,1],u[:,0]))) % 360.
        else:
            QoI = (np.radians(180+90)-np.arctan2(u[:,1],u[:,0])) % (2*np.pi)
    elif name == 'alphaProfiles':
        QoI = np.degrees(np.arctan2(u[:,2],np.linalg.norm(u[:,:2],axis=1)))
    return QoI

def load_vtk(filename):
    if filename[-3:] == 'vts':
        reader = vtkXMLStructuredGridReader()
    elif filename[-3:] == 'vtu':
        reader = vtkXMLUnstructuredGridReader()
    else:
        print('Error: Not implemented')

    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()

def createVTKcurve(points):
    vpoints = vtkPoints()
    vpoints.SetData(numpy_to_vtk(points, deep=True))

    vgrid = vtkStructuredGrid()
    vgrid.SetDimensions(points.shape[0],1,1)
    vgrid.SetPoints(vpoints)

    return vgrid

def interpolatePoints(nodes,noPts=1000,z_min=0.0,z_max=120.0):
    nodes = np.vstack([np.append(nodes[0,:-1],z_max), nodes,np.append(nodes[-1,:-1],z_min)])
    num_true_pts = nodes.shape[0]

    tck, u = interpolate.splprep(np.transpose(nodes),k=num_true_pts-1)
    x_knots, y_knots, z_knots = interpolate.splev(tck[0], tck)
    u_fine = np.linspace(0,1,noPts)
    x_fine, y_fine, z_fine = interpolate.splev(u_fine, tck)

    points = np.transpose([x_fine, y_fine, z_fine])
    return createVTKcurve(points)

def interpolateVTK(orig_data, sub_grid):
    pfilter = vtkProbeFilter()
    pfilter.SetSourceData(orig_data)
    pfilter.SetInputData(sub_grid)
    pfilter.Update()
    return pfilter.GetStructuredGridOutput()


@click.command()
@click.option('--savefigure/--no-savefigure', default=True, type=bool, help='Export results to png file')
@click.option('--showplots/--no-showplots', default=False, type=bool, help='Show plots from matplotlib')
@click.option('--i_date', default=0, type=int, help='Index of date (0-29) in November')
def main(savefigure,showplots,i_date):
    # Collect metainformation (path,name,date) from SIMRA files
    home = expanduser("~")
    simraResultsFolder = home+'/results/simra/Sula/'
    openFoamResultsFolder = home+'/results/openfoam/Sula/'
    fileNamesOrg = [
        openFoamResultsFolder+'2020111906/steady_3.vtu',
        openFoamResultsFolder+'2020111906/steadyIsothermal_3.vtu',
        openFoamResultsFolder+'2020111906/unSteady_3.vtu',
        simraResultsFolder+'met/2020111906/metmesh_3.vtu',
    ]
    metFiles = glob(simraResultsFolder+"met_new/20*/*.vts")

    fileNamesOrg.extend(metFiles)
    df = pd.DataFrame({'path':fileNamesOrg})
    df['name'] = [path.split('/')[-1].split('.')[0].split('_')[0] for path in df.path]
    df['addtime'] = [path.split('/')[-1].split('.')[0].split('_')[-1] for path in df.path]
    baseDate = [datetime.strptime(path.split('/')[-2],'%Y%m%d%H') for path in df.path]
    df['baseDate'] = baseDate 
    df['date'] = [date+timedelta(hours=float(addtime)) for date,addtime in zip(baseDate,df.addtime)]
    df = df.sort_values(by='date')
    df = df.reset_index()
    uniqueDates = df['date'].map(pd.Timestamp).unique()

    # Get sensor locations
    sensorLoc,CoordUTM32,mastb,masth,Sensorh = computeSensorLoc(originx=-200,originy=6899800)
    noMasts = len(sensorLoc)
    curves = [''] * noMasts 
    for i in range(noMasts):
        curves[i] = interpolatePoints(sensorLoc[i])

    # Initiate plot arrays
    dataTypes = ['raw','rawMid','rawNew','rawMidNew']
    colorsData = [[0,0,0],[0.4,0.4,0.4],[0.8,0.4,0.8],[0.4,0.8,0.8]] 
    mastNames = ['Kvitneset', 'Traelboneset','Langeneset','Kaarsteinen']
    layoutNames = ['VelocityProfiles', 'WindDirProfiles', 'alphaProfiles']
    xArrayNames = ['u_mag', 'meandir', 'alpha']
    bottomAxisRangeMinimums = [0.0, 0.0, -30.0]
    bottomAxisRangeMaximums = [30.0, 2*np.pi, 30.0]
    bottomAxisTitles = ['$u$ [m/s] (magnitude)', 'Wind dir', 'Angle of Attack $[^\circ]$']
    
    #fig2 = plt.figure(2)
    #ax3d = fig2.add_subplot(111, projection='3d')
    
    noPlots = len(layoutNames)
    noDataTypes = len(dataTypes)

    print('no dates = '+str(len(uniqueDates)))
    uniqueDates = [uniqueDates[i_date]]

    # Loop over all dates and generate plots
    for date in uniqueDates:
        # Initiate figures
        fig = [''] * noPlots
        ax = [''] * noPlots
        datestr = pd.to_datetime(date).strftime('%Y-%m-%d %H:%M')
        print('Running '+datestr)
        for i_l in range(noPlots):
            usePolar = layoutNames[i_l] == 'WindDirProfiles'
            if usePolar:
                fig[i_l], ax[i_l] = plt.subplots(1, noMasts, sharey=True,figsize=(40,10), subplot_kw=dict(projection="polar"))
            else:
                fig[i_l], ax[i_l] = plt.subplots(1, noMasts, sharey=True,figsize=(20,10))
            suptitle = layoutNames[i_l]+' '+datestr
            fig[i_l].suptitle(suptitle)

        # Loop over all masts and data types
        for i in range(noMasts):
            for j in range(noDataTypes):
                # Collect data from all sensors
                df_obs = pd.DataFrame()
                noSensors = len(Sensorh[i])
                for k in range(noSensors):
                    z = np.floor(Sensorh[i][k]).astype(int)
                    filename = simraResultsFolder+'measurements/'+dataTypes[j]+'/10hz_'+mastNames[i]+'_60mnturbulence_statistics_'+str(z)+'_202011.csv'
                    df_all = pd.read_csv(filename)
                    if j == 0:
                        df_obs = df_obs.append(df_all[df_all.date==0])
                    indices = np.array(df_all.date,dtype='datetime64[m]') == np.array(date,dtype='datetime64[m]')
                    if df_all[indices].empty:
                        df_obs = df_obs.append(pd.Series(dtype='object'), ignore_index=True)
                    else:
                        df_obs = df_obs.append(df_all[indices])
                if df_obs.empty:
                    continue

                # Plot experimental data
                points = sensorLoc[i] 
                for i_l in range(noPlots):
                    QoI = df_obs[xArrayNames[i_l]]
                    if layoutNames[i_l] == 'WindDirProfiles':
                        QoI = np.radians(QoI)
                    ax[i_l][i].plot(QoI,points[:,2],color = colorsData[j],marker='.',label = 'Exp. '+dataTypes[j])
                    ax[i_l][i].set_xlim([bottomAxisRangeMinimums[i_l],bottomAxisRangeMaximums[i_l]])
                    ax[i_l][i].set_ylim([0,120])
                    ax[i_l][i].set_title(mastNames[i])
                    if j == 0:
                        if i == 0 and not layoutNames[i_l] == 'WindDirProfiles':
                            ax[i_l][i].set(xlabel=bottomAxisTitles[i_l], ylabel='$z$ [m]')
                        else:
                            ax[i_l][i].set(xlabel=bottomAxisTitles[i_l])
                    
                        if layoutNames[i_l] == 'WindDirProfiles':
                            ax[i_l][i].set_theta_zero_location('N')
                            ax[i_l][i].set_theta_direction(-1)

        # Iterate over all SIMRA files corresponding to the given date
        df_sub = df[df.date == date]
        i_df = 0
        dfSensor = pd.DataFrame()
        colorsCases = cm.jet(range(256))[0::256//len(df_sub)]
        dfSensor = pd.DataFrame()
        for _, df_date in df_sub.iterrows():
            # Load vtk files
            vtkData = load_vtk(df_date.path)
            datestr = pd.to_datetime(df_date['baseDate']).strftime('%Y%m%d%H')
            
            # Plot SIMRA results
            for i in range(noMasts):
                resampled = interpolateVTK(vtkData, curves[i])
                uUTM = vtk_to_numpy(resampled.GetPointData().GetAbstractArray('u')).copy()
                u_mag = np.linalg.norm(uUTM,axis=1)
                nodes = sensorLoc[i]
                points = vtk_to_numpy(curves[i].GetPoints().GetData()).copy()
                points = points[u_mag > 0,:]
                uUTM = uUTM[u_mag > 0,:]
                u_mag = u_mag[u_mag > 0]

                # Sample also at sensor locations
                resampledSensor = interpolateVTK(vtkData, createVTKcurve(sensorLoc[i]))
                uUTMSensor = vtk_to_numpy(resampledSensor.GetPointData().GetAbstractArray('u')).copy()
                u_magSensor = np.linalg.norm(uUTMSensor,axis=1)
                for i_l in range(noPlots):
                    QoI = getQoI(layoutNames[i_l],points,uUTM,u_mag)
                    ax[i_l][i].plot(QoI,points[:,2],color = colorsCases[i_df],label = df_date['name']+' '+pd.to_datetime(df_date['baseDate']).strftime('%Y%m%d%H')+'+'+df_date['addtime'])
                    if layoutNames[i_l] == 'WindDirProfiles':
                        ax[i_l][i].legend(loc='lower right')
                    else:
                        ax[i_l][i].legend(loc='upper left')

                    QoI = getQoI(layoutNames[i_l],points,uUTM,u_mag,useDeg=True)
                    dfCurve = pd.DataFrame({'z': points[:,2], xArrayNames[i_l]: QoI})
                    resultFileName = simraResultsFolder+'profileResults/'+layoutNames[i_l]+'_'+mastNames[i]+'_'+df_date['name']+'_'+datestr+'+'+df_date['addtime']+'.csv'
                    dfCurve.to_csv(resultFileName,index=False)

                for i_s in range(uUTMSensor.shape[0]):
                    row = pd.DataFrame({'date': [pd.to_datetime(date).strftime('%Y-%m-%d %H:%M')]})
                    row['name'] = df_date['name']
                    row['location'] = mastNames[i]
                    row['baseDate'] = df_date['baseDate']
                    row['addtime'] = df_date['addtime']
                    row['z'] = sensorLoc[i][i_s,-1]
                    for i_l in range(noPlots):
                        row[xArrayNames[i_l]] = getQoI(layoutNames[i_l],sensorLoc[i][[i_s],:],uUTMSensor[[i_s],:],u_magSensor[[i_s]],useDeg=True)

                    dfSensor = dfSensor.append(row)

            i_df += 1
            del vtkData

        resultFileName = simraResultsFolder+'sampledResults.csv'
        if isfile(resultFileName):
            dfSensor.to_csv(resultFileName,mode='a',index=False,header=False)
        else:
            dfSensor.to_csv(resultFileName,mode='a',index=False)


        if showplots:
            plt.show()

        if savefigure:
            for i_l in range(noPlots):
                fig[i_l].savefig(simraResultsFolder+'profileResults/'+layoutNames[i_l]+'_'+pd.to_datetime(date).strftime('%Y%m%d%H')+'.pdf', bbox_inches='tight',pad_inches = 0)

        plt.close()

if __name__ == '__main__':
    main()
