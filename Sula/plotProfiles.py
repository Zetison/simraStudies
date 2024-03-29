import numpy as np
import matplotlib.pylab as plt
import netCDF4
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
from siso.coord import find_system, conversion_path
from siso.filter.coordtransform import CoordTransformCache
from siso.util import FieldData

def getAromeData(filename,Sensorh,mastb,date=False):
    # Get Arome data
    try:
        nc_Arome = netCDF4.Dataset(filename)
    except Exception:
        print("Data is lacking for "+mastNames[i])
        return pd.DataFrame()

    times = nc_Arome.variables['time']
    tAll = netCDF4.num2date(times[:], times.units,only_use_cftime_datetimes=False).data
    if date:
        indices = np.array(tAll,dtype='datetime64[m]') == np.array(date,dtype='datetime64[m]')
    else:
        indices = slice(len(times))

    w = nc_Arome.variables['upward_air_velocity'][0].data[indices]

    df0 = pd.DataFrame()
    df0['date'] = tAll[indices]
    df0['meanU'] = nc_Arome.variables['windspeed'][0].data[indices]
    meanU = df0['meanU'] 
    df0['u_mag'] = np.sqrt(meanU**2 + w**2) 
    df0['alpha'] = np.degrees(np.arctan2(w,meanU))
    df0['meandir'] = nc_Arome.variables['winddirection'][0].data[indices]

    return df0

def getQoI(name,utmPoints,uUTM,u_mag,useDeg=False):
    src = find_system('utm:33n')
    tgt = find_system('geodetic')
    path = conversion_path(src, tgt)
    my_points = FieldData(np.array(utmPoints))
    my_vectors = FieldData(np.array(uUTM))
    cache = CoordTransformCache(path)
    out_points = cache.convert_geometry(my_points, 0).numpy()
    u = cache.convert_vectors(my_vectors, 0).numpy()
    
    if name == 'VelocityProfiles':
        QoI = u_mag
    elif name == 'WindDirProfiles':
        if useDeg:
            QoI = (180+90-np.degrees(np.arctan2(u[:,1],u[:,0]))) % 360.
        else:
            QoI = (np.radians(180+90)-np.arctan2(u[:,1],u[:,0])) % (2*np.pi)
    elif name == 'alphaProfiles':
        QoI = np.degrees(np.arctan2(u[:,2],np.linalg.norm(u[:,:2],axis=1)))
    elif name == 'meanU':
        QoI = np.linalg.norm(u[:,:2],axis=1)
    elif name == 'meanu':
        QoI = u[:,0]
    elif name == 'meanv':
        QoI = u[:,1]
    elif name == 'meanw':
        QoI = u[:,2]
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
@click.option('--savefigure/--no-savefigure', default=False, type=bool, help='Export results to png file')
@click.option('--showplots/--no-showplots', default=True, type=bool, help='Show plots from matplotlib')
@click.option('--loadvtk/--no-loadvtk', default=False, type=bool, help='Reload vtk files or the existing generated .csv files')
@click.option('--i_date', default=0, type=int, help='Index of date (0-29) in November')
@click.option('--plot_simraresults/--no-plot_simraresults', default=False, type=bool, help='Plot SIMRA results')
@click.option('--plot_wrfresults/--no-plot_wrfresults', default=True, type=bool, help='Plot WRF results')
def main(savefigure,showplots,loadvtk,i_date,plot_simraresults,plot_wrfresults):
    # Collect metainformation (path,name,date) from SIMRA files
    home = expanduser("~")
    simraResultsFolder = home+'/results/simra/Sula/'
    wrfResultsFolder = home+'/results/WRF/Sula/'
    openFoamResultsFolder = home+'/results/openfoam/Sula/'
    fileNamesOrg = [
        openFoamResultsFolder+'2020111906/steady_3.vtu',
        openFoamResultsFolder+'2020111906/steadyIsothermal_3.vtu',
        openFoamResultsFolder+'2020111906/unSteady_3.vtu',
        simraResultsFolder+'met/2020111906/metmesh_3.vtu',
    ]
    fileNamesOrg = []
    if plot_simraresults:
        metFiles = glob(simraResultsFolder+"met_new/20*/*.vts")
    if plot_wrfresults:
        wrfFiles = glob(wrfResultsFolder+"20*/*.vts")

    if plot_simraresults:
        fileNamesOrg.extend(metFiles)
    
    if plot_wrfresults:
        fileNamesOrg.extend(wrfFiles)

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
    if plot_wrfresults:
        originx = 0
        originy = 0
    else:
        originx = -200
        originy = 6899800

    sensorLoc,CoordUTM32,mastb,masth,Sensorh = computeSensorLoc(originx,originy)
    noLocs = len(sensorLoc)
    curves = [''] * noLocs
    for i in range(noLocs):
        curves[i] = interpolatePoints(sensorLoc[i])

    # Initiate plot arrays
    #dataTypes = ['raw','rawMid','rawNew','rawMidNew']
    dataTypes = ['rawMidNew']
    colorsData = [[0,0,0],[0.4,0.4,0.4],[0.8,0.4,0.8],[0.4,0.8,0.8]] 
    mastNames = ['Kvitneset', 'Traelboneset','Langeneset','Kaarsteinen', 'Bridgecenter']
    layoutNames = ['VelocityProfiles', 'WindDirProfiles', 'alphaProfiles']
    xArrayNames = ['u_mag', 'meandir', 'alpha', 'meanU']
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
                fig[i_l], ax[i_l] = plt.subplots(1, noLocs, sharey=True,figsize=(40,10), subplot_kw=dict(projection="polar"))
            else:
                fig[i_l], ax[i_l] = plt.subplots(1, noLocs, sharey=True,figsize=(20,10))
            suptitle = layoutNames[i_l]+' '+datestr
            fig[i_l].suptitle(suptitle)

        # Loop over all masts and data types
        for i in range(noLocs):
            for j in range(noDataTypes):
                # Collect data from all sensors
                df_obs = pd.DataFrame()
                df_arome = pd.DataFrame()
                noSensors = len(Sensorh[i])
                for k in range(noSensors):
                    z = np.floor(Sensorh[i][k]).astype(int)
                    filename = simraResultsFolder+'measurements/'+dataTypes[j]+'/10hz_'+mastNames[i]+'_60mnturbulence_statistics_'+str(z)+'_202011.csv'
                    df_all = pd.read_csv(filename)
                    #if j == 0:
                    #    df_obs = df_obs.append(df_all[df_all.date==0])
                    indices = np.array(df_all.date,dtype='datetime64[m]') == np.array(date,dtype='datetime64[m]')
                    if df_all[indices].empty:
                        df_obs = pd.concat([df_obs,pd.Series(dtype='object')], ignore_index=True)
                    else:
                        df_obs = pd.concat([df_obs,df_all[indices]])

                    absHeight = str(Sensorh[i][k]+mastb[i]).rstrip('0').rstrip('.')
                    filename = home+'/results/simra/Sula/measurements/202011_%s_%sm_mepsdetml.nc' % (mastNames[i],absHeight)
                    df_arome = pd.concat([df_arome,getAromeData(filename,Sensorh,mastb,date)])

                if df_obs.empty:
                    continue

                # Plot experimental data
                points = sensorLoc[i] 
                for i_l in range(noPlots):
                    QoI = df_obs[xArrayNames[i_l]]
                    if layoutNames[i_l] == 'WindDirProfiles':
                        QoI = np.radians(QoI)
                    # ax[i_l][i].plot(QoI,points[:,2],color = colorsData[j],marker='.',label = 'Exp. '+dataTypes[j])
                    if mastNames[i] == 'Bridgecenter':
                        ax[i_l][i].plot(QoI,points[:,2],color = colorsData[j],marker='.',label = 'Lidar')
                    else:
                        ax[i_l][i].plot(QoI,points[:,2],color = colorsData[j],marker='.',label = 'Anemometers')

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

                    if layoutNames[i_l] == 'alphaProfiles' or mastNames[i] == 'Bridgecenter':
                        continue
                    else:
                        QoI = df_arome[xArrayNames[i_l]]
                        if layoutNames[i_l] == 'WindDirProfiles':
                            QoI = np.radians(QoI)

                        ax[i_l][i].plot(QoI,points[:,2],color = [1.0,0.0,0.0],marker='x',label = 'Arome')

        # Iterate over all SIMRA files corresponding to the given date
        df_sub = df[df.date == date]
        i_df = 0
        dfSensor = pd.DataFrame()
        colorsCases = cm.jet(range(256))[0::256//len(df_sub)]
        dfSensor = pd.DataFrame()
        for _, df_date in df_sub.iterrows():
            if loadvtk:
                # Load vtk files
                vtkData = load_vtk(df_date.path)

            datestr = pd.to_datetime(df_date['baseDate']).strftime('%Y%m%d%H')
            
            # Plot SIMRA results
            for i in range(noLocs):
                if loadvtk:
                    resampled = interpolateVTK(vtkData, curves[i])
                    try:
                        uUTM = vtk_to_numpy(resampled.GetPointData().GetAbstractArray('u')).copy()
                    except:
                        uUTM = vtk_to_numpy(resampled.GetPointData().GetAbstractArray('WIND')).copy()
                    u_mag = np.linalg.norm(uUTM,axis=1)
                    nodes = sensorLoc[i]
                    points = vtk_to_numpy(curves[i].GetPoints().GetData()).copy()
                    points = points[u_mag > 0,:]
                    uUTM = uUTM[u_mag > 0,:]
                    u_mag = u_mag[u_mag > 0]

                    # Sample also at sensor locations
                    resampledSensor = interpolateVTK(vtkData, createVTKcurve(sensorLoc[i]))
                    try:
                        uUTMSensor = vtk_to_numpy(resampledSensor.GetPointData().GetAbstractArray('u')).copy()
                    except:
                        uUTMSensor = vtk_to_numpy(resampledSensor.GetPointData().GetAbstractArray('WIND')).copy()
                    u_magSensor = np.linalg.norm(uUTMSensor,axis=1)

                for i_l in range(noPlots):
                    if loadvtk:
                        QoI = getQoI(layoutNames[i_l],points,uUTM,u_mag)
                        z = points[:,2]
                    else:
                        resultFileName = simraResultsFolder+'profileResults/'+layoutNames[i_l]+'_'+mastNames[i]+'_'+df_date['name']+'_'+datestr+'+'+df_date['addtime']+'.csv'
                        dfCurve = pd.read_csv(resultFileName)
                        z = dfCurve['z']
                        QoI = dfCurve[xArrayNames[i_l]]

                    ax[i_l][i].plot(QoI,z,color = colorsCases[i_df],label = df_date['name']+' '+pd.to_datetime(df_date['baseDate']).strftime('%Y%m%d%H')+'+'+df_date['addtime'])

                    if layoutNames[i_l] == 'WindDirProfiles':
                        ax[i_l][i].legend(loc='lower right')
                    else:
                        ax[i_l][i].legend(loc='upper left')

                    if loadvtk:
                        resultFileName = simraResultsFolder+'profileResults/'+layoutNames[i_l]+'_'+mastNames[i]+'_'+df_date['name']+'_'+datestr+'+'+df_date['addtime']+'.csv'
                        QoI = getQoI(layoutNames[i_l],points,uUTM,u_mag,useDeg=True)
                        dfCurve = pd.DataFrame({'z': points[:,2], xArrayNames[i_l]: QoI})
                        dfCurve.to_csv(resultFileName,index=False)

                if loadvtk:
                    for i_s in range(uUTMSensor.shape[0]):
                        row = pd.DataFrame({'date': [pd.to_datetime(date).strftime('%Y-%m-%d %H:%M')]})
                        row['name'] = df_date['name']
                        row['location'] = mastNames[i]
                        row['baseDate'] = df_date['baseDate']
                        row['addtime'] = df_date['addtime']
                        row['z'] = sensorLoc[i][i_s,-1]
                        layoutNamesOut = layoutNames + ['meanU']
                        for i_l in range(len(layoutNamesOut)):
                            row[xArrayNames[i_l]] = getQoI(layoutNamesOut[i_l],sensorLoc[i][[i_s],:],uUTMSensor[[i_s],:],u_magSensor[[i_s]],useDeg=True)

                        dfSensor = pd.concat([dfSensor,row])

            i_df += 1
            if loadvtk:
                del vtkData

        resultFileName = simraResultsFolder+'sampledResults.csv'
        if loadvtk:
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
