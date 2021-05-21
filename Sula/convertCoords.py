from pyproj import Transformer
import pandas as pd
import vtk
from pyevtk.hl import pointsToVTK
from vtk.util.numpy_support import vtk_to_numpy
import vtk.util.numpy_support as vtknp
import click
import numpy as np
@click.command()
@click.option('--originx', default=0.0, type=float, help='x-coordinate of origin')
@click.option('--originy', default=0.0, type=float, help='y-coordinate of origin')
@click.option('--date', default='2020-11-19 06:00', type=str, help='Date for data extraction')
@click.option('--measurementfolder', default='measurements', type=str, help='Folder for the experimental data')
@click.option('--resultsfolder', default='measurements', type=str, help='Folder for the results')
def main(originx,originy,date,measurementfolder,resultsfolder):
    proj32 = "+proj=utm +zone=32K, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    proj33 = "+proj=utm +zone=33K, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    transformer = Transformer.from_crs(proj32,proj33)

    height = 70 # Height of bridge at middle point
    # Data for 'SulaNW (Kvitneset)', 'SulaNE (Trælboneset)', 'SulaSW (Langeneset)', 'SulaSE (Kårsteinen)'
    Boom1 = np.array([6.1,6.1,4.4,3.6])
    Sensorh = np.array([[92.5,71.5,44.5],[76.8, 48.3, 27.3],[94.8, 75.0, 50.0, 27.0],[62.8, 40.0, 13.4]],dtype=object)
    BoomOrient = np.array([[72,74,74],[289, 290, 290],[81, 81, 81, 81],[223, 223, 223]],dtype=object)
    masth = np.array([100.5, 78.0,97.0,63.0])
    if True:
        # Based on (lat,lon) coords from .nc file: (62.421595,6.00112),(62.42763,6.062626),(62.386301,6.031318), (62.400134,6.119395)
        CoordUTM32 = np.array([[6924741.06, 345141.99],[6925267.00,348347.02],[6920740.04,346519.99],[6922073.99,351139.99]])
        mastb = np.array([6.0,14.0,6.0,12.0])
    else:
        CoordUTM32 = np.array([[6924741, 345142],[6925267, 348347],[6920740, 346520],[6922074, 351140]])
        CoordUTM32 = np.array([[6924741, 345142],[6925263.08,348346.72],[6920736.2,346525.95],[6922076.42,351140.03]])
        mastb = np.array([8.0,12.2,1.8,13.9])

    noMasts = len(masth)
    sensorLoc = [''] * noMasts
    f = open("POINT_FILE", "w")
    noPoints = sum([len(listElem) for listElem in Sensorh])+len(Sensorh)//2
    f.write("%d\n" % noPoints)
    for i in range(noMasts):
        noSensors = len(Sensorh[i])
        sensorLoc[i] = np.zeros((noSensors,3))
        for j in range(noSensors):
            alpha = np.radians(BoomOrient[i][j])
            x, y = transformer.transform(CoordUTM32[i,1]+Boom1[i]*np.sin(alpha),CoordUTM32[i,0]+Boom1[i]*np.cos(alpha)) # This is probably not the exact formula due to the UTM transformation?
            sensorLoc[i][j,:] = [x-originx,y-originy,mastb[i]+Sensorh[i][j]]
            f.write("%f %f %f\n" % tuple(sensorLoc[i][j,:]))

    for i in range(0,CoordUTM32.shape[0]//2):
        x1, y1 = transformer.transform(CoordUTM32[2*i,1],CoordUTM32[2*i,0])
        x2, y2 = transformer.transform(CoordUTM32[2*i+1,1],CoordUTM32[2*i+1,0])
        f.write("%f %f %f\n" % ((x1+x2)/2-originx,(y1+y2)/2-originy,height))

    f.close()

    f = open("mastLoc.csv", "w")
    for i in range(len(masth)):
        x, y = transformer.transform(CoordUTM32[i,1],CoordUTM32[i,0])
        f.write("%f,%f,%f,%f\n" % (x-originx,y-originy,mastb[i],masth[i]))

    f.close()
    
    mastNames = ['Kvitneset', 'Traelboneset','Langeneset','Kaarsteinen']
    for i in range(noMasts):
        df = pd.DataFrame()
        noSensors = len(Sensorh[i])
        for j in range(noSensors):
            z = np.floor(Sensorh[i][j]).astype(int)
            filename = measurementfolder+'/'+'10hz_'+mastNames[i]+'_60mnturbulence_statistics_'+str(z)+'_202011.csv'
            df_all = pd.read_csv(filename)
            if j == 0:
                df = df.append(df_all[df_all.date==0])

            if df_all[df_all.date==date].empty:
                df = df.append(pd.Series(dtype='object'), ignore_index=True)
            else:
                df = df.append(df_all[df_all.date==date])
        df['coordsZ'] = mastb[i]+Sensorh[i]
        print(date)
        #grid = vtk.vtkStructuredGrid()
        #grid.SetDimensions(noSensors,1,1)
        #points = vtk.vtkPoints()
        #points.SetData(vtknp.numpy_to_vtk(sensorLoc[i]))
        #grid.SetPoints(points)
        #for column in df:
        #    if not column == 'index' and not column == 'date':
        #        data = vtknp.numpy_to_vtk(df[column])
        #        data.SetName(column)
        #        grid.GetPointData().AddArray(data)

        #writer = vtk.vtkStructuredGridWriter()# or vtkXMLStructuredGridReader, or vtkUnstructuredReader, or vtkStructuredReader
        #writer.SetFileName('measurements/VelocityProfile_'+mastNames[i]+'_'+date[11:13]+'.vtk')
        #writer.SetInputData(grid)
        #writer.Write()
        X = np.array(sensorLoc[i][:,0])
        Y = np.array(sensorLoc[i][:,1])
        Z = np.array(sensorLoc[i][:,2])
        data = df.drop(['date'], axis=1).to_dict(orient='list')
        for column in data:
            data[column] = np.array(data[column])

        pointsToVTK(resultsfolder+'/VelocityProfile_'+mastNames[i]+'_'+str(int(date[11:13])), X, Y, Z, data = data)


if __name__ == '__main__':
    main()
