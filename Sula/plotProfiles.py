import numpy as np
import matplotlib.pylab as plt
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import sys
import click
from vtk import vtkXMLStructuredGridReader, vtkPoints, vtkStructuredGrid, vtkProbeFilter
from vtk import vtkXMLStructuredGridWriter
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy
import utm
from glob import glob
from os.path import expanduser
from convertCoords import computeSensorLoc 

def load_vts(filename):
    reader = vtkXMLStructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()


def interpolatePoints(nodes,noPts=1000,z_min=0.0,z_max=120.0):
    nodes = np.vstack([np.append(nodes[0,:-1],z_max), nodes,np.append(nodes[-1,:-1],z_min)])
    num_true_pts = nodes.shape[0]

    tck, u = interpolate.splprep(np.transpose(nodes),k=num_true_pts-1)
    x_knots, y_knots, z_knots = interpolate.splev(tck[0], tck)
    u_fine = np.linspace(0,1,noPts)
    x_fine, y_fine, z_fine = interpolate.splev(u_fine, tck)

    points = np.transpose([x_fine, y_fine, z_fine])

    vpoints = vtkPoints()
    vpoints.SetData(numpy_to_vtk(points, deep=True))

    vgrid = vtkStructuredGrid()
    vgrid.SetDimensions(noPts,1,1)
    vgrid.SetPoints(vpoints)

    return vgrid

def interpolateVTK(orig_data, sub_grid):
    pfilter = vtkProbeFilter()
    pfilter.SetSourceData(orig_data)
    pfilter.SetInputData(sub_grid)
    pfilter.Update()
    return pfilter.GetStructuredGridOutput()


@click.command()
#@click.option('--originx', default=0.0, type=float, help='x-coordinate of origin')
#@click.option('--originy', default=0.0, type=float, help='y-coordinate of origin')
#@click.option('--date', default='2020-11-19 06:00', type=str, help='Date for data extraction')
#@click.option('--measurementfolder', default='measurements', type=str, help='Folder for the experimental data')
#@click.option('--resultsfolder', default='measurements', type=str, help='Folder for the results')
def main():
    home = expanduser("~")
    simraResultsFolder = home+'/results/simra/Sula/'
    openFoamResultsFolder = home+'/results/openfoam/Sula/'
    fileNamesOrg = [
        openFoamResultsFolder+'2020111906_OF_unStdy.pvd',
        openFoamResultsFolder+'2020111906_OF_stdy.pvd',
        openFoamResultsFolder+'2020111906_OF_stdyIso.pvd',
        simraResultsFolder+'met/cont_met.pvd',
    ]
    metFiles = glob(simraResultsFolder+"met_new/2*.pvd")

    fileNamesOrg.extend(metFiles)

    caseNamesOrg = [
    '2020111906_OF_unStdy',
    '2020111906_OF_stdy',
    '2020111906_OF_stdyIso',
    '2020111906_met',
    ]
    caseNamesOrg.extend([f[-17:-4] for f in metFiles])

    #indices = [2,5,6,7,8]
    indices = range(9) 
    #indices = [0,3]
    indices = [0,1,6]
    #indices = range(len(fileNamesOrg))
    fileNames = [fileNamesOrg[i] for i in indices]
    caseNames = [caseNamesOrg[i] for i in indices]
    
    # Get sensor locations
    sensorLoc,CoordUTM32,mastb,masth = computeSensorLoc(originx=-200,originy=6899800)
    noMasts = len(sensorLoc)
    curves = [''] * noMasts 
    for i in range(noMasts):
        curves[i] = interpolatePoints(sensorLoc[i])

    # Load vts files
    #orig_data = load_vts(home+'/results/simra/Sula/met_new/2020111906/M0_0.vts')
    orig_data = load_vts(home+'/kode/simraStudies/Sula/small.vts')
    resampled = [''] * noMasts
    for i in range(noMasts):
        resampled[i] = interpolateVTK(orig_data, curves[i])
    
    # Load experimental data

    vpcsv[i][j] = PVDReader(registrationName=mastNames[i], FileName=fileName+mastNames[i]+'.pvd')

    # Plot data
    #dataTypes = ['filtered','raw','rawMid']
    dataTypes = ['raw','rawMid']
    colorsData = [[0,0,0],[0.4,0.4,0.4],[0.8,0.8,0.8]] 
    mastNames = ['Kvitneset', 'Traelboneset','Langeneset','Kaarsteinen']
    layoutNames = ['VelocityProfiles', 'WindDirProfiles', 'alphaProfiles']
    layoutNames = ['VelocityProfiles']
    xArrayNames = ['u_mag', 'meandir', 'alpha']
    bottomAxisRangeMinimums = [0.0, 0.0, -10.0]
    bottomAxisRangeMaximums = [30.0, 360.0, 10.0]
    bottomAxisTitles = ['$u$ [m/s] (magnitude)', 'Wind dir', 'Angle of Attack $[^\circ]$']
    #fig2 = plt.figure(2)
    #ax3d = fig2.add_subplot(111, projection='3d')

    
    noPlots = len(layoutNames)
    fig = [''] * noPlots
    ax = [''] * noPlots
    for i_l in range(noPlots):
        fig[i_l], ax[i_l] = plt.subplots(1, noMasts, sharey=True)
        fig[i_l].suptitle(layoutNames[i_l]) 

    for i in range(noMasts):
        u = vtk_to_numpy(resampled[i].GetPointData().GetAbstractArray('u')).copy()
        nodes = sensorLoc[i]
        points = vtk_to_numpy(curves[i].GetPoints().GetData()).copy()
        for i_l in range(noPlots):
            ax[i_l][i].plot(np.linalg.norm(u,axis=1),points[:,2],color = [0,1,0],label = 'SIMRA')
            ax[i_l][i].set_xlim([bottomAxisRangeMinimums[i_l],bottomAxisRangeMaximums[i_l]])
            ax[i_l][i].set_ylim([0,120])
            ax[i_l][i].set_title(mastNames[i])
            if i == 0:
                ax[i_l][i].set(xlabel=bottomAxisTitles[i_l], ylabel='$z$ [m]')
            else:
                ax[i_l][i].set(xlabel=bottomAxisTitles[i_l])
            ax[i_l][i].legend()
            #ax3d.plot(nodes[:,0],nodes[:,1],nodes[:,2], 'ro')
            #ax3d.plot(points[:,0], points[:,1], points[:,2], 'r')

    
    plt.show()

if __name__ == '__main__':
    main()
