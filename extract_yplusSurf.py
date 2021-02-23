from os.path import expanduser,exists
import click
import sys
import numpy as np
import vtk
from splipy.io import G2
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy
import vtk.util.numpy_support as vtknp
import splipy.surface_factory as sfac
from splipy import SplineObject,BSplineBasis

@click.command()
@click.option('--resultsfolder', type=str, help='Path of the folder containing the results')
@click.option('--nelu', type=int, help='Number of elements in u-dir')
@click.option('--nelv', type=int, help='Number of elements in v-dir')
@click.option('--nelw', type=int, help='Number of elements in w-dir')
@click.option('--nstep', type=int, help='Number of elements in w-dir')
@click.option('--savehist', type=int, help='Number of elements in w-dir')
@click.option('--pvdname', default='HuntHill', type=str, help='Name of pvd file')

def main(resultsfolder: str, nelu: int, nelv: int, nelw: int, nstep: int, savehist: int, pvdname: str):
    n_u = nelu+1
    n_v = nelv+1
    n_w = nelw+1
    
    reader = vtk.vtkXMLUnstructuredGridReader()  # or vtkXMLStructuredGridReader, or vtkUnstructuredReader, or vtkStructuredReader
    reader.SetFileName(resultsfolder+'/'+pvdname+'.pvd-data/data-'+str(nstep//savehist)+'.vtu')
    reader.Update()
    grid = reader.GetOutput()
    yplus = vtk_to_numpy(grid.GetPointData().GetArray('y+'))
    tk = vtk_to_numpy(grid.GetPointData().GetArray('tk'))
    X = vtk_to_numpy(grid.GetPoints().GetData())
    ptsarray = X[n_u*n_v:2*n_u*n_v]
    
    points = vtk.vtkPoints()
    points.SetData(vtknp.numpy_to_vtk(ptsarray))
    
    grid = vtk.vtkStructuredGrid()
    grid.SetDimensions(n_u,n_v,1)
    grid.SetPoints(points)
    
    data = vtknp.numpy_to_vtk(yplus[n_u*n_v:2*n_u*n_v])
    data.SetName('y+')
    grid.GetPointData().AddArray(data)
    data_tk = vtknp.numpy_to_vtk(tk[n_u*n_v:2*n_u*n_v])
    data_tk.SetName('tk')
    grid.GetPointData().AddArray(data_tk)
    
    X = np.reshape(X,(n_w,n_v,n_u,3))
    yplus = np.reshape(yplus,(n_w,n_v,n_u))
    controlpoints = np.copy(X[1,:,:,:])
    controlpoints[:,:,2] = yplus[1,:,:]
    controlpoints = np.transpose(controlpoints, (1,0,2))
    
    p = 1
    U = np.pad(np.linspace(0,n_u-p,n_u-(p-1)),(p,p),'constant',constant_values=(0,n_u-p))
    V = np.pad(np.linspace(0,n_v-p,n_v-(p-1)),(p,p),'constant',constant_values=(0,n_v-p))
    basis_u = BSplineBasis(knots = U, order=p+1)
    basis_v = BSplineBasis(knots = V, order=p+1)
    patch = sfac.interpolate(controlpoints, [basis_u, basis_v])
    i = 0
    while True:
        yplusfilename = 'yplus'+str(i)+'.g2'
        if exists(yplusfilename):
            i += 1
        else:
            print('Writing '+yplusfilename)
            with G2(yplusfilename) as my_file:
                my_file.write(patch)
    
            break
    
    writer = vtk.vtkStructuredGridWriter()
    writer.SetFileName('yplus'+str(i)+'.vtk')
    writer.SetInputData(grid)
    writer.Write()

if __name__ == '__main__':
    main()
