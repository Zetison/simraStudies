import click
from typing import IO
import vtk

@click.command()
@click.argument('infile', type=click.File('r'))

def main(infile: IO):
    reader = vtk.vtkUnstructuredGridReader()  # or vtkXMLStructuredGridReader, or vtkUnstructuredReader, or vtkStructuredReader
    reader.SetFileName(infile.name)
    reader.Update()
    grid = reader.GetOutput()
    for oldName, newName in zip(['k','U'],['tk','u']): 
        data_tk = grid.GetPointData().GetArray(oldName)
        data_tk.SetName(newName)
    
    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetFileName(infile.name[:-4]+'.vtk')
    writer.SetInputData(grid)
    writer.Write()

if __name__ == '__main__':
    main()
