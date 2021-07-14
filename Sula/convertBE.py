from scipy.io import FortranFile
import numpy as np
import click

@click.command()
@click.option('--mesh_input', default='mesh_new.dat-wSula', type=str, help='Input simra mesh.dat file in big endian format to be converted')
@click.option('--cont_input', default='contwSula.res_00+3', type=str, help='Input simra cont.res file in big endian format to be converted')
@click.option('--mesh_output', default='mesh.dat', type=str, help='Output simra mesh.dat file in little endian format to be converted')
@click.option('--cont_output', default='cont.res', type=str, help='Output simra cont.res file in little endian format to be converted')
def main(mesh_input,cont_input,mesh_output,cont_output):
    m = FortranFile(mesh_input, 'r', header_dtype='>u4')
    grid_ints = m.read_ints(dtype='>u4')
    grid_coords = m.read_reals(dtype='>f4')
    grid_elems = m.read_ints(dtype='>u4')
    om = FortranFile(mesh_output, 'w')
    print(grid_ints)
    print(len(grid_coords))
    print(len(grid_elems))
    om.write_record(np.array(grid_ints,dtype=np.int32))
    om.write_record(np.array(grid_coords,dtype=np.float32))
    om.write_record(np.array(grid_elems,dtype=np.int32))

    f = FortranFile(cont_input, 'r', header_dtype='>u4')
    of = FortranFile(cont_output, 'w')
    while True:
        try:
            data = f.read_reals('>f4')
            of.write_record(np.array(data, dtype=np.float32))
        except Exception:
            break

if __name__ == '__main__':
    main()
