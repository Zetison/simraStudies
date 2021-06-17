from scipy.io import FortranFile
import numpy as np

m = FortranFile('mesh_new.dat-wSula', 'r', header_dtype='>u4')
grid_ints = m.read_ints(dtype='>u4')
grid_coords = m.read_reals(dtype='>f4')
grid_elems = m.read_ints(dtype='>u4')
om = FortranFile('mesh_converted.dat', 'w')
print(grid_ints)
print(len(grid_coords))
print(len(grid_elems))
om.write_record(np.array(grid_ints,dtype=np.int32))
om.write_record(np.array(grid_coords,dtype=np.float32))
om.write_record(np.array(grid_elems,dtype=np.int32))

f = FortranFile('contwSula.res_00+3', 'r', header_dtype='>u4')
of = FortranFile('cont_met.res', 'w')
while True:
    try:
        data = f.read_reals('>f4')
        of.write_record(np.array(data, dtype=np.float32))
    except FortranFormattingError:
        break
