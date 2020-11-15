from itertools import zip_longest, chain
import numpy as np
from splipy.io import G2
import splipy.surface_factory as sfac
from splipy import SplineObject
from splipy import BSplineBasis

# Parameters from Utnes2007mos
H = 0.229
r1 = 0.203 / H
a1 = 0.076 / H

xmin = -6 * H
xmax = 7 * H
ymin = -5 * H
ymax = 5 * H

# Number of points in each direction
Ny = 89
#Ny = 11
Nx = round(Ny*(xmax-xmin)/(ymax-ymin))
#Nx = 89

p = 3
Ncx = Nx-1+p
Ncy = Ny-1+p

# Number of points for least square/interpolate
interpolate = True
if interpolate:
	NcX = Ncx
	NcY = Ncy
else:
	NcX = 10*Ncx
	NcY = 10*Ncy
U = np.pad(np.linspace(0,Ncx-p,NcX-(p-1)),(p,p),'constant',constant_values=(0,Ncx-p))
V = np.pad(np.linspace(0,Ncy-p,NcY-(p-1)),(p,p),'constant',constant_values=(0,Ncy-p))
basis_U = BSplineBasis(knots = U, order=p+1)
basis_V = BSplineBasis(knots = V, order=p+1)
u = basis_U.greville()
v = basis_V.greville()

def u2x(u,xmin,xmax):
	return xmin + (xmax-xmin)*(np.array(u)-u[0])/(u[-1]-u[0])

x = u2x(u,xmin,xmax)
y = u2x(v,ymin,ymax)
xx, yy = np.meshgrid(x, y)

def HuntHill(x,y):
	r = np.sqrt(x**2 + y**2) / H
	return H * (1.04 / (1 + r**4) - 0.083/(1 + (r-r1)**2 / a1**2) - 0.03)

#zz = HuntHill(xx,yy)
zz = HuntHill(xx,yy)-HuntHill(xmin,ymin)
#zz = np.zeros(xx.shape)

# Generate Greville points for least square

# Create patch
U = np.pad(np.linspace(0,Ncx-p,Ncx-(p-1)),(p,p),'constant',constant_values=(0,Ncx-p))
V = np.pad(np.linspace(0,Ncy-p,Ncy-(p-1)),(p,p),'constant',constant_values=(0,Ncy-p))
basis_u = BSplineBasis(knots = U, order=p+1)
basis_v = BSplineBasis(knots = V, order=p+1)
if interpolate:
	patch = sfac.interpolate(np.swapaxes(np.array([xx,yy,zz]),0,2), [basis_u, basis_v])
else:
	patch = sfac.least_square_fit(np.swapaxes(np.array([xx,yy,zz]),0,2), [basis_u, basis_v], [u,v])

with G2('HuntHill.g2') as my_file:
	my_file.write(patch)
