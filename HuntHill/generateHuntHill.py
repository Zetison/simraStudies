from itertools import zip_longest, chain
import click
import matplotlib.pyplot as plt
import numpy as np
from splipy.io import G2
import splipy.surface_factory as sfac
from splipy import SplineObject,BSplineBasis,Curve

def HuntHill(x,y,x0,y0,height):
    r1 = 0.203/height
    a1 = 0.076/height
    r = np.sqrt((x-x0)**2 + (y-y0)**2) / height
    return height * (1.04 / (1 + r**4) - 0.083/(1 + (r-r1)**2 / a1**2) - 0.03)

@click.command()
@click.option('--height', default=0.229, help='Height of Hunthill.')
@click.option('--ny', default=401,   help='Number of grid points in y-dir')
@click.option('--filename', default='HuntHill.g2', help='Output file name')
@click.option('--plotdistr/--no-plotdistr',     default=False, help='Plot distribution of stretching')
@click.option('--stretch/--no-stretch',         default=True, help='Toggle horizontal stretching of nodes')
@click.option('--interpolate/--no-interpolate', default=True, help='Use interpolation of Hunthill or a least square')
@click.option('--dcpts', default='0,0.6,0.6,0.6', help='Control points of the spline stretching distribution')
#def main(height, ny, plotdistr, stretch, interpolate, dcpts):
def main(height: float, ny: int, filename: str, plotdistr: bool, stretch: bool, interpolate: bool, dcpts: str):
    # Parameters from Utnes2007mos
    ymin = -5*height
    ymax = 5*height
    xmin = -6*height
    xmax = 7*height
    x0 = 0
    y0 = 0
    
    Nx = round((ny-1)*(xmax-xmin)/(ymax-ymin))+1
    
    p = 3
    Ncx = Nx-1+p
    Ncy = ny-1+p
    
    # Number of points for least square/interpolate
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
    def distr(xi):
        basis = BSplineBasis(order=3)
        basis.insert_knot(0.8)
        controlpoints = np.zeros((4,2))
        controlpoints[:,0] = np.linspace(0,1,4)
        controlpoints[:,1] = np.array(list(map(float, dcpts.split(','))))
        
        C = Curve(basis,controlpoints=controlpoints)
        if plotdistr:
            xiplot = np.linspace(0,1,1000)
            P = C(xiplot)
            plt.plot(controlpoints[:,0],controlpoints[:,1],'ko-',markerfacecolor='red')
            plt.plot(P[:,0],P[:,1])
            plt.show()

        shape = xi.shape
        xi = np.reshape(xi,shape[0]*shape[1])
        xi[xi<0] = 0
        return np.reshape(C(xi)[:,1],shape)
    
    def f(x,y):
        r_se = [xmax,ymin]
        r_ne = [xmax,ymax]
        r_nw = [xmin,ymax]
        r_sw = [xmin,ymin]
        r_0 = [x0,y0]
        maxD = -1
        for r in [r_se,r_ne,r_nw,r_sw]:
            D = np.linalg.norm(np.array(r)-r_0)
            if maxD < D:
                maxD = D
    
        Rx = x0-x
        Ry = y0-y
        temp = (1-np.sqrt(Rx**2+Ry**2)/maxD)*(x-xmin)*(xmax-x)*(y-ymin)*(ymax-y)/np.max([xmax-x0,x0-xmin])**2/np.max([ymax-y0,y0-ymin])**2
        temp = distr(temp)
        return x + Rx*temp, y + Ry*temp
    
    if stretch:
        xx,yy = f(xx,yy)
    
    zz = HuntHill(xx,yy,x0,y0,height)
    
    # Generate Greville points for least square
    U = np.pad(np.linspace(0,Ncx-p,Ncx-(p-1)),(p,p),'constant',constant_values=(0,Ncx-p))
    V = np.pad(np.linspace(0,Ncy-p,Ncy-(p-1)),(p,p),'constant',constant_values=(0,Ncy-p))
    basis_u = BSplineBasis(knots = U, order=p+1)
    basis_v = BSplineBasis(knots = V, order=p+1)
    if interpolate:
        patch = sfac.interpolate(np.swapaxes(np.array([xx,yy,zz]),0,2), [basis_u, basis_v])
    else:
        patch = sfac.least_square_fit(np.swapaxes(np.array([xx,yy,zz]),0,2), [basis_u, basis_v], [u,v])
    
    with G2(filename) as my_file:
        my_file.write(patch)
    
if __name__ == '__main__':
    main()
