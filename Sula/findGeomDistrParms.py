import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
from scipy import optimize
import sys
from os.path import expanduser
home = expanduser("~")
sys.path.insert(1, home+'/kode/simramesh/simramesh')
from geom_distr import distr 

data = pd.read_csv('metVertDistr.csv')
z = data['Points_2'].to_numpy()
xz = data['Point ID'].to_numpy().astype('float')
z /= z[-1] # Normalize levels
xz /= xz[-1] # Normalize argument


def f(x):
    return np.linalg.norm(z-distr(xz,x[0],x[1]))

res = optimize.minimize(f,x0=(0.8,5),bounds=((0,1),(1,None)))
x_g = res.x[0] 
alpha = res.x[1] 
print("Optimal values are x_g = %f and alpha = %f" % (x_g,alpha))
print("Linear transition at G(x_g) = %f" % distr(x_g,x_g,alpha))

npts = 1000
x = np.linspace(0,1,npts)
plt.plot(x,distr(x,1,alpha),label = 'Geometric')
plt.plot(x,x,label = 'Linear')
plt.plot(xz,z,label = 'Old simra algorithm')
plt.plot(x,distr(x,x_g,alpha),label = 'Combined geometric/linear')
plt.legend()
if False:
    plt.savefig('distr.pdf', dpi=400)
else:
    plt.show()


