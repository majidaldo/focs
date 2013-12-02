import numpy as np


d={}
for atf in set(['vx','vy','coords','phi'
           # ,'d','mc' diagnostics
            #,'ml','mlm1','mc','rx','ry','d'
            ]):
    fd=np.loadtxt(atf,dtype=[('m',int),('n',int),('v',float)])
    m=max(fd['m']+1)
    n=max(fd['n']+1)
    a=np.empty([m,n])
    for i in xrange(len(fd)):a[fd[i]['m']][fd[i]['n']]=fd[i]['v']
    d[atf]=a


x = d['coords'].T[0]
y = d['coords'].T[1]
phi = d['phi'].T[0]
vx=d['vx'].T[0]
vy=d['vy'].T[0]

import matplotlib
from matplotlib.mlab import griddata
from matplotlib import pyplot as plt
# griddata and contour.
xi = np.linspace(min(x),max(x),100)
yi = np.linspace(min(y),max(y),100)
zi = griddata(x,y,phi,xi,yi,interp='linear')

def pltphi():
    plt.contour(xi,yi,zi,15
                 ,cmap=plt.cm.gray
             #,Normalize=plt.normalize(vmax=abs(zi).max(), vmin=-abs(zi).max())
             )

    
def pltv():
    plt.quiver(x,y,vx,vy)
