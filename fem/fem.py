import numpy as np


d={}
for atf in set(['vx','vy','coords','phi','p'
           # ,'d','mc' diagnostics
            #,'ml','mlm1','mc','rx','ry','d'x
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
p=d['p'].T[0]

import matplotlib
from matplotlib.mlab import griddata
from matplotlib import pyplot as plt
# griddata and contour.
xi = np.linspace(min(x),max(x),1e3)
yi = np.linspace(min(y),max(y),1e3)

phizi = griddata(x,y,phi,xi,yi,interp='linear')
pzi = griddata(x,y,p,xi,yi,interp='linear')
circm=np.zeros_like(phizi,dtype=bool)


from numba import autojit
@autojit
def maskphi():
    for ix,xx  in enumerate(xi):
        for iy,yy in enumerate(yi):
            if(xx**2+yy**2<.3**2):phizi[ix,iy]=True
maskphi();
                
phima=np.ma.masked_array(phizi,mask=circm)

def pltphi():
    return plt.contourf(xi,yi,phima
                        ,25
                 ,cmap=plt.cm.bone
             #,Normalize=plt.normalize(vmax=abs(zi).max(), vmin=-abs(zi).max())
             )
   
def pltv():#not used
    return plt.quiver(x,y,vx,vy,(vx**2+vy**2)**.25#pzi
                      ,cmap=plt.cm.hot#RdBu
            )

def pltpmm():
    with open('pmax') as f: ipmax=int(f.read() )
    with open('pmin') as f: ipmin=int(f.read())
    for (al,ax,ay,clr) in [ ('min p='+str((p[ipmin])),x[ipmin],y[ipmin],'b')
                           ,('max p='+str(p[ipmax]),x[ipmax],y[ipmax],'r')]:
        plt.scatter([ax],[ay]
                       ,s=40#,markersize=20
                       ,color=clr#
                       ,label=al
                       )
#    plt.legend(lns,['min'+str(p[ipmin]),'max'+str(p[ipmax])])


def pltp():
    return plt.contourf(xi,yi,pzi,50
                 )


def pltall():
    pltphi();pltv();pltpmm();
    plt.xlim(min(x),max(x))
    plt.ylim(min(y),max(y))
    plt.gcf().set_size_inches(3,3)
    plt.legend();
