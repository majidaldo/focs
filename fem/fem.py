import numpy as np


d={}
for atf in ['vx','vy','coords','phix','phiy']:
    fd=np.loadtxt(atf,dtype=[('m',int),('n',int),('v',float)])
    m=max(fd['m']+1)
    n=max(fd['n']+1)
    a=np.empty([m,n])
    for i in xrange(len(fd)):a[fd[i]['m']][fd[i]['n']]=fd[i]['v']
    d[atf]=a
    

    
# pl.quiver(d['coords'].T[0],d['coords'].T[1],d['phix'].T[0],d['phiy'][0])

#pl.scatter(d['coords'].T[0],d['coords'].T[1],marker='.')
