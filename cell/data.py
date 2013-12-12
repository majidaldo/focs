
import numpy as np

mui=np.uint32
mf=np.float16
fd=np.loadtxt("coords",dtype=[('r',mui),('c',mui),('v',mf)])
m=max(fd['r']+1)
n=max(fd['c']+1)
coords=np.empty([m,n])
for i in xrange(len(fd)):coords[fd[i]['r']][fd[i]['c']]=fd[i]['v']
fd=  np.loadtxt("conn",  dtype=[('r',mui),('c',mui),('v',mui)])
m=max(fd['r']+1)
n=max(fd['c']+1)
conn=np.empty([m,n])
for i in xrange(len(fd)):conn[fd[i]['r']][fd[i]['c']]=fd[i]['v']
del m;del n;del fd

cellframes= np.loadtxt("cells", dtype=np.uint8)

ct2num={'normal':0,'cancer':1,'complex':2,'necrotic':3}
nt2clr={0:'white',1:'yellow',2:'green',3:'black'}

ncells=len(conn)
nframes=len(cellframes)

