dfl=open('t.txt')
import numpy as np
from StringIO import StringIO
import matplotlib as mpl
from matplotlib import pylab

dfs=[]#data frames
ds=''
for al in dfl:
    if 'FRAME' in al:
        if ds!='':
            dfs.append(np.loadtxt(StringIO(ds)
            ,dtype=[('xi','uint'),('yi','uint'),('T','float')]))
        ds=''
        continue
    else:
        ds=ds+al
dfs.append(np.loadtxt(StringIO(ds)
            ,dtype=[('xi','uint'),('yi','uint'),('T','float')]))


def pia(adf):#put in  array
    adf.sort(order=['xi','yi'])
    nx=adf['xi'][-1]-adf['xi'][0]+1
    ny=adf['yi'][-1]-adf['yi'][0]+1 
    a=np.array(adf['T'],dtype='float').reshape(nx,ny)\
    .swapaxes(0,1)
    return a
#d=np.loadtxt('t.txt',dtype=[('xi','uint'),('yi','uint')
#,('T','float')])

#mpl.pyplot.contourf(pia(dfs[950]),cmap=mpl.cm.Greys)


