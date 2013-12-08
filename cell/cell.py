import data
from data import (conn,coords,cells)

import numpy as np
import matplotlib.pyplot as plt
#import mayavi
#from mayavi import tools
#from mayavi.tools import pipeline
#tm=pipeline.triangular_mesh_source(x, y, z, triangles, **kwargs)
#tm(coords.T[0],coords.T[1],coords.T[2],conn)
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.colors import colorConverter
cc = lambda arg: colorConverter.to_rgba(arg, alpha=0.6)

ax=Axes3D(plt.figure())

# x = [0,1,1,4]
# y = [0,0,1,1]
# z = [0,1,0,3]
# verts = [zip(x, y,z)]


verts=[   ((coords[atri[0]][0],coords[atri[0]][1],coords[atri[0]][2],)
           ,(coords[atri[1]][0],coords[atri[1]][1],coords[atri[1]][2],)
           ,(coords[atri[2]][0],coords[atri[2]][1],coords[atri[2]][2],)) for atri in conn[:10]  ]

#verts=[[[1,2,3],[4,3,2],[3,3,5]]]
p3dcol=Poly3DCollection(verts,facecolors= [cc('r')] )
ax.add_collection3d(p3dcol)
ax.autoscale()

p3dcol.set_linewidth(0)#dont want edge




