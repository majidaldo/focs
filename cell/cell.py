import data
from data import (conn,coords,cellframes,nt2clr,nframes,ncells)


import numpy as np
import matplotlib.pyplot as plt




from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.colors import colorConverter
cc = lambda arg: colorConverter.to_rgba(arg, alpha=0.6)#color2rgba
ct=lambda nt: cc(nt2clr[nt])#celltype2color



# ax=Axes3D(plt.figure())

# verts=[   ((coords[atri[0]][0],coords[atri[0]][1],coords[atri[0]][2],)
#            ,(coords[atri[1]][0],coords[atri[1]][1],coords[atri[1]][2],)
#            ,(coords[atri[2]][0],coords[atri[2]][1],coords[atri[2]][2],)) for atri in conn[:]  ]


# fc= lambda : [ct(act) for act in cellframes[0][:len(verts)]] 

# #verts=[[[1,2,3],[4,3,2],[3,3,5]]]
# p3dcol=Poly3DCollection(verts                             #frame0
#                         ,facecolors= fc() )
# ax.add_collection3d(p3dcol)
# ax.set_xlim3d(min(coords.T[0]),max(coords.T[0]))
# ax.set_ylim3d(min(coords.T[1]),max(coords.T[1]))
# ax.set_zlim3d(min(coords.T[2]),max(coords.T[2]))

# p3dcol.set_linewidth(0)#dont want edge
# ax.set_frame_on(False)
# ax.get_xaxis().set_visible(False)#doesn't work with axes3d apparently
# ax.get_yaxis().set_visible(False)



#plt.show

# import mayavi
# #from mayavi import tools
# from mayavi.tools import pipeline
# tm=pipeline.triangular_mesh_source(    coords.T[0],coords.T[1],coords.T[2],conn
#                                     , cell_scalars=cellframes[100]
#                                     ,representation='wireframe'
#                                     ,opacity=0
#                                     )
# # #tm(coords.T[0],coords.T[1],coords.T[2],conn)

from mayavi import mlab
# Plot it
cellframes=np.array(cellframes,dtype=float)


mlab.clf()
mesh = mlab.triangular_mesh(coords.T[0],coords.T[1],coords.T[2],conn
                            ,representation='wireframe'
                            ,opacity=0
                            )

cell_data = mesh.mlab_source.dataset.cell_data
cell_data.scalars = cellframes[0]
cell_data.scalars.name = 'Cell data'
cell_data.update()

mesh2 = mlab.pipeline.set_active_attribute(mesh,
                cell_scalars='Cell data')
cs=mlab.pipeline.surface(mesh2)

#execfile this file in mayavi
#take a moment to config cell data "clrs and lgnds"legend w the gui then anim

def an(arange):# inefficient but idk another way
    for af in arange:

        cell_data = mesh.mlab_source.dataset.cell_data
        cell_data.scalars = cellframes[af]
        cell_data.scalars.name = 'Cell data'
        cell_data.update()

        mesh2 = mlab.pipeline.set_active_attribute(mesh,
            cell_scalars='Cell data')
        cs=mlab.pipeline.surface(mesh2)
        mlab.show()
        mlab.savefig(str(af)+'face.jpg')
        #mlab.clf(mlab.gcf())
        #mlab.close()
        mesh2.remove()
        cs.remove()
        del mesh2, cs, cell_data
        #mesh.remove()
        #mlab.close()
#mlab.show()

#mlab.close() everythong closes
