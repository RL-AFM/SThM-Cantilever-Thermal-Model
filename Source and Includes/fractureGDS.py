%matplotlib

import gdspy
import numpy as np
import gdsmanipulation as gdsm

"""Load a gds file and fracture into 1um/25nm segments"""

'''Uncomment below to chose the file via GUI'''
filename=gdsm.choosefile()
'''Else hardcode it'''
#folderpath='C:/Users/me/OneDrive - University of Glasgow/PhD/2018/Cantilever-Modeller-GDS-Edition-LINUX FORK/gds generated for testing/'
#filename=folderpath + 'Dual dx tested/Original and quick cut/OriginalCant_unfractured.gds'

theDir='/'.join(filename.split('/')[:-1])+'/'
theFile=''.join(filename.split('/')[-1][:-4])+'_Fractured.gds'

gdsii = gdspy.GdsLibrary()

#read in the gds file to an object named layout
layout=gdsii.read_gds(filename)

#ensure the user only inputs files with one cell
if len(layout.cell_dict)>1: raise ImportError('Too many cells in input file')

#extract the cell by its name
#get the 0th key in the dictionary of cells, which = the name of the cell as a string
c=gdsii.extract(layout.cell_dict.keys()[0])

#filter the polygons by their layer. Result is a dictionary of form {Material:[polygons]}
#probably fails if one layer cotains more than one polygon
filtered=gdsm.filter_elements_by_layer(c.elements)

#define an output cell
cut=gdspy.Cell('Fractured')
#define an x array of 25nm sampling width for the tip. Omit the last sample, where width would = 0
tipArray=np.linspace(139., 149.950, 439)


layernum={'SiN':0,'Au':1,'Pd':2,'Au passive':3,'Pd passive':4}

#for all materials, perform a slice operation every 1um along the polygon.
#Three polygons result per slice operation, (left, slice, right), we only care
#about the slice (element[1])
for key, val in filtered.iteritems():
    for n in range(0,139):
        cut.add(gdspy.slice(val, [n,n+1], 0, layer=layernum[key])[1])


#as above but with 25nm segments
for key, val in filtered.iteritems():
    for n in tipArray:
        cut.add(gdspy.slice(val, [n,n+0.025], 0, layer=layernum[key])[1])

OutputCell=gdspy.GdsLibrary()
OutputCell.add(cut)
OutputCell.write_gds(theDir+theFile)
