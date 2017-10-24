# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 11:50:51 2017

@author: me
"""
import os.path
from gdsCAD import *
import numpy as np
import matplotlib.pyplot as plt
from math import log10, floor


# Create some things to draw:
#dont know how to get the layers dict to work properly
probe_GDS = core.GdsImport('fractured10nm.gds',layers={'SiN':0,'Au':2,'Pd':3})
#assign the bundle of elements to a variable so we are not constantly typing...

#get the name of the cell being used
probe_GDS.keys()[0]

elements=probe_GDS['fineFracture10nm'].elements


#get the number of elements in the file.
len(elements)

#sorting the group of elements by their x coordinate position (I suppose this must be done by material).

#helper function for sorting. Tells an element to look for x coordinate of its points
def getX0(element):
	return element.points[0][0]

sorted_elements=sorted(elements, key=getX0)

for entry in sorted_elements:
	print entry.points

#in the layout file probeGDS, in cell name fractures, in the first element, tell me the x coordinate of the first pair
getX0(elements[0])

#should be called on the sorted group of elements

def getAreas(set_of_elements):
	areas=np.zeros(len(set_of_elements))
	if areas.size == 0:
		raise ValueError ('No data in one of the arrays - one material has no data on its layer')
	for count,entry in enumerate(set_of_elements):
		areas[count]=entry.area()*1e-12				#scaling factor! um*um = 1e-12
	return areas

myAreas=getAreas(sorted_elements)
widths=myAreas/1e-6
print widths





def filter_elements_by_layer(set_of_elements):
	lst1=[]
	lst2=[]
	lst3=[]
	for entry in set_of_elements:
		if entry.layer==0: lst1.append(entry)
		elif entry.layer==1: lst2.append(entry)
		elif entry.layer==2: lst3.append(entry)
		else:
			raise ValueError ("Layer outside acceptable range (0-2). Please use Layer 0 for SiN, 1 for Au, 2 for Pd")
			break
	return {'SiN':lst1,'Au':lst2,'Pd':lst3}

filtered=filter_elements_by_layer(sorted_elements)
filtered['SiN'][0].points



widths={'SiN':None,'Au':None,'Pd':None}

#for all of the materials, sort the elements by x and find the widths
for material, value in filtered.items():
	print material
	widths[material]=sorted(value,key=getX0)
	widths[material]=getAreas(widths[material])

round(0.12845, 2)

#rounds an input number to the desired number of significant figures
round_to_n = lambda x, n: round(x, -int(floor(log10(x))) + (n - 1))

#check that all elements have appropriate dimensions (5,2) and sampling width
def check_dimensions_of_elements(set_of_elements):
	#use the first two x coordinates of the first point in the set to calculate a baseline dx, rounded to nearest sig fig
	dx=round_to_n(set_of_elements[0].points[1][0]-set_of_elements[0].points[0][0],1)
	print dx
	for entry in set_of_elements:
		if entry.points.shape!=np.zeros((5,2)).shape:		#if not a 5,2 polygon then break
					raise ValueError ("Element coordinates not correct shape (5,2). Element is {0} at position{1}".format(entry,entry.points))
		elif round_to_n(entry.points[1][0]-entry.points[0][0],1)!=dx:	#if subsequent entries dont match our initial dx, raise this
			raise ValueError ('Some elements do not have equal spacing')
			continue
	print 'Size check complete: all elements correct dimensions and have  sampling width ({0})'.format(dx)
	return dx

#@10nm dx is reported as 0.01 therefore we can assume um units
dx=check_dimensions_of_elements(elements)





"""

#get the points of the boundary in the cell in the file.
def get_cantilever_points():
	pointList=[]
	for element in probe_GDS['root']:
		pointList.append(element.points)


	return pointList#probe_GDS['root'][0].points

verticeList=get_cantilever_points()

print verticeList

#SiN=np.rec.array(verticeList[0], dtype=['x0','w-'),('knee','w-'),('tip','y0'),('knee','w+'),('x0','w+'),('x0','w-')})
SiNx=np.array(verticeList[0][:,0])
SiNx=np.rec.array(SiNx, dtype=[('p0','f4'),('knee','f4'),('tip','f4'),('p3','f4'),('p4','f4'),('p5','f4')])
print SiNx.knee
print SiNx.tip

SiNy=np.array(verticeList[0][:,1])
SiNy=np.rec.array(SiNy, dtype=[('wn','f4'),('p1','f4'),('tip','f4'),('wp','f4'),('p4','f4'),('p5','f4')])

print SiNy.wn
print SiNy.wp


global dx
dx=0.1
result = probe_GDS.slice()



def samples(value):
    return (value/dx)+1


def setup(SiN):
    end=SiN[2,0]
    x = np.linspace(0,end,samples(end),endpoint=True)
    y=np.zeros(len(x))
    return x,y

x,y=setup(verticeList[0])

verticeList=np.array(verticeList)
verticeList_samples=samples(verticeList)
m=SiNy.wp/(SiNx.tip-SiNx.knee)
c=-1*m*SiNx.tip
y[0:verticeList_samples[0][1][0]]=verticeList[0][3][1]
y[verticeList_samples[0][1][0]:]=m*x[verticeList_samples[0][1][0]:]+c
plt.plot(x,y)
plt.show()

#furthest_point(verticeList[0],0.1)
#print np.linspace(0,150,151,endpoint=True)
"""
