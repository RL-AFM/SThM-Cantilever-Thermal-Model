# -*- coding: utf-8 -*-
%matplotlib
"""
Created on Fri Feb 17 11:50:51 2017

@author: me
"""
import os.path
from gdsCAD import *
import numpy as np
import matplotlib.pyplot as plt
from math import log10, floor
import matplotlib
from matplotlib import style
from matplotlib import cm
style.use('ggplot')
from matplotlib import rc
rc('text', usetex=True)
matplotlib.rcParams['font.serif'] = 'CMU Serif'
matplotlib.rcParams['font.family'] = 'serif'

matplotlib.rcParams['font.size'] = 22
LARGE_FONT= ("Segoe UI Semilight", 12)
BUTTON_FONT = ("Segoe UI", 10)

#-------------------------------------------------------------------------------
def get_elements_from_GDS(myfile,plot=False):
	#start by opening the GDS file
	probe_GDS = core.GdsImport(myfile,layers={'SiN':0,'Au':2,'Pd':3})

	#get the list of cells
	cellname=probe_GDS.keys()
	#we only want one cell in the list - if there are more, break
	if len(cellname)>1: raise ValueError ('Too many cells in GDS file. Please include one cell only')
	#change the name from single item list to str
	cellname=cellname[0]

	if plot==True:drawProbeFromGDS(probe_GDS[cellname])
	#get the elements contained within this cell
	elements=probe_GDS[cellname].elements
	print '--> Sucessfully loaded elements from GDS file'
	return elements

#-------------------------------------------------------------------------------
#rounds an input number to the desired number of significant figures
round_to_n = lambda x, n: round(x, -int(floor(log10(x))) + (n - 1))

#-------------------------------------------------------------------------------
#check that all elements have appropriate dimensions (5,2) and sampling width
def check_dimensions_of_elements(set_of_elements):
	#use the first two x coordinates of the first point in the set to calculate a baseline dx, rounded to nearest sig fig
	dx=round_to_n(set_of_elements[0].points[1][0]-set_of_elements[0].points[0][0],1)
	#print dx
	for entry in set_of_elements:
		if entry.points.shape!=np.zeros((5,2)).shape:		#if not a 5,2 polygon then break
					print 'element has unusual shape'
					#raise ValueError ("Element coordinates not correct shape (5,2). Element is {0} at position{1}".format(entry,entry.points))
		elif round_to_n(entry.points[1][0]-entry.points[0][0],1)!=dx:	#if subsequent entries dont match our initial dx, raise this
			raise ValueError ('Some elements do not have equal spacing')
			continue
	print '--> Size check complete: all elements correct dimensions and have equal sampling width ({0}um)'.format(dx)
	dx*=1e-6
	print '----> dx = {0}'.format(dx)
	return dx

#-------------------------------------------------------------------------------
#puts elements different layers of the cell into seperate dictionary entries
def filter_elements_by_layer(set_of_elements):
	myDict={'SiN':[],'Au':[],'Pd':[],'Au_passive':[],'Pd_passive':[]}

	for entry in set_of_elements:
		if entry.layer==0: myDict['SiN'].append(entry)
		elif entry.layer==1: myDict['Au'].append(entry)
		elif entry.layer==2: myDict['Pd'].append(entry)
		elif entry.layer==3: myDict['Au_passive'].append(entry)
		elif entry.layer==4: myDict['Pd_passive'].append(entry)

		else:
			raise ValueError ("Layer outside acceptable range (0-4). Please use Layer 0 for SiN, 1 for Au, 2 for Pd, 3 for passive Au, 4 for passive Pd")
			break

	if len(myDict['SiN'])==0: raise ValueError ('No data in SiN layer')
	if len(myDict['Au'])==0: raise ValueError ('No data in Au layer')
	if len(myDict['Pd'])==0: raise ValueError ('No data in Pd layer')
	if len(myDict['Au_passive'])==0: print 'No data found for passive Au'
	if len(myDict['Pd_passive'])==0: print 'No data found for passive Pd'

	print '--> Successfully split layers into unique dict entries'
	return myDict

#-------------------------------------------------------------------------------
#helper function for sorting. Tells an element to look for x coordinate of its points
def getX0(element):
	return element.points[0][0]



def sort_in_x(dict_of_elements):
	dictSorted=newDictSameKeys(filtered)
	for name,item in dict_of_elements.iteritems():
		print 'Sorting {0}'.format(name)
		dictSorted[name]=sorted(item, key=getX0)
	print '--> Sucessfully sorted all materials by x position'
	return dictSorted

def newDictSameKeys(any_dict):
	return {key:[] for key in any_dict.keys()}
#-------------------------------------------------------------------------------

def simpleWidth(oneMaterial):
	#prep a new array for the result of width calculation
	width=np.zeros(len(oneMaterial))
	for count, entry in enumerate(oneMaterial):
		#the width of any element is its area divided by sampling width. x2 for symmetry.
		width[count]=((entry.area()*1e-12)/dx)*2.
	return width


def getWidths(dict_of_elements):
	#get the length of each material's array
	lengths = {key:len(value) for key,value in dict_of_elements.iteritems()}

	#if one of them is zero, you dont have any data there
	if any(x == 0 for x in lengths):
		raise ValueError ('No data in one of the arrays - one material has no data on its layer')
	#construct a dictionary of widths using helper function simpleWidth
	widths = {key:simpleWidth(value) for key,value in dict_of_elements.iteritems()}
	print '--> Successfully calculated widths of all elements for all materials'
	return widths

#-------------------------------------------------------------------------------
#prints out a dictionary with keys and micron scaled values
def micronPrint(anyDict):
	for key,value in anyDict.iteritems():
		print key,["{0:0.1f}".format(i*1e6) for i in value]

#-------------------------------------------------------------------------------

#Zero Pad
def zeroPad(width_dict):
	#dictionary comprehemsion. make a new dict with same keys containing length of arrays
	lengths = {key:len(value) for key,value in width_dict.iteritems()}
	#get the difference in length between the longest array (SiN) and the others
	differences={key:(lengths['SiN']-value) for key, value in lengths.iteritems()}

	#we can't pad with a dict comprehension, since Pd and Au need different pad types (front and rear)
	width_dict['Au']=np.pad(width_dict['Au'],(0,differences['Au']),'constant')
	width_dict['Au_passive']=np.pad(width_dict['Au_passive'],(0,differences['Au_passive']),'constant')
	width_dict['Pd']=np.pad(width_dict['Pd'],(differences['Pd'],0),'constant')
	width_dict['Pd_passive']=np.pad(width_dict['Pd_passive'],(differences['Pd_passive'],0),'constant')
	print '-->Successfully zero padded metal arrays. All arrays equal length'
	return width_dict


#-------------------------------------------------------------------------------

def plotArray(arr,log=False):
	xarr=np.linspace(0,len(arr)*dx,len(arr))
	plt.plot(xarr*1e6,arr,'.')
	plt.xlabel('Length ($\\mu m$)')
	plt.ylabel('Width ($\\mu m$)')
	if log==True:plt.yscale('log')
	plt.show()

#-------------------------------------------------------------------------------

def getElementalResistance(padded_dict):
	r={}	#resistances
	g={}	#conductances
	n=len(padded_dict['SiN'])

	#working thicknesses and thermal resistances
	k={'SiN':3.,'Au':153.,'Pd':31.,'Au_passive':153.,'Pd_passive':31.}
	t={'SiN':400e-9,'Au':150e-9,'Pd':40e-9, 'Au_passive':150e-9, 'Pd_passive':40e-9}
	print '***Calculating elemental thermal resistances. Have you checked Pd/Pt???'

	#setup blank resistance dict
	for keys, vals in padded_dict.iteritems():
		r[keys]=np.zeros(n)
		#calculate the elemental resistance of each material, ignoring placeholder zeros
		for count,entry in enumerate(vals):
			if entry==0.: continue
			#if on the flat we good
			elif count*dx<=139e-6:
				r[keys][count]=dx/(entry*k[keys]*t[keys])
			#on the pyramid dx is elongated and the thickness of metal layers is reduced
			#SiN thickness is the same though, which is a nuisance. We have to catch this case
			else:
				if keys=='SiN': r[keys][count]=tipTilt(dx)/(entry*k[keys]*t[keys])
				else: r[keys][count]=tipTilt(dx)/(entry*k[keys]*correctEvapOnPyramid(t[keys]))



		#print r[keys]

	for keys, vals in r.iteritems():
		g[keys]=reciprocal(vals)
	r['Parallel']=np.zeros(n)
	r['Parallel']=(np.sum((g['SiN'],g['Au'],g['Pd'],g['Au_passive'],g['Pd_passive']),axis=0))**-1

	return r
#-------------------------------------------------------------------------------
#get the reciprocal of any nonzero elements
def reciprocal(array):
	output=np.zeros(len(array))
	for n in range(0,len(output)):
		if array[n]!=0:
			output[n]=array[n]**-1
	return output

#-------------------------------------------------------------------------------
#plots all elements of a given dictionary in log or linspace y-axis. x axis is always length in um
def plotDict(anyDict,log=False):
	xarr=np.linspace(0,len(anyDict['SiN'])*dx,len(anyDict['SiN']))
	try:
		for key, entry in anyDict.iteritems():
			plt.plot(xarr*1e6, entry,label=key,c=colour_dict[key])
	except KeyError:
		pass

	if log==True:plt.yscale('log')
	plt.legend()
	plt.xlabel('Length ($\\mu m$)')
	plt.show()

#the width of the platinum sensor is not the same in x as it is for the direciton of current flow
#this can make a big difference in the power calculation, so must be corrected for.
#we also have the concern of elongation and reduced thickness up the pyramid end,
#but perhaps that can wait till later. We dont necessarily want to do this correction
#in place, rather it would be better if this makes a new array. We can put into
#a dict later if required
def transform_heater(Pd_width):
	output=Pd_width*np.sin(45)
	print '--> Transforming width of heater for power calculation'
	#print '---->Width of platinum in direction of current flow = {0}m'.format(output[-1])
	return output

#calculate the generated thermopower of the metal elements. Pd needs a correction for the
#45 degree angle (its narrower in x' (direction of current flow) than in x).
#this is already handled in constructMetalDicts
def getJouleHeat(metal_dict):
	n=len(metal_dict['Pd']['width'])
	I=2.2e-3
	for material, dictionary in metal_dict.iteritems():
		metal_dict[material]['power']=np.zeros(n)
		for s in range(0,n):
			if dictionary['width'][s]!=0:
				w=dictionary['width'][s]/2.
				#if before the pyramid region, L=dx as normal
				if s*dx<=139:
					dictionary['power'][s]=2.*((I**2)*dictionary['resistivity']*dx)/(dictionary['thickness'][s]*w)
				#when on the pyramid, the drawn shape is actually dx/cos46.5 degrees longer
				#and the deposited metal layers are thinner
				#reduced thickness of metal arrays is already handled in constructMetalDicts()
				else:
					dictionary['power'][s]=2.*((I**2)*dictionary['resistivity']*tipTilt(dx))/(dictionary['thickness'][s]*w)
	print '--> Successfully calculated elemental power generation for all metal layers'



#creates a nested dict of metal layer parameters. Thickeness, width and resistivity as their own entires
def constructMetalDicts(padded_dict):
	Pd={}
	Au={}
	Metals={'Au':Au,'Pd':Pd}
	for key,value in padded_dict.iteritems():
		if key=='Pd':
			Pd['width']=transform_heater(value)	#correct width in direction of current flow for 45 degree angle
			Pd['resistivity']=0.238e-6
			#thickness = 40nm anywhere there is non-zero width
			Pd['thickness']= [40e-9 if x!=0. else 0. for x in iter(value)]
			#now reduce the thickness at any point that was deposited on the pyramid
			Pd['thickness']= [x if (dx*count)<=139e-6 else correctEvapOnPyramid(x) for count,x in enumerate(Pd['thickness'])]

		elif key =='Au':
			Au['width']=value
			Au['resistivity']=47.847e-9
			Au['thickness']= [150e-9 if x!=0. else 0. for x in iter(value)]
			Au['thickness']= [x if (dx*count)<=139e-6 else correctEvapOnPyramid(x) for count,x in enumerate(Au['thickness'])]
	return Metals

#takes the generated power from Au and Pd and sums them (needed for getTemp calc)
def powerSum(metal_dict):
	return np.sum((metal_dict['Pd']['power'],metal_dict['Au']['power']),axis=0)

#no air in this version
def getTemp(parallel_resistance,summed_power,sample_resistance=np.inf):
	n=parallel_resistance.size
	if summed_power.size!=n:
		raise ValueError ('Power and resistance arrays are not equal size')

	if sample_resistance!=(np.inf): print '-->Calculating probe contact temperature\n---->Sample Rth={0}'.format(sample_resistance)
	else: print '-->Calculating probe out of contact temperature'
	#we actually want thermal conductance, not resistance
	g=reciprocal(parallel_resistance)

	#construct the conductance matrix
	A=np.zeros((n,n))

	#fill in the matrix
	for r in range (0,n-1):
		A[r,r]=g[r]+g[r+1] #diagonals
		A[r,r+1]=A[r+1,r]=g[r+1]*-1 #either side of diagonal
		A[n-1,n-1]=g[n-1]+(sample_resistance**-1) #final diagonal position (n,n)

	Ainv=np.linalg.inv(A)
	V=np.dot(Ainv,summed_power)

	return V

#****************************
#value from Kamil's Thesis
def tipTilt(num):
	return num/(np.cos(np.deg2rad(46.5)))

def correctEvapOnPyramid(thickness):
	return thickness*np.cos(np.deg2rad(43.5))

#******************************
#draws the full probe from the half probe cell provided by GDS
def drawProbeFromGDS(cell):
	cell.show()		#show initial cell
	#reflect the elements of the cell along x axis () - for some reason this doesnt work
	#on full cells, only boundaries.
	other_half=[i.reflect('x') for i in cell.elements]
	for b in other_half: b.show()	#plot all the flipped elements





colour_dict={'SiN':'green','Au':'gold','Pd':'grey','Parallel':'red'}
#load the file
elements=get_elements_from_GDS('4Tdesign1_f500nm.gds')
#check shapes and dimensions
dx=check_dimensions_of_elements(elements)
#filter by layer type
filtered=filter_elements_by_layer(elements)
#sort the layers by x position
filteredAndSorted=sort_in_x(filtered)

#)()()()()()()()()()()()()()()()()()()()()()()()()()()()()()

def collapseStackedElements(sorted_dict):
	output=newDictSameKeys(sorted_dict)
	for key,value in sorted_dict.iteritems():
		print 'calculating {0}'.format(key)
		x_pos=[]
		for n in value:
			x_pos.append(getX0(n))
		x_pos=np.unique(x_pos)
		#print len(x_pos)
		#"""we are not getting unique values here"""
		#x_pos=np.unique(getX0(value))	#the x array should be the number of unique positions of elements in the array
		print '{0} unique elements found'.format(len(x_pos))
		list2d=[[] for _ in range(len(x_pos))]
		#iterate over the x positions array once inserting the Boundaries present at the locations in the relevant sample of the 2dlist
		for counter,number in enumerate(x_pos):
			for element in value:
				if getX0(element)==number:
					list2d[counter].append(element)

		#iterate over it again performing the calculation
		for n in range(len((x_pos))):
			placeholder=[]
			for entry in list2d[n]:
				if type(entry)!=np.float32:
					placeholder.append(entry.area())
			output[key].append(2*(((sum(placeholder))*1e-12)/dx))

	return output







out=collapseStackedElements(filteredAndSorted)

for n,v in out.iteritems():
	print n,len(v)

#get widths
#widths=getWidths(filteredAndSorted)
#do zero padding
padded=zeroPad(out)
#calculate elemental resistances
thermal_resistances=getElementalResistance(padded)
#plot the resistance graph to ensure it is correct
plotDict(thermal_resistances,True)
#before getting the Joule heating, prepare the metal layers
Metal=constructMetalDicts(padded)



#get the Joule heating power. Just append it to the metals' dictionary under key 'power'
getJouleHeat(Metal)
#plot the powers for each layer to make sure
for key,val in Metal.iteritems():
	plt.legend(key)
	plotArray(val['power'],True)

#just checking the power sum worked
plotArray(powerSum(Metal),True)
#get the temperature rise given the selected current (hardcoded)
T=getTemp(thermal_resistances['Parallel'],powerSum(Metal))

plotArray(T)
