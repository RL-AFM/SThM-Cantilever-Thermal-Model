# -*- coding: utf-8 -*-
from gdsCAD import *
from Tkinter import Tk
from tkFileDialog import askopenfilename
from math import log10, floor
import numpy as np
#------------------------FILE CHOOSER-------------------------------------------
def choosefile():
    Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
    filename = askopenfilename() # show an "Open" dialog box and return the path to the selected file
    print(filename)
    return filename
#----------------------------------------------
def newDictSameKeys(any_dict):
    return {key:[] for key in any_dict.keys()}

#------Given a GDSII file with one cell, return all the elements present--------
def get_elements_from_GDS(myfile,plot=False):
    #start by opening the GDS file
    print '--> Importing data from file: {0}'.format(myfile)
    probe_GDS = core.GdsImport(myfile)

    #get the list of cells
    cellname=probe_GDS.keys()
    #we only want one cell in the list - if there are more, break
    if len(cellname)>1: raise ValueError ('Too many cells in GDS file. Please include one cell only')
    #change the name from single item list to str
    cellname=cellname[0]

    #get the elements contained within this cell
    elements=probe_GDS[cellname].elements
    print '--> Sucessfully loaded elements from GDS file'
    return elements

#-------------------------------------------------------------------------------
#puts elements different layers of the cell into seperate dictionary entries
def filter_elements_by_layer(set_of_elements):
    myDict={'SiN':[],'Au':[],'Pd':[],'Au passive':[],'Pd passive':[]}

    for entry in set_of_elements:
        if entry.layer==0: myDict['SiN'].append(entry)
        elif entry.layer==1: myDict['Au'].append(entry)
        elif entry.layer==2: myDict['Pd'].append(entry)
        elif entry.layer==3: myDict['Au passive'].append(entry)
        elif entry.layer==4: myDict['Pd passive'].append(entry)

        else:
            raise ValueError ("Layer outside acceptable range (0-4). Please use Layer 0 for SiN, 1 for Au, 2 for Pd, 3 for passive Au, 4 for passive Pd")
            break

    if len(myDict['SiN'])==0: raise ValueError ('No data in SiN layer')
    if len(myDict['Au'])==0: raise ValueError ('No data in Au layer')
    if len(myDict['Pd'])==0: raise ValueError ('No data in Pd layer')
    if len(myDict['Au passive'])==0: print 'No data found for passive Au'
    if len(myDict['Pd passive'])==0: print 'No data found for passive Pd'

    print '--> Successfully split layers into unique dict entries'
    return myDict

#-------------------------------------------------------------------------------
#helper function for sorting. Tells an element to look for x coordinate of its points
def getX0(element):
    return element.points[0][0]

def sort_in_x(dict_of_elements):
    print '--> Sorting elements by their x position'
    dictSorted=newDictSameKeys(dict_of_elements)
    for name,item in dict_of_elements.iteritems():
        print '---->Sorting {0}'.format(name)
        dictSorted[name]=sorted(item, key=getX0)
    print '--> Sucessfully sorted all materials by x position'
    return dictSorted


#-------------------------------------------------------------------------------
#check that all elements have appropriate dimensions (5,2), i.e a 4 point polygon
def check_dimensions_of_elements(myDict):
    for name,data in myDict.iteritems():
        print '--> Checking element shape of {0} layer'.format(name)
        for entry in data:
            if entry.points.shape!=np.zeros((5,2)).shape:        #if not a 5,2 polygon then break
                print '----> Element coordinates not correct shape (5,2). Element is {0} at position:\n{1}'.format(entry,entry.points)
    print '--> Element checking complete\n'


#-------------------------------------------------------------------------------
#rounds an input number to the desired number of significant figures
round_to_n = lambda x, n: round(x, -int(floor(log10(x))) + (n - 1))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#DEPRECATED
#-------------------------------------------------------------------------------
#given an element object, find its horizontal length (dx)
#supply with [0] and [-1] to find the first and last dx's
def getdx(segment):
    return round_to_n(segment.points[1][0]-segment.points[0][0],2)


def get_dual_dx(SiN_ordered):
    lowdx=getdx(SiN_ordered[0])*1e-6
    hidx=getdx(SiN_ordered[-1])*1e-6
    print '--> Low resolution dx = {0:.2}um'.format(lowdx*1e6)
    print '--> High resolution dx = {0:.2}um'.format(hidx*1e6)
    return lowdx, hidx
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#DEPRECATED
#-------------------------------------------------------------------------------
