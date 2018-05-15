# -*- coding: utf-8 -*-
%matplotlib

"""
Created on Fri Feb 17 11:50:51 2017

@author: me
"""

# -*- coding: utf-8 -*-
import os.path
import datetime

import numpy as np
import numpy.polynomial.polynomial as poly

import matplotlib.pyplot as plt
import matplotlib
from matplotlib import style
from matplotlib import cm
from matplotlib import rc
style.use('ggplot')

import matplotlib.ticker
rc('text', usetex=True)
rc('figure', figsize=(22,14.2))
matplotlib.rcParams['font.serif'] = 'CMU Serif'
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.size'] = 22
import tabulate

import fitsig2 as fs2
from gdsmanipulation import *
import add_power

lowdx=1e-6
hidx=25e-9

fig_of_merit={}

plt.ioff()
f,ax = plt.subplots(2,2)

fdict={}

img='Schematic'
r='Resistance $K/W$'
p='Power $W$'
t='Temperature rise $K$'

fdict[img]=ax[0][0]
fdict[r]=ax[0][1]
fdict[p]=ax[1][0]
fdict[t]=ax[1][1]

for key,plot in fdict.iteritems():
    plot.set_xlabel('Length ($\mu m$) ')

colour_dict={'SiN':'darkgreen','Au':'orange','Pd':'darkgrey','Au passive':'gold','Pd passive':'lightblue','Parallel':'red'}




#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------

def plotArray(arr,plotref,xax,log=False):
    ax=fdict[plotref]
    ax.plot(xax,arr)
    ax.set_xlabel('Length ($\\mu m$)')
    ax.set_ylabel(plotref)
    if log==True:ax.set_yscale('log')

def plotPower(arr,plotref,xax,log=False):
    ax=fdict[plotref]
    scaled=[arr[n] if xax[n]<=139 else arr[n]*(lowdx/hidx) for n,v in enumerate(xax)]
    ax.plot(xax,scaled)
    ax.set_xlabel('Length ($\\mu m$)')
    ax.set_ylabel(plotref)
    ax.text(0,1e-5, 'Total power generated = {0:.2f}mW'.format(np.sum(arr)*1e3))
    if log==True:ax.set_yscale('log')



#-------------------------------------------------------------------------------

def getElementalResistance(padded_dict):
    r={}    #resistances
    g={}    #conductances
    n=padded_dict['SiN'].shape[0]

    #working thicknesses and thermal resistances
    k={'SiN':3.,'Au':153.,'Pd':31.,'Au passive':153.,'Pd passive':31.}
    t={'SiN':400e-9,'Au':150e-9,'Pd':40e-9, 'Au passive':150e-9, 'Pd passive':40e-9}
    print '***Calculating elemental thermal resistances. Have you checked Pd/Pt???'

    #setup blank resistance dict
    for keys, vals in padded_dict.iteritems():
        r[keys]=np.zeros(n)
        #calculate the elemental resistance of each material, ignoring placeholder zeros
        for count,(xpos,w) in enumerate(vals):
            #if no width, there is no element present, ignore
            if w==0.: continue
            #if its on the flat
            elif xpos<=139:
                r[keys][count]=lowdx/(w*k[keys]*t[keys])
            else:
                if keys=='SiN': r[keys][count]=tipTilt(hidx)/(w*k[keys]*t[keys])
                else: r[keys][count]=tipTilt(hidx)/(w*k[keys]*correctEvapOnPyramid(t[keys]))

        #print r[keys]

    for keys, vals in r.iteritems():
        g[keys]=reciprocal(vals)
    r['Parallel']=np.zeros(n)
    r['Parallel']=(np.sum((g['SiN'],g['Au'],g['Pd'],g['Au passive'],g['Pd passive']),axis=0))**-1
    reportRth(r)

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
def plotDict(anyDict,plotref,xax,log=False):

    ax=fdict[plotref]
    for key, entry in anyDict.iteritems():
        scaled=[v if xax[n]<=139 else v*(lowdx/hidx) for n,v in enumerate(entry)]
        try:
            ax.plot(xax, scaled,label=key,c=colour_dict[key])
        except KeyError:
            print 'there was a key error'
            continue
    ax.set_ylabel(plotref)
    ax.text(0,1e6, 'Total thermal resistance={0:.2E}K/W'.format(np.sum(thermal_resistances['Parallel'])))
    fig_of_merit['Parallel Thermal Resistance (K/W)'] = ['{0:.2E}'.format(np.sum(thermal_resistances['Parallel']))]
    if log==True:ax.set_yscale('log')
    f.legend(ax.lines, colour_dict.keys(),'upper left')



#the width of the platinum sensor is not the same in x as it is for the direciton of current flow
#this can make a big difference in the power calculation, so must be corrected for.
#we also have the concern of elongation and reduced thickness up the pyramid end,
#but perhaps that can wait till later. We dont necessarily want to do this correction
#in place, rather it would be better if this makes a new array. We can put into
#a dict later if required
def transform_heater(Pd_width):
    output=Pd_width*np.sin(np.deg2rad(45))
    print '--> Transforming width of heater for power calculation'
    #try: print '----> [array] Width of Palladium in direction of current flow = {0}m'.format(output)
    #except IndexError: print'----> [scalar] Width of Palladium in direction of current flow = {0}um'.format(output*1e6)
    return output

#ask the user for a current value
def askCurrent():
    valid=False
    while not valid:
        try:
            I=float(raw_input('Please enter the desired operating current in mA'))*1e-3
            valid=True
        except ValueError: print 'not valid, try again'
    return I

def askSampleConductivity():
    valid=False
    while not valid:
        try:
            I=float(raw_input('Please enter the sample conductance to be used in the temperature screen '))
            valid=True
        except ValueError: print 'not valid, try again'
    return I





def sumFinite(i):
    if type(i) is not np.ndarray:
        i=np.asarray(i)
    return np.sum(i[np.isfinite(i)])

def replaceNonFinite(i):
    if type(i) is not np.ndarray:
        i=np.asarray(i)
    out=np.zeros(len(i))
    for num,val in enumerate(i):
        if np.isfinite(val):
            out[num]=val
    return out


#creates a nested dict of metal layer parameters. Thickeness, width and resistivity as their own entires
def constructMetalDicts(padded_dict,current):
    dxprime=hidx/np.cos(np.deg2rad(45))
    Pd={}
    Au={}
    Metals={'Au':Au,'Pd':Pd}
    for key,value in padded_dict.iteritems():
        if key=='Pd':
            print 'Doing Pd'
            #function is fed effective widths, so we have to half these for elec calcs
            Pd['width']=transform_heater(value[:,1]/2.)    #correct width in direction of current flow for 45 degree angle
            Pd['resistivity']=0.238e-6
            #thickness = 40nm anywhere there is non-zero width
            Pd['thickness']= [40e-9 if x!=0. else 0. for x in iter(value[:,1])]

            #now reduce the thickness at any point that was deposited on the pyramid
            Pd['thickness']=[Pd['thickness'][count] if n<=139 else correctEvapOnPyramid(Pd['thickness'][count]) for count, n in enumerate(value[:,0])]
            print '----> Getting the fixed width of Pd in x\' direction'
            Pd['fixedwidth']= estimatePdWidth(Pd['width'],doTransform=False,effectiveWidth=False) #set to False as we already halved these widths
            #this is rho/A, like it says. Just splitting it up into seperate lines
            rhoOverA=(np.array(Pd['resistivity']))/(Pd['fixedwidth']*np.array(Pd['thickness']))
            #R=rhoL/A, but L varies depending on whether or not the sample is on the pyramid
            Pd['elec_resistance']=[rhoOverA[count]*lowdx if n<=139 else rhoOverA[count]*tipTilt(hidx) for count, n in enumerate(value[:,0])]
            #rho/A has some infs. We can sum all but infs using numpy, but we have to convert to an array from a list
            Pd['elec_resistance']=np.asarray(Pd['elec_resistance'])
            #perform the summation (excluding infs) and round to an int.

            PdResistance=sumFinite(Pd['elec_resistance'])*2.
            print '----> The total electrical resistance of the Pd sensor is {0:.2f}Ω'.format(PdResistance)
            fig_of_merit['Palladium Resistance']= ['{0:.2f}'.format(PdResistance)]
            #replace with user current
            Pd['power']=(current**2.)*Pd['elec_resistance']*(Pd['width']/Pd['fixedwidth'])*2.
            Pd_power_mW=sumFinite(Pd['power'])*1e3
            fig_of_merit['Palladium Power'] = ['{0:.2f}mW'.format(Pd_power_mW)]
            #numpy has a NaNsum <3. Wish I knew about this earlier
            print '----> Total power generated by Pd region; {0:.2f}mW'.format(Pd_power_mW)

        elif key =='Au':
            print 'Doing Au'
            Au['width']=value[:,1]/2.
            Au['resistivity']=47.847e-9
            Au['thickness']= [150e-9 if x!=0. else 0. for x in iter(value[:,0])]
            Au['thickness']=[Au['thickness'][count] if n<=139 else correctEvapOnPyramid(Au['thickness'][count]) for count, n in enumerate(value[:,0])]
            #calculate the x direction invariant parts first

            #THIS BIT IS COMPLETELY WRONG
            _k=(2e-3**2)*Au['resistivity']


            numerator = np.asarray([Au['resistivity']*lowdx if n<=139 else Au['resistivity']*tipTilt(dxprime) for n in iter(value[:,0])])
            denominator = np.asarray([Au['width'][n]*t for n,t in enumerate(Au['thickness'])])
            Au['elec_resistance'] = numerator/denominator
            fig_of_merit['Gold Resistance - Active']= ['{0:.2f}'.format(sumFinite(Au['elec_resistance'])*2.)]
            print '----> The total electrical resistance of Au active wiring is {0}Ω'.format(sumFinite(Au['elec_resistance'])*2.)
            Au['power']=current**2.*Au['elec_resistance']
            fig_of_merit['Gold power'] = ['{0:.2f}'.format(sumFinite(Au['power'])*1e3)]

            print '----> The total power generated by Au is {0}mW'.format(sumFinite(Au['power'])*1e3)
    totalpower=replaceNonFinite(Au['power'])+replaceNonFinite(Pd['power'])
    print '----> The total (lumped) power generation of all materials is {0:.4}mW'.format(np.sum(totalpower)*1e3)

    return totalpower,PdResistance

#takes the generated power from Au and Pd and sums them (needed for getTemp calc)
def powerSum(metal_dict):
    print '--> Summing distributed Joule heating contributions'
    totalpower=np.sum((metal_dict['Pd']['power'],metal_dict['Au']['power']),axis=0)
    print '----> Total (lumped) power generation: {0}mW'.format(np.sum(totalpower)*1e3)
    return totalpower

#no air in this version
def getTemp(parallel_resistance,summed_power,sample_resistance=np.inf,user_input=False,contact=False):
    #contact conductance, reciprocal of contact resistance, from 8.33e5 in YF's thesis
    if contact==False:
        gc = 0
    else:
        gc=1.2e-6

    if user_input==True: sample_resistance=askSampleConductivity()**-1
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
        A[n-1,n-1]=g[n-1]+(sample_resistance**-1)+gc #final diagonal position (n,n)


    Ainv=np.linalg.inv(A)
    V=np.dot(Ainv,summed_power)

    return V


#****************************
#value from Kamil's Thesis.
# Corrects the difference due to the writing plane being different to the projection plane
def tipTilt(num):
    return num/(np.cos(np.deg2rad(46.5)))

def correctEvapOnPyramid(thickness):
    return thickness*np.cos(np.deg2rad(46.5))


#******************************

def drawProbeNiceLayers(filtered, skip=False):
    if skip==True:
        print '--> Skipping probe schemtaic drawing'
        return
    else:
        print '--> Drawing schematic'

        #these two for loops are seperate to ensure that the SiN layer is drawn first!

        for key in filtered.keys():
            if len(filtered[key])==0: continue

            if 'SiN' in key:
                print '----> Drawing SiN layer first'
                for item in filtered[key]:

                    item.show(fdict['Schematic'])   #show item on 'img' axis
                    item.reflect('x').show(fdict['Schematic']) #reflect and show
                break

        for key in filtered.keys():
            if len(filtered[key])==0: continue

            if 'SiN' not in key:
                print '----> Drawing metal layers on top'
                for item in filtered[key]:
                    item.show(fdict['Schematic'])
                    item.reflect('x').show(fdict['Schematic'])

def collapseStackedElements(sorted_dict):
    print '--> Collapsing stacks of elements & calculating effective widths'
    output=newDictSameKeys(sorted_dict)

    #for every material in the dictionary
    for key,value in sorted_dict.iteritems():

        x_pos=[]
        #for each segment, get the lower left coordinate (X0)
        for n in value:
            x_pos.append(getX0(n))
        #report ther number of unique entries
        x_pos=np.unique(x_pos)
        print '----> Processing {0}: {1} unique elements found'.format(key,len(x_pos))


        #this is required in the case that we have two elements of the same material
        #occupying the same x position but different y positions. We make a list
        #of all the distinct positions a material occupies, then sum the areas of
        #all the boundaries that occur at that position.
        list2d=[[] for _ in range(len(x_pos))]
        #iterate over the x positions array once inserting the Boundaries present at the locations in the relevant sample of the 2dlist
        for counter,number in enumerate(x_pos):
            for element in value:
                if getX0(element)==number:
                    list2d[counter].append(element)

        #iterate over it again performing the calculation
        for n, v in enumerate(x_pos):

            placeholder=[]
            for entry in list2d[n]:
                #if the entry is a boundary object, rather than a float position indicator
                if type(entry)!=np.float32:
                    placeholder.append(entry.area())
            if v<=139.:
                output[key].append([v,2*(((sum(placeholder))*1e-12)/lowdx)])


            else:
                output[key].append([v,2*(((sum(placeholder))*1e-12)/hidx)])


    #turn all the lists into numpy arrays
    for key,value in output.iteritems():
        output[key]=np.array(output[key])


    return output


#for this we need the contact radius/ radius of heating, let us assume this is equal to the width of the last sample in SiN
def conductivityToResistance(k,contact_area):
    return 1./(4.*k*contact_area)


def getAverageTemp(Tdist, position):
    #find the index of the position supplied
    marker=np.where(collapsed['SiN'][:,0]>position)[0][0]
    print np.average(Tdist[marker:])
    avgT=round_to_n(np.average(Tdist[marker:]),4)
    print '-->Average temperature increase over tip region = {0:.2f}K'.format(avgT)
    return avgT

def plotSweep(probeR,probeP,sensor_pos,contact_area,ax):
    #set up ten log spaced samples of sample thermal conductivity (spans most materials)
    sample_conductivities=np.logspace(-3,3,8)
    #convert them to thermal resistances for adding to the probe network
    sample_Rth=conductivityToResistance(sample_conductivities,contact_area)
    #empty list for resulting probe temperatures
    probe_response=[]

    for sample in sample_Rth:
        #add the average tip temperature to the list
        probe_response.append(getAverageTemp(getTemp(probeR,probeP,sample,contact=True),sensor_pos))

    #plot the points. Take a note of the line so as its color can be grabbed for later
    pointsLine=ax.plot(sample_conductivities, probe_response,'o',label='Temperature Simulations')
    linecolor=pointsLine[0].get_color()
    #get the coefficients of the fit
    c1=fs2.getCoeffs(sample_conductivities,probe_response,doPrint=True)

    #'define an x array for the plotting of the curve (gives higher res than otherwise)'
    fit_x=np.logspace(-3,3,100)

    #'apply the curve function to the x array'
    y2 = fs2.CMIcurve(fit_x, *c1)

    #'plot the fitted curve'
    ax.plot(fit_x,y2,color=linecolor,label='Sigmoid Fit')

    #'calculate the derivative and plot'
    twinax=ax.twinx()
    twinax.plot(fit_x,(fs2.CMIderivative(fit_x,c1[0],c1[1])),ls='dashed',color=linecolor, label='Derivative of Fit')

    ################

    #trick the legend element into adding the entry for the derivative in the ax1 list
    #plots nothing, but defines a red line with the derivative label
    ax.plot(0,0,'--',color=linecolor,label='Derivative of Fit')
    nticks=7
    ax.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(nticks))
    twinax.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(nticks))

    #set the ticks to be the same for both axes by defining 5 tick marks for each curve

    #ax.set_yticks(np.linspace(ax.get_ybound()[0], ax.get_ybound()[1], 5))
    #ax.set_yticks(np.linspace(0, 250, 5))
    #twinax.set_yticks(np.linspace(-50, 0, 5))
    peakLoc=fit_x[np.argmin((fs2.CMIderivative(fit_x,c1[0],c1[1])))]
    peakVal=np.min((fs2.CMIderivative(fit_x,c1[0],c1[1])))
    fig_of_merit['Maximum Sensitivity Magnitude (K/W/mK)'] = ['{0:.2f}'.format(peakVal)]
    fig_of_merit['Maximum Sensitivity Location (W/mK)'] = ['{0:.2f}'.format(peakLoc)]
    ax.axvline(peakLoc,ls='dashed',color='darkgrey',label='Maximum Sensitivity of {0:.2f}K per W/mK at {1:.2f}W/mK'.format(peakVal,peakLoc))
    ax.set_ylim(100,160)
    #ax.text(10., 150, 'Maximum Sensitivity
    twinax.set_ylim(-5,1)
    ax.set_ylabel('Temperature Change, dT (K)')
    ax.set_xlabel('Thermal Conductivity, log(k) (W/mK)')
    twinax.set_ylabel('Sensitivity (${\delta T}/ \delta k $)')
    ax.legend(loc=3)
    ################

    ax.set_xscale('log')

def annotateTemp(plot,sensor_pos,avgT):
    plot.axhline(avgT,color=colour_dict['Pd'],ls='--')
    plot.axvline(sensor_pos,color=colour_dict['Pd'],ls='--')
    plt.text(0, avgT,'Average sensor temperure rise = {0:.2f}K'.format(avgT))

#using the experimentally determined value for maximum sustained current density
#calculate the driving current in mA recommended for probe operation.
def checkCurrentDensity(Pd_widths):
    J=8.56e10  #from Yunfeis experiments
    #get the x-section, assuming the minimum Pd thickness (found on pyramid)
    print '--> Calculating maximum allowable current given Pd dimensions'
    a=estimatePdWidth(Pd_widths)*correctEvapOnPyramid(40e-9)
    maxI=J*a    #current is just current density time x-section
    print 'Maximum operating current: {0:.2f}mA'.format(maxI*1e3)

def estimatePdWidth(Pd_widths,doTransform=True,effectiveWidth=True):
    #implement an option for dividing the width in here
    scaling = 2. if effectiveWidth else 1.
    #get the unique numbers, and how many of each, discounting zeros
    vals,counts=np.unique(Pd_widths[Pd_widths.nonzero()],return_counts=True)
    #find the indice of the most freqently occuring number (like mode of the width)
    ind = np.argmax(counts)
    #by default, do the 45degree transform on this value/2. (assumes the input is an effective width)
    if doTransform==True:
        print transform_heater(vals[ind]/scaling)
        return transform_heater(vals[ind]/scaling)
    #optionally skip the transformation (for when the function is applied to pre-transformed data)
    else:
        print vals[ind]/scaling
        return vals[ind]/scaling

def checkPdElectricalResistance(Pd_widths):
    dxprime=dx/np.cos(np.deg2rad(45))
    rhoPd=238e-9
    #get length along xprime axis
    l=tipTilt(2*dxprime*len(Pd_widths.nonzero()[0])) #element [1] of nonzero is the dtype?

    R=(rhoPd*l)/(estimatePdWidth(Pd_widths)*correctEvapOnPyramid(40e-9))
    print 'Total electrical resistance of Pd region: {0:.2f}'.format(R)

#################cantilever gold resistance#####################################
def getElementalElectricalResistance(metal_dict):
    elecR=newDictSameKeys(metal_dict)
    for key,metal in metal_dict.iteritems():
        if 'Pd' in key:
            L=tipTilt(dx/np.cos(np.deg2rad(45)))
        else: L=dx

        elecR[key]=np.zeros(len(metal['width']))
        A=(metal['width']/2.)*metal['thickness']
        elecR[key]=np.array(metal['resistivity']*L/(A)) #get most frequent
        elecR[key]=np.sum(elecR[key][elecR[key]<np.inf])*2.
    return elecR
################################################################################

def reportRth(dictionary):
    print '----> Reporting thermals resistances of all materials'
    for key,val in dictionary.iteritems():
        print '-------->{0}: {1:.2f}x10^6K/W'.format(key,np.sum(val)/1e6)

#c is the dictionary of collapsed widths
def metalPlacement(c):
    #prepare an output dictionary with the same keys
    output=newDictSameKeys(c)

    output['SiN']=c['SiN']
    #whilst iterating, make sure to ignore SiN, because its fine.
    for key,value in c.iteritems():
        if key=='SiN':
            print 'hit {0}, ignoring'.format(key)
            continue
        #if there is nothing in the array, then ignore it
        if len(value)==0:
            print 'nothing in {0}, ignoring'.format(key)
            continue
        print 'padding {0}'.format(key)

        #if none of the above, then its a metal layer to be padded.
        xax=np.copy(c['SiN']) #copy the x axis from the SiN layer, as this is continuous
        xax[:,1]=0          #set widths to be zero

        #scan through the x axis by index and value simultaneously
        #if at any time the value you find is equal to the value of
        #the first entry of a metal layer, c[key][:,0][0], note down the array index where this
        #occurred. this is all because of irregular sample spacing on the x axis
        #we are trying to relate the position in the arrays to the true x position
        for n,v in enumerate(xax[:,0]):
            if v==c[key][:,0][0]:
                pointer=n

        #c[key][:,0][0]
        #in the dictionary c, under entry 'key', on the first column (real x values), report the first (0th) number

        #the widths axis of xax, which is a blank copy, is to be filled with
        #the width data from the metal 'key', in the position beginning at 'pointer'
        #and etending to the length of the metal array.
        xax[:,1][pointer:pointer+c[key].shape[0]]=c[key][:,1]

        output[key]=np.copy(xax)

    return output

#%%<AssySim>
############################# ASSY Experiment###################################
def AssySim():
    Isweep=np.array((0.25,0.50,0.75,1.00,1.25,1.50,1.75,2.00))*1e-3
    Psweep=np.array(())
    Tsweep=[]

    for I in Isweep:
        sum_power,R=constructMetalDicts(padded,I)
        Psweep=np.append(Psweep,(I**2)*R)
        T=getTemp(thermal_resistances['Parallel'],sum_power,user_input=False)
        Tsweep.append(getAverageTemp(T,sensor_pos))
        f3,ax3=plt.subplots(1)
    coefs= poly.polyfit(Psweep,Tsweep,1)
    ax3.plot(Psweep*1e3,Tsweep)
    ax3.set_xlabel('Power (mW)')
    ax3.text(0,140, 'Thermal Resistance={0:.2E} K/W'.format(coefs[1]))
    ax3.set_ylabel('Average Sensor Temperature Rise (K)')
    f3.show()

    print 'Assy style thermal resistance = {0:.2E} K/W'.format(coefs[1])
    return coefs[1]
#%%<Main:file prep>
################################# MAIN #########################################


print 'Welcome to Cantilever Modeller GDSII Edition!'
print '--> Beginning data import and preparation'


#load the file


#filename='.gds generated for testing/out5.gds'
filename=choosefile().encode('ascii', 'ignore')
elements=get_elements_from_GDS(filename)
#filter by layer type
filtered=filter_elements_by_layer(elements)


check_dimensions_of_elements(filtered)
#sort the layers by x position
filteredAndSorted=sort_in_x(filtered)

#get the dx's
#lowdx,hidx=get_dual_dx(filteredAndSorted['SiN'])
####MIGRATED TO HARDCODED DX's in cmodel_funcs

#transforms all geometry in a given dx slice to an effective width
#this includes slices which have two distinct elements
collapsed=collapseStackedElements(filteredAndSorted)

#put all the metal elements in the right place, surrounded by zeros
padded=metalPlacement(collapsed)


#get the x axis as we will need this for plotting and external heating
xax=padded['SiN'][:,0]


"""
#########Plot effective widths##########
for val in padded.itervalues():
    if len(val)!=0:
        plt.plot(val[:,0],val[:,1])
plt.show()
#######################################
"""

#prepare the metal layers and calculate the total (elemental) power gen
sum_power,R=constructMetalDicts(padded,2e-3)


#external_heat=add_power.generate_external_heat(1e-3, 30, 40, xax)

#if len(sum_power)==len(external_heat):
    #sum_power +=external_heat

"""
transform_heater(1.5)
J=2.5e-3/(transform_heater(1.5e-6)*correctEvapOnPyramid(40e-9))
print J/1e10
checkCurrentDensity(padded['Pd'])
"""

print '--> Data processing complete. Beginning simulation'



#%%<Main:simulation>
############################### SIMULATION #####################################

#calculate elemental resistances
thermal_resistances=getElementalResistance(padded)
#before inputting a current value, sanity check the maximum allowable current
#checkCurrentDensity(collapsed['Pd'])
#get the Joule heating power. Just append it to the metals' dictionary under key 'power'
#getJouleHeat(Metal)
#work out the contact radius before doing any temperature distributions
#this is important for calculation of the sample spreading resistance Rs=1/4ka
contact_area=2.*hidx
print '--> Calculated contact area as {0}nm'.format(contact_area*1E9)
print '--> Solving elemental temperature matrix'
#get the temperature rise given the selected current (hardcoded)
T=getTemp(thermal_resistances['Parallel'],sum_power,user_input=False)

print '--> Temperature distribution solved, drawing plots'

#%%<Main:plotting>
############################# PLOT 1 ##########################################

sensor_pos=145. #sensor position in x, in um
avgT=getAverageTemp(T,sensor_pos)
fig_of_merit['Out of contact sensor temperature rise (K)']=['{0:.2f}'.format(avgT)]

drawProbeNiceLayers(filteredAndSorted, skip=False)
plotPower(sum_power,p,xax,True)
plotDict(thermal_resistances,r,xax,True)
plotArray(T,t, xax)
annotateTemp(fdict[t],sensor_pos,avgT)
f.show()

print '--> Attempting to show first plot'

############################## PLOT 2 #########################################

f2,ax2=plt.subplots(1)
plotSweep(thermal_resistances['Parallel'],sum_power,sensor_pos,contact_area,ax2)
f2.show()



print 'Simulation Complete'


Rth=AssySim()
fig_of_merit['Assy Style Thermal Resistance (K/W)'] = ['{0:.2E}'.format(Rth)]


############################Table of useful values#############################

theDir='/'.join(filename.split('/')[:-1])+'/'
theFile=''.join(filename.split('/')[-1][:-4])


fig_of_merit['GDS Name']=[theFile]
table= tabulate.tabulate(fig_of_merit, headers='keys', tablefmt='html',stralign='center',numalign='center')
f = open(theDir+theFile+'_Results.txt', 'w')
f.write(table)
f.close()

theFile
############################Save the results####################################

simresult={'k':thermal_resistances['Parallel'], 'p':sum_power, 'x':xax, 't':T}

np.save(theDir+theFile+'_Data.npy', simresult)
