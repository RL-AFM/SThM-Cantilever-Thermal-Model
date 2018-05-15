%matplotlib
import numpy as np
import os
from myStyle import *
from scipy.interpolate import spline
import numpy.polynomial.polynomial as poly
matplotlib.rcParams['font.size'] = 22
def plotPower(arr,plotref,xax,ylabel,color='red',log=True, scaling=True):
    ax=plotref
    if scaling==True: scaled=[arr[n] if xax[n]<=139 else arr[n]*(lowdx/hidx) for n,v in enumerate(xax)]
    else: scaled=arr
    ax.plot(xax,scaled,color=color)
    ax.set_xlabel('Length ($\\mu m$)')
    ax.set_ylabel(ylabel)
    #ax.text(0,1e-5, 'Total power generated = {0:.2f}mW'.format(np.sum(arr)*1e3))
    if log==True:ax.set_yscale('log')

lowdx=1e-6
hidx=25e-9


theDir='\\'.join(os.getcwd().split('\\')[:-1])+'\\gds generated for testing\\180312 Post JWPD\\'

cutout=[]

cutout.append(np.load(theDir+'Cut0_Fractured_Data.npy').item())
cutout.append(np.load(theDir+'Cut4_Fractured_Data.npy').item())
cutout.append(np.load(theDir+'Cut6_Fractured_Data.npy').item())
cutout.append(np.load(theDir+'Cut8_Fractured_Data.npy').item())
cutout.append(np.load(theDir+'Cut10_Fractured_Data.npy').item())


f,ax= plt.subplots()

n=len(cutout)
color=iter(cm.afmhot(np.linspace(0.25,0.75,n)))

for val in cutout:
    c=next(color)
    plotPower(val['t'], ax, val['x'],'Temperature Rise (K)',log=False,scaling='False',color=c)


f.show()


f2,ax=plt.subplots()
assyRes=[355000,355000,360000,384000,441000,443000,448000,455000,462000]
avgT=[152.9,152.1,155.1,165.4,190.0,190.8,193.2,196.3,198.9]
cutsize=[0,0.5,1,2,3,4,6,8,10]
ax.plot(cutsize,avgT,'o',label='Simulated Points')

ax.plot(cutsize,assyRes,'o',label='Simulated Points')
ax.axvline(4,label='End of sensor')
ax.set_xlabel('Cutout Length ($\mu m$)')
ax.set_ylabel('Effective Thermal Resistance ($K/W$)')
plt.legend(loc='best')


#polyfits for assysim. Didnt really work
"""
coefs=poly.polyfit(cutsize[:4], assyRes[:4], 2)
ffit = poly.polyval(np.linspace(0,3,20), coefs)
plt.plot(np.linspace(0,3,20), ffit)
coefs2=poly.polyfit(cutsize[5:],assyRes[5:],1)
ffit2=poly.polyval(np.linspace(3,10,20), coefs2)
plt.plot(np.linspace(3,10,20), ffit2,c='lightblue',label='Fit')
"""
plt.show()
