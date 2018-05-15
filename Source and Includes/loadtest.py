%matplotlib

"""Test Module for loading numpy data dumps. The idea being that the model can
save the result arrays to an external file. This can then be loaded by a custom
script for plotting comparisons between datasets"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import style
from matplotlib import cm
from matplotlib import rc
style.use('ggplot')
rc('text', usetex=True)
rc('figure', figsize=(22,14.2))
matplotlib.rcParams['font.serif'] = 'CMU Serif'
matplotlib.rcParams['font.family'] = 'serif'
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

foo=np.load('test.npy').item()
bar=np.load('testcold.npy').item()
foo.keys()
f,ax= plt.subplots()
foo.keys()
bar.keys()
plotPower(foo['t'], ax, foo['x'],'Temperature Rise (K)',log=False,scaling='False')
plotPower(bar['t'],ax,bar['x'],'power2','blue',log=False,scaling='False')

f.show()
