%matplotlib

#%% <cell>
import os.path
import datetime

#%%
import gdspy
import numpy as np
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly
from math import log10, floor
import matplotlib
from matplotlib import style
from matplotlib import cm
style.use('ggplot')
from matplotlib import rc
import fitsig2 as fs2
import matplotlib.ticker
rc('text', usetex=True)
rc('figure', figsize=(22,14.2))
matplotlib.rcParams['font.serif'] = 'CMU Serif'
matplotlib.rcParams['font.family'] = 'serif'

matplotlib.rcParams['font.size'] = 22
LARGE_FONT= ("Segoe UI Semilight", 12)
BUTTON_FONT = ("Segoe UI", 10)

print 'hello'

gdsii = gdspy.GdsLibrary()
#read in the gds file to an object named layout
layout=gdsii.read_gds('./gds generated for testing/Dual dx tested/Original and quick cut/OriginalCant_unfractured_hole.gds')
#print the cell names in the gds file
print layout.cell_dict
#extract the cell
c=gdsii.extract('Half_cant_copy')

SiN_poly=[]
Au_poly=[]
Pd_poly=[]
for n in c.elements:
    if n.layer==0:
        SiN_poly.append(n)
    if n.layer==1:
        Au_poly.append(n)
    if n.layer==2:
        Pd_poly.append(n)


cut=gdspy.Cell('GAH')
for n in range(0,139):
    SiN_cut=gdspy.slice(SiN_poly, [n,n+1], 0, layer=0)
    Au_cut=gdspy.slice(Au_poly, [n,n+1], 0, layer=1)
    Pd_cut=gdspy.slice(Pd_poly, [n,n+1], 0, layer=2)
    cut.add(SiN_cut[1])
    cut.add(Au_cut[1])
    cut.add(Pd_cut[1])

tipArray=np.linspace(139., 150., 441)

for n in tipArray:
    SiN_cut=gdspy.slice(SiN_poly, [n,n+0.025], 0, layer=0)
    Au_cut=gdspy.slice(Au_poly, [n,n+0.025], 0, layer=1)
    Pd_cut=gdspy.slice(Pd_poly, [n,n+0.025], 0, layer=2)
    cut.add(SiN_cut[1])
    cut.add(Au_cut[1])
    cut.add(Pd_cut[1])

anotherCell=gdspy.GdsLibrary()
anotherCell.add(cut)
anotherCell.write_gds('out4.gds')
