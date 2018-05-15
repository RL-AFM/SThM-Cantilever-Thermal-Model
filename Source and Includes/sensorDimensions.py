import numpy as np
import tabulate
import os
from IPython.display import HTML, display
t=21e-9
rho=238e-9
Jmax=1.1e11
res=110

def R(L,w):
    A=w*1e-6*t
    R=(rho*L*1e-6)/(A)
    print 'Resistance calculated as {0}Ohms'.format(R)
    return R

def P(i,R):
    print 'Power = {0:.2f}mW'.format((i*1e-3)**2*R*1e3)


PurpPd=R(4,1.175)
OuterPd=R(5.8, 0.290)
SensorPd=R(2.24,0.290)

P(0.67, PurpPd)
P(0.67, OuterPd)
P(0.67, SensorPd)

def getWidth(L):
    width=[]
    i=0
    for length in L:
        Lcor=length/np.cos(np.deg2rad(45))*2
        width.append((rho*Lcor*1e3)/(t*res))
        print 'For Length {0}um in x direction ({2:.3f}um in 45degrees), Resistor width = {1:.0f}nm'.format(length,width[i],Lcor/2)
        i+=1
    return np.array(width)


def maxcurrent(w):
    A=w*1e-9*t
    Imax=Jmax*A*1e3
    #print 'Maximum current calculated as {0}mA'.format(Imax)
    return Imax

Ls=np.arange(0.5,6,0.5)
LsCorrected=(Ls/np.cos(np.deg2rad(45)))*2
w=getWidth(Ls)
I=maxcurrent(w)



data=np.array((Ls,np.round(w),np.round(I,2)))
data=data.transpose()
tableobj=tabulate.tabulate(data, headers=('Length in x (um)','Width (nm)', 'Maximum Current (mA)'), tablefmt='html')
display(HTML(tableobj))
theDir='\\'.join(os.getcwd().split('\\')[:-1])+'\\gds generated for testing\\180312 Post JWPD\\'
f = open(theDir+'sensorDimsTable'+'_Results.html', 'w')
f.write(tableobj)
f.close()
