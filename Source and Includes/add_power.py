import numpy as np

#given a lumped power and an extent, spread it out across the defined range.
#this function takes into account the coarse and fine sampling change
def generate_external_heat(lumped, start, end, xarr):
    print 'Generating an extrnal heat source  of total power {0}mW'.format(lumped*1e3)
    scaledpower=lumped/(end-start)      #get the power per micron
    print 'Given the extents provided, this corresponds to a distributed power of {0}uW/um'.format(scaledpower*1e6)
    dist_power=[]                       #output array for the distributed power
    for n,v in enumerate(xarr):
        if not start<v<end:
            dist_power.append(0)
            continue
        else:
            dist_power.append(scaledpower)
    return dist_power
