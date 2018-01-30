import numpy as np
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly

R=278

I=np.array((0.5,1.0,1.5,2.0))*1e-3
P=(I**2)*R

T=np.array((21.47,85.87,193.20,343.50))+23

plt.plot(P,T)
plt.show()
coefs= poly.polyfit(P,T,1)
print '{0:.2E} K/W'.format(coefs[1])
