import numpy as np
from matplotlib import pyplot


#b = np.concatenate((np.ones(50),np.cos(np.radians(np.linspace(0,180,25)))*0.5 +0.5, np.zeros(25) ))
#r = b[::-1]
#g = g = 1./np.cosh(np.radians(np.linspace(-180,180,100)) )

g = np.concatenate((np.cos(np.radians(np.linspace(28,180,70)))*0.5 + 0.5, np.zeros(25) ))
#g = b[::1]
#g = np.concatenate((np.cos(np.radians(np.linspace(10,180,95)))*0.375 +0.625 ))
#r = np.concatenate((1./np.cosh(np.radians(np.linspace(10,90,95)) )))
b = np.concatenate((np.cos(np.radians(np.linspace(40,180,30)))*0.5 + 0.5, np.zeros(65) ))#0.375 +0.725 -0.03
r = 1./np.cosh(np.radians(np.linspace(20,90,95)) )


np.savetxt('mfi-colortable.txt',np.transpose(np.array([r,g,b])))
pyplot.plot(r,'r')
pyplot.plot(b,'b')
pyplot.plot(g,'g')
pyplot.show()
