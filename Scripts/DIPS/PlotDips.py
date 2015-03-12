import scipy.interpolate as ip
import numpy as np
from matplotlib import pyplot
import jdcal

def StellaEnviro(fname):

    humdl = np.loadtxt(fname,dtype='string',skiprows=0,delimiter=',')
    hums = np.array(humdl[:,1],dtype='float')
    
    datetime = np.array([p.split(' ') for p in humdl[:,0]])
    
    date = np.array( [day.split('.') for day in datetime[:,0]] ,dtype='float')
    time = np.array( [t.split(':')   for t   in datetime[:,1]] ,dtype='float')
    
    jd  = np.array([jdcal.gcal2jd(d[2],d[1],d[0]) for d in date])

    mjd = jd[:,1] + time[:,0]/24. + time[:,1]/(60.*24.)

    huminterp = ip.interp1d(mjd,hums)
    return huminterp


jd0 = 56244.
data = np.loadtxt('DipAbsorption.txt')

huminterp = StellaEnviro('stella-environment.txt')

gd = np.where(data[-1,:]+jd0 < 57020)[0]

#pyplot.plot(data[-1,gd],huminterp(data[-1,gd]+jd0)/100.,'o')

print np.median(data[7,gd]),np.median(data[6,gd]),np.median(data[3,gd]),np.median(data[2,gd])
print np.median(data[7,gd])/np.sin(30.*np.pi/180.),np.median(data[6,gd])/np.sin(30.*np.pi/180.),np.median(data[3,gd])/np.sin(30.*np.pi/180.),np.median(data[2,gd])/np.sin(30.*np.pi/180.)


colors = np.random.uniform(size=gd.size,low=0.4)
hums = huminterp(data[-1,gd]+jd0)
goodday = np.where(hums < 3.)[0]
badday  = np.where(hums > 70.)[0]
freqs = [11.,13.,17.,19.]
vals  = [7+8,6+8,3+8,2+8]
cvals  = [5+8,4+8,1+8,0+8]

print gd[goodday]
test = data[:,gd[goodday]]
test = test[vals,:]
ctest = data[:,gd[goodday]]
ctest = ctest[cvals,:]

print test.shape,np.median(test,axis=1)-2.73
#data[vals,gd[badday[0]]]
pyplot.scatter(freqs,np.median(test,axis=1)-2.73,c='#F22000',s=110,lw=1,label='Uncorrelated Channels')
pyplot.scatter(freqs,np.median(ctest,axis=1)-2.73,c='#59F200',s=80,lw=1,label='Correlated Channels',alpha=0.8)
pyplot.xlabel('Frequency (GHz)',size=16)
pyplot.ylabel('System Temperature (K)',size=16)
pyplot.ylim(0.5,54.)
pyplot.xticks(size=14)
pyplot.yticks(size=14)
pyplot.legend(frameon=False,numpoints=1)
pyplot.show()



hums = huminterp(data[-1,gd]+jd0)
goodday = np.where(hums < 3.)[0]
badday  = np.where(hums > 70.)[0]
freqs = [11.,13.,17.,19.]
vals  = [7,6,3,2]

print gd[goodday]
test = data[:,gd[goodday]]
test = test[vals,:]
btest = data[:,gd[badday]]
btest = btest[vals,:]

print test.shape,np.median(btest,axis=1)
#data[vals,gd[badday[0]]]
pyplot.scatter(freqs,np.median(btest,axis=1),c='#F22000',s=80,lw=1,label='High Humidity')
pyplot.scatter(freqs,np.median(test,axis=1),c='#59F200',s=80,lw=1,label='Low Humidity')
pyplot.xlabel('Frequency (GHz)',size=16)
pyplot.ylabel('Percentage Absorption',size=16)
pyplot.ylim(0.5,4.)
pyplot.xticks(size=14)
pyplot.yticks(size=14)
pyplot.legend(frameon=False,numpoints=1)
pyplot.show()


#Plot Time Vs Absorption
colors = np.random.uniform(size=gd.size,low=0.4)
pyplot.scatter(data[-1,gd],data[2,gd],c=colors,cmap='gray',s=80,lw=0)
pyplot.xlabel('Observation Day',size=16)
pyplot.ylabel('Percentage Absorption',size=16)
pyplot.ylim(0.5,4.)
pyplot.xticks(size=14)
pyplot.yticks(size=14)
pyplot.show()


#Plot Humidity Vs Absorption
colors = np.random.uniform(size=gd.size,low=0.4)
pyplot.scatter(huminterp(data[-1,gd]+jd0),data[6,gd],c=colors,cmap='gray',s=80,lw=0)
pyplot.xlabel('Relative Humidity (%)',size=16)
pyplot.ylabel('Percentage Absorption',size=16)
pyplot.ylim(0.5,4.)
pyplot.xlim(0,100)
pyplot.xticks(size=14)
pyplot.yticks(size=14)
pyplot.show()
