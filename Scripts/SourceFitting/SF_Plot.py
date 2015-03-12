                            fig = pyplot.figure()
                            samps = np.arange(volts.size)
                            isGaussian   = (x**2 + y**2 < beam['X'][0]*0.8**2) 
                            isBackground = (x**2 + y**2 < beam['X'][0]*6.**2) & (x**2 + y**2 > beam['X'][0]*3.**2)
                            mid = (samps > 81850) & (samps < 81950)
                            low = (samps < 81850)
                            high = (samps > 81950)

                            pyplot.plot(volts*1000.,'-',color='#ADADAD',label = 'Measured Voltage')
                            pyplot.plot(SourceFitting.lm_SourceFit.ModelFit(P1,x,y)*1000.,'-',lw=3,alpha=0.5,color='#0059F2',label='Model Fit')
                            pyplot.plot(samps[((isBackground | isGaussian) & mid)],
                                        SourceFitting.lm_SourceFit.ModelFit(P1,x[((isBackground | isGaussian) & mid)],y[((isBackground | isGaussian) & mid)])*1000.,
                                        '--',lw=3,alpha=0.8,color='#F29900',label='Model Fit Region')
                            pyplot.legend(frameon=False)
     
                            pyplot.plot(samps[((isBackground | isGaussian) & low)],
                                        SourceFitting.lm_SourceFit.ModelFit(P1,x[((isBackground | isGaussian) & low)],y[((isBackground | isGaussian) & low)])*1000.,
                                        '--',lw=3,alpha=0.8,color='#F29900')
                            pyplot.plot(samps[((isBackground | isGaussian) & high)],
                                        SourceFitting.lm_SourceFit.ModelFit(P1,x[((isBackground | isGaussian) & high)],y[((isBackground | isGaussian) & high)])*1000.,
                                        '--',lw=3,alpha=0.8,color='#F29900')

                            
                            pyplot.xticks(size=14)
                            pyplot.yticks(size=14)
                            pyplot.xlabel('Sample',size=16)
                            pyplot.ylabel('Voltage (mV)',size=16)
                            pyplot.show()
