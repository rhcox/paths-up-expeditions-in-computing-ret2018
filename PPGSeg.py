import time
import math
from scipy.signal import find_peaks_cwt
import numpy as np
#from pypeaks import Data, Intervals
import csv
#import plotly
#import plotly.plotly as py
#import plotly.graph_objs as go
import matplotlib.pyplot as plt
from numpy import sin, linspace, pi
from scipy import fft, arange
import PPGWave
import PPGfunc

    
class PPGSeg(object):
    def make(self,PPGData):
        self.ppg=PPGData
#       self.smoothS=smoothSpan
#       self.smooth= [[0 for x in range(columns)] for y in range(rows)]
#       self.smoothN=[

        self.InflectionPoints=[]
        self.PulseWaveSystolicPeakPos=-1
        self.PulseWaveEndPos=-1
        self.DicroticNotchPos=-1
        self.DiastolicPeakPos=-1
        self.HalfHeight=-1
        self.HalfTimeLow=-1
        self.HalfTimeHigh=-1                        
        self.PulseWaveAmplitude=-1
        self.SystolicPhase=-1
        self.DiastolicPhase=-1
        self.PulseWaveDuration=-1
        self.PulseWavePropogationTme=-1
        self.PulseWidthTime=-1
        self.RiseTime=-1

        self.wf1=self.ppg.getWFnormP()
        #print "Waveform0", self.wf1[1]
        #data_obj = Data(self.wf1[0], self.wf1[1], smoothness=11)
        #print "Waveform1", self.wf1[1]
        #data_obj.get_peaks(method='slope')
        #print "Waveform2", self.wf1[1]
        self.wfD=self.ppg.getWFderivative()
        self.wfCP=self.ppg.getWFCriticalP()
        self.AverageWaveForm=[]
        #print "Waveform3", self.wf1[1]
        

        #print data_obj.peaks
        #print data_obj.peaks['peaks'][0]
        #print data_obj.peaks['peaks'][1]
        #print data_obj.peaks['valleys']

        #print 'Critical Points',self.wfCP


#        Pulse = go.Scatter(x=self.wf1[0],y=self.wf1[1],mode = 'markers',marker = dict(size = 1,color = 'rgba(255, 0, 0, .8)'),name='Pulse',)
#        Trace4 = go.Scatter(x=self.wfCP[0],y=self.wfCP[1],mode = 'markers',marker = dict(size = 3,color = 'rgba(255, 0, 0, .8)'),name='Pulse',)
        #Valleys = go.Scatter(x=data_obj.peaks['peaks'][0],y=data_obj.peaks['peaks'][1],mode = 'markers',marker = dict(size = 5,color = 'rgba(0, 255, 100, .8)'),name='Valleys',)
        #Peaks = go.Scatter(x=data_obj.peaks['valleys'][0],y=data_obj.peaks['valleys'][1],mode = 'markers',marker = dict(size = 5,color = 'rgba(0, 0, 255, .8)'),name='Peaks',)
#        data = [Pulse,Trace4]
#        layout = go.Layout(title='PPG',hovermode='closest',
#        xaxis=dict(title='Time [ms]',),
#        yaxis=dict(title='Pulse Waveform [Arb. Units]',),)
#        fig = go.Figure(data=data, layout=layout)
        #print data
#        py.plot(fig, filename='basic-scatter')
        #print "Out"



        self.AverageWaveForm=self.ppg.AverageWF()
        self.InflectionPoints=PPGfunc.inflection(self.AverageWaveForm[0],self.AverageWaveForm[1])
        self.PulseWaveSystolicPeakPos=int(self.AverageWaveForm[2][0])
        self.PulseWaveEndPos=int(self.AverageWaveForm[2][1])
        #print "Inflection Points= ", self.InflectionPoints, 'Length= ',len(self.InflectionPoints)
        if(len(self.InflectionPoints)>0):
            temp=self.InflectionPoints.pop()
            self.DicroticNotchPos=temp[0]
            self.DiastolicPeakPos=temp[1]
        else:
            self.DicroticNotchPos=-1
            self.DiastolicPeakPos=-1
        self.HalfHeight=0.5*self.AverageWaveForm[1][self.PulseWaveSystolicPeakPos]
        print 'Half Height= ',self.HalfHeight
        for i in range(0,self.PulseWaveSystolicPeakPos-1):
            if self.AverageWaveForm[1][i]>self.HalfHeight:
#                print 'Point i=(',self.AverageWaveForm[0][i],self.AverageWaveForm[1][i],')'
#                print 'Point i-1=(',self.AverageWaveForm[0][i-1],self.AverageWaveForm[1][i-1],')'
                DeltaX=self.AverageWaveForm[0][i]-self.AverageWaveForm[0][i-1]
#                print "Delta X=", DeltaX
                DeltaY=self.AverageWaveForm[1][i]-self.AverageWaveForm[1][i-1]
#                print "Delta y=", DeltaY
                DeltaV=self.HalfHeight-self.AverageWaveForm[1][i-1]
#                print "Delta V=", DeltaV
                self.HalfTimeLow=self.AverageWaveForm[0][i-1]+DeltaX*DeltaV/DeltaY
#                print "Half Time Low=", self.HalfTimeLow
                break
        for i in range(self.PulseWaveSystolicPeakPos+1,self.PulseWaveEndPos-1):
            if self.AverageWaveForm[1][i]<self.HalfHeight:
#               print 'Point i=(',self.AverageWaveForm[0][i],self.AverageWaveForm[1][i],')'
#                print 'Point i-1=(',self.AverageWaveForm[0][i-1],self.AverageWaveForm[1][i-1],')'
                DeltaX=self.AverageWaveForm[0][i]-self.AverageWaveForm[0][i-1]
#                print "Delta X=", DeltaX
                DeltaY=-self.AverageWaveForm[1][i]+self.AverageWaveForm[1][i-1]
#                print "Delta y=", DeltaY
                DeltaV=self.HalfHeight-self.AverageWaveForm[1][i]
#                print "Delta V=", DeltaV
                self.HalfTimeHigh=self.AverageWaveForm[0][i-1]+DeltaX*DeltaV/DeltaY
#                print "Half Time Low=", self.HalfTimeLow
                break
                        
        self.PulseWaveAmplitude=self.AverageWaveForm[1][self.PulseWaveSystolicPeakPos]
        self.SystolicPhase=self.RiseTime=self.AverageWaveForm[0][self.PulseWaveSystolicPeakPos]
        self.DiastolicPhase=self.AverageWaveForm[0][self.PulseWaveEndPos]-self.AverageWaveForm[0][self.PulseWaveSystolicPeakPos]
        self.PulseWaveDuration=self.AverageWaveForm[0][self.PulseWaveEndPos]
        self.PulseWavePropogationTme=self.AverageWaveForm[0][self.DicroticNotchPos]-self.AverageWaveForm[0][self.PulseWaveSystolicPeakPos]
        self.PulseWidthTime=self.HalfTimeHigh-self.HalfTimeLow


        


        
        #help(Data)
    def PlotWaveForm(self):
        x1 = self.AverageWaveForm[0]
        y1 = self.AverageWaveForm[1]
        plt.scatter(x1, y1,label = "PPG Wave Form", color= "red", marker= "*", s=10)
        plt.xlabel('Time [ms]')
        plt.xlim([0,self.PulseWaveDuration])
        plt.ylim([-.5,self.PulseWaveAmplitude])
        plt.annotate('Dicrotic Notch: '+str(int(self.AverageWaveForm[0][self.DicroticNotchPos]))+' ms',
        xy=(self.AverageWaveForm[0][self.DicroticNotchPos], self.AverageWaveForm[1][self.DicroticNotchPos]), xytext=(100, 40),
        textcoords='offset points', ha='right', va='bottom',arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
        plt.annotate('Diastolic Peak: '+str(int(self.AverageWaveForm[0][self.DiastolicPeakPos]))+' ms',
        xy=(self.AverageWaveForm[0][self.DiastolicPeakPos], self.AverageWaveForm[1][self.DiastolicPeakPos]), xytext=(150, 20),
        textcoords='offset points', ha='right', va='bottom',arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
        plt.plot([self.HalfTimeLow, self.HalfTimeHigh], [self.HalfHeight, self.HalfHeight], color='k', linestyle='-', linewidth=2)
        text='Pulse Width'+ str(int(self.PulseWidthTime))+' ms'
        plt.annotate(text,
        xy=(self.HalfTimeLow+(-self.HalfTimeLow+ self.HalfTimeHigh)/2,self.HalfHeight), xytext=(100, -20),
        textcoords='offset points', ha='right', va='bottom',arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
 
        plt.show()


    def PlotPPG(self):
        x1 = self.wf1[0]
        y1 = self.wf1[1]
        plt.scatter(x1, y1,label = "Pulse Wave Form", color= "red", marker= "*", s=1)
        x2 = self.wfCP[0]
        y2 = self.wfCP[1]
        plt.scatter(x2, y2,label = "Critical Points", color= "green", marker= "*", s=20)
        x3 = self.wfD[0]
        y3 = self.wfD[1]
        plt.plot(x3, y3,label = "Derivative", color= "red")
        plt.xlabel('Time [ms]')
        plt.ylabel('Pulse Waveform [Arb. Units]')
        plt.legend(loc='lower right')
        plt.title('PGP Critical Points')
        plt.xlim([5000,18000])
        plt.ylim([-.1,1])
        plt.show()

    def __init__(self,PPG):
        self.make(PPG)

    def getPWBegin(self):
        return self.PWF[self.getPWBegin]

    def getPWSystolicP(self):
        return self.PWF[self.PWPWSystolicP]

    def getPWDicroticNotch(self):
        return self.PWF[self.PWDicroticNotch]

    def getPWDicroticNotch(self):
        return self.PWF[self.PWDicroticNotch]
    
    def getPWDiastolicPeak(self):
        return self.PWF[self.PWDiastolicPeak]
    
    def getPWEnd(self):
        return self.PWF[self.PWEnd]

    def getPWamplitude(self):
        return self.PWF[self.PWSystolicP]-self.PWF[self.PWBegin]
        
    def getRiseTime(self):
        return self.PWF[self.PWSystolicP]-self.PWF[self.PWBegin]

    def getPWDuration(self):
        return self.PWF[self.PWEnd]-self.PWF[self.PWBegin]

    def getPWPropTime(self):
        return self.PWF[self.PWDiastolicPeak]-self.PWF[self.PWSystolicP]

    def PrintAllParameters(self):
        print 'Pulse Wave Systolic Peak Position= ', self.PulseWaveSystolicPeakPos
        print 'Pulse Wave End Position= ', self.PulseWaveEndPos
        print 'Dicrotic Notch Position Position= ', self.DicroticNotchPos
        print 'Diastolic Peak Position= ',self.DiastolicPeakPos
        print  'Half Height= ', self.HalfHeight
        print  'Half Time Low= ', self.HalfTimeLow
        print  'HalfTimeHigh= ', self.HalfTimeHigh                      
        print  'Pulse Wave Amplitude= ', self.PulseWaveAmplitude
        print  'Systolic Phase= ', self.SystolicPhase, ' ms'
        print  'Rise Time= ', self.RiseTime, ' ms'
        print  'Diastolic Phase= ', self.DiastolicPhase, ' ms'
        print  'Pulse Wave Duration= ', self.PulseWaveDuration, ' ms'
        print  'Interbeat Interval= ', self.PulseWaveDuration, ' ms'
        print  'PulseWavePropogationTme= ', self.PulseWavePropogationTme, ' ms'
        print  'PulseWidthTime= ', self.PulseWidthTime, ' ms'

#plotly.tools.set_credentials_file(username='rhcox', api_key='OWm9FUQ3FgSdgiuKs8Cq')
#p1=PPGWave.PPGWave()
#wf=PPGSeg(p1)
#p1.Spectrum()
#wf.PrintAllParameters()
