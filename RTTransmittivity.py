import time
import numpy as np
import pandas as pd
import scipy as sp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import sys
import os
from matplotlib.figure import Figure

class Analysis:
    def __init__(self, fileName1=None, fileName2=None, thresh=0.03, dataStart=342, dataEnd=1298, offset=1e-9, p0=[0.02, 1e-10, -5e-10, 5e5, 5e-4, 0]):
        self.reference = True # determines self.fitPulse() behaviour

        self.fileName1 = fileName1
        self.fileName1Tail = None # tail because it's the second part of the filename before it gets concatenated together into a proper file
        self.fileName2 = fileName2
        self.fileName2Tail = None # tail because it's the second part of the filename before it gets concatenated together into a proper file
        self.dataStart = dataStart # starting index of data
        self.dataEnd = dataEnd # end index of data
        self.thresh = thresh # minimum value for pulse to be considered peak
        self.offset = offset # deals with drift 

        self.x = None # time axis
        self.y = None # voltage axis

        # used for method rightleft2
        self.amp = None
        self.b = None
        self.c1 = None
        self.lam = None
        self.d1 = None
        self.d2 = None

        self.h1 = None
        self.h2 = None
        # peak() outputs
        self.peakNum = None
        self.peakPos = None
        self.peakVal = None
        
        # bases() outputs
        self.baseNum = None
        self.basePos = None
        self.baseVal = None

        self.baseNum2 = None
        self.basePos2 = None
        self.baseVal2 = None

        # to be used to plot sets of 4 peaks
        self.every4x = None
        self.every4y = None

        # used as starting point to begin overlapping points
        self.newx = None

        # variables used in superimposing pulse
        self.xTot = np.zeros(0)
        self.yTot = np.zeros(0)
        self.maxx = None

        self.xFit = np.zeros(0)
        self.yFit = np.zeros(0)
        self.popt = None # popt contains all the variables used in rightleft2()
        self.pcov = None
        self.p0 = p0
        
        # fittedX and fittedY will hold all fitted x and y data points for the entire burst respectively
        # xfitpeaks and yfitpeaks hold all locations and amplitudes respectively
        self.fittedX = None
        self.fittedY = None
        self.xFitPeaks = None
        self.yFitPeaks = None

        self.k = None # used in piecewise1() and set in fitAllPulses()

        # used in __finaleReference() and __finaleTransmission
        self.newPeak = None
        self.peakPositionsRefer = None
        self.peakValuesRefer = None
        self.maxPeak = None

        # used in final ratio of transmission to reference plot
        self.peakValues1 = None
        self.peakValues2 = None
        self.peakPositions1 = None
        self.peakPositions2 = None
        self.reflectivity = None 

        self.outputPlots = [None] * 5
        self.xmax = None
        self.curveX = None # for fitted pulses
        self.curveY = None # for fitted pulses

        # these slightly alter the behaviour of self.runAnalysis() 
        self.data1FromArray = False
        self.data2FromArray = False

        # variable that contains x,y values for the ratio plot 
        self.ratioPlotToCSV = []

#  public methods / user methods
    def setFileName1Tail(self, fileNameTail):
        self.fileName1Tail = fileNameTail # sets file number ex. C1trace01234, this method sets the 1234 part

    def setFileName2Tail(self, fileNameTail): # sets file number ex. C1trace01234, this method sets the 1234 part
        self.fileName2Tail = fileNameTail

    def __getData1(self): # reads data from .csv
        filename = '{}{:04d}.csv'.format(self.fileName1, self.fileName1Tail)
        data = pd.read_csv(os.path.join('', filename), skiprows=6).to_numpy().T
        self.x, self.y = data[0,:], data[1,:]

    def getData1FromArray(self, x, y): # read data from array
        self.x, self.y = x, y
        self.data1FromArray = True

    def __getData2(self): # reads data from .csv
        filename = '{}{:04d}.csv'.format(self.fileName2, self.fileName2Tail)
        data = pd.read_csv(os.path.join('', filename), skiprows=6).to_numpy().T
        self.x, self.y = data[0,:], data[1,:]
        # bottom part accounts for "drift" that may occur due to differing BNC connector. 
        for i in range(np.size(self.x)):
            if self.x[i] >= self.x[self.dataStart] + self.offset:
                lowlimit = i
                break
        for i in range(np.size(self.x)):
            if self.x[i] >= self.x[self.dataEnd] + self.offset:
                highlimit = i
                break

        self.x, self.y = self.x[lowlimit:highlimit], self.y[lowlimit:highlimit]

    def getData2FromArray(self, x, y): # reads data from array
        self.x, self.y = x, y
        for i in range(np.size(self.x)):
            if self.x[i] >= self.x[self.dataStart] + self.offset:
                lowlimit = i 
                break
        for i in range(np.size(self.x)):
            if self.x[i] >= self.x[self.dataEnd] + self.offset:
                highlimit = i 
                break
        
        self.x, self.y = self.x[lowlimit:highlimit], self.y[lowlimit:highlimit]
        self.data2FromArray = True

    def __grabPeakValues(self, reference): # reference should be a boolean value, it determines whether or not these peaks are for reference or transmission photodiode
        # gets the number of peaks, time position of each peak, and amplitude of each peak
        if (reference not in [True, False]):
            raise Exception("Reference should be boolean")
        self.reference = reference
        if (reference):    
            self.peakNum, self.peakPos, self.peakVal = self.__peaks(self.x[self.dataStart:self.dataEnd], self.y[self.dataStart:self.dataEnd], self.thresh,updown=6)
        else:
            self.peakNum, self.peakPos, self.peakVal = self.__peaks(self.x[0::], self.y[0::], self.thresh, updown=8)

    def __setEveryFourth(self):
        # creates an array of every fourth peak
        self.every4x = self.__everyFourth(self.peakPos[0:self.peakPos.size])
        self.every4y = self.__everyFourth(self.peakVal[0:self.peakVal.size])
    
    # allows user to set reference photodiode parameters
    def setPhotodiodeParams(self, amp=None, b=None, c1=None, lam=None, d1=None, d2=None):
        if amp is not None:
            self.amp = amp
        if b is not None:
            self.b = b
        if c1 is not None:
            self.c1 = c1
        if lam is not None:
            self.lam = lam
        if d1 is not None:
            self.d1 = d1
        if d2 is not None:
            self.d2 = d2





    def __rightleft2(self, x, amp, b, c1, lam, d1, d2):
        h1 = amp-d1
        t = 200 # should change
        h2 = (h1*np.exp(-(-t*10**(-12))**2/(2*c1**2))+d1)-d2
        y = (h1*np.exp(-(x-b)**2/(2*c1**2))+d1)*(np.heaviside(-(x-b-t*10**(-12)), amp/2)) + (h2*np.exp(-(x-b-t*10**(-12))*lam) + d2)*(np.heaviside(x-b-t*10**(-12), amp/2))
        return y
    
    #THIS FUNCTION OUTPUTS THE NUMBER OF PEAKS AND THEIR X AND Y VALUES
    def __peaks(self, x,y, threshold, updown, maxpeaks=100000, xpos=-1, timespace=3e-9):
        count = 0
        peakpos = np.zeros(0)
        peakval = np.zeros(0)
        for i in range(x.size):
            
            true = 1 # don't like this
            
            for j in range(updown):
                if (i-j) < 0 or (i+j)>x.size-1 or y[i]<threshold:
                    true = 0
                elif (y[i]<y[i-j] or y[i]<y[i+j]):
                    true = 0
                    
                        
            if y[i] == y[i-1] or count == maxpeaks or x[i]<xpos:
                true = 0    
            if count != 0 and x[i]<peakpos[count-1]+timespace: 
                        true = 0
            if true == 1:
                count += 1
                peakpos = np.append(peakpos, x[i])
                peakval = np.append(peakval, y[i])
        return (count, peakpos, peakval)

    # Outputs every fourth peak in burst
    def __everyFourth(self, arr):
        cycle = 0
        every4 = np.zeros(0)
        for i in range(arr.size):
            if cycle == 0: # should change
                every4 = np.append(every4, arr[i])     
            if cycle == 3:
                cycle = 0
            else:
                cycle = cycle+1
        return (every4)   
    # specifies the exact point to start overlapping points (artificial high sample)
    def __calculateNewx(self):
        self.newx = 0
        for k in range(self.y.size):
            if self.y[k] == max(self.every4y):
                for h in range(self.peakPos.size):
                    if self.x[k] == self.peakPos[h]:
                        self.newx = self.x[k]  
    # finds average pulse shape
    def __pulseAverageShape(self):
        peakxVal = np.zeros(10)
        self.xTot = np.zeros(0)
        self.yTot = np.zeros(0)
        for i in range(peakxVal.size):
            peakxVal[i] = (self.newx)-i*18.8835616438*10**(-9)
            datax = np.zeros(0)
            datay = np.zeros(0) 
            for j in range(self.x.size):
                if np.abs(peakxVal[i]-self.x[j])<2.6*10**(-9):
                    datax = np.append(datax, self.x[j])
                    datay = np.append(datay, self.y[j])
            datax = datax+i*18.8835616438*10**(-9) 
            self.xTot = np.append(self.xTot, datax)
            self.yTot = np.append(self.yTot, datay)
            if peakxVal[i] == 0:
                print(i)
    # arranges averaged pulses in ascending order
    def __arrangeSamplesLowHigh(self):
        for i in range(self.xTot.size):
            lowest = 100
            for j in range(self.xTot.size-i):
                h = j+i
                if self.xTot[h]<lowest:
                    lowest = self.xTot[h]
                    num3 = h
                lowesty = self.yTot[num3] 
            self.xTot[num3] = self.xTot[i]
            self.yTot[num3] = self.yTot[i]
            self.xTot[i] = lowest
            self.yTot[i] = lowesty
        for i in range(self.yTot.size):
            if self.yTot[i] == max(self.yTot):
                self.maxx = self.xTot[i]
    # fits the superimposed pulse created from 50 pulses
    def __fitPulse(self):
        for i in range(self.xTot.size):
            if -2*10**(-9) < (self.xTot[i]-self.maxx) < 2.5*10**(-9): #should change
                self.xFit = np.append(self.xFit, self.xTot[i])
                self.yFit = np.append(self.yFit, self.yTot[i])

        self.xFit = self.xFit - self.maxx
        if (self.reference):
            self.popt, self.pcov = curve_fit(self.__rightleft2, self.xFit, self.yFit, sigma=np.zeros(self.xFit.size)+0.0005, absolute_sigma=True, p0=self.p0)
    # updates h1 and h2 in preperation for fitting each of the 50 pulses individually
    def __updateH(self):
        self.h1 = self.amp-self.d1
        self.h2 = (self.h1*np.exp(-(-200*10**(-12))**2/(2*self.c1**2))+self.d1)-self.d2  # should change 200
    # piecewise function used in the fitting of each pulse
    def __piecewise1(self, x, a, b):
        if self.reference:
            y = a*((self.h1*np.exp(-(x-b)**2/(2*self.c1**2)))*(np.heaviside(-(x-b-200*10**(-12)),(self.h2+self.d2-self.d1)/2)) + (self.h2*np.exp(-(x-b-200*10**(-12))*self.lam) + self.d2-self.d1)*(np.heaviside(x-b-200*10**(-12), (self.h2+self.d2-self.d1)/2)))+self.baseVal[self.k]
        elif not self.reference:
            y = a*((self.h1*np.exp(-(x-b)**2/(2*self.c1**2)))*(np.heaviside(-(x-b-200*10**(-12)),(self.h2+self.d2-self.d1)/2)) + (self.h2*np.exp(-(x-b-200*10**(-12))*self.lam) + self.d2-self.d1)*(np.heaviside(x-b-200*10**(-12), (self.h2+self.d2-self.d1)/2)))+self.baseVal2[self.k]
        return y  # should change 200

    # modified peaks()
    def __bases(self, x, y, updown, maxPeaks=1e5, xPos=-1, timeSpace=3e-9):
        count = 0
        peakPos = np.zeros(0)
        peakVal = np.zeros(0)
        for i in range(x.size):
            true = 1 # don't like this, but im a lil eepy right now and don't want to think variable names and this is the best i could think of
            for j in range(updown):
                if (i - j) < 0 or (i + j) > x.size - 1:
                    true = 0  
                elif y[i] > y[i - j] or y[i] > y[i + j]:
                    true = 0  
            if y[i] == y[i-1] or count == maxPeaks or x[i] < xPos:
                true = 0
            if (count != 0 and x[i] < peakPos[count-1] + timeSpace):
                true = 0
            if true == 1:
                count += 1
                peakPos = np.append(peakPos, x[i])
                peakVal = np.append(peakVal, y[i])
        return (count, peakPos, peakVal)
    # filters out useless oscilloscope data
    def __sliceXY(self):
        self.x=self.x[self.dataStart:self.dataEnd]
        self.y=self.y[self.dataStart:self.dataEnd]
    # fits each of the 50 pulses individually
    def __fitAllPulse(self):
        self.fittedX = np.zeros(0)
        self.fittedY = np.zeros(0)
        self.xFitPeaks = np.zeros(0)
        self.yFitPeaks = np.zeros(0)
        self.curveX = []
        self.curveY = []
        for k in range(self.peakPos.size): 

            self.k = k
            fitDataX = np.zeros(0)
            fitDataY = np.zeros(0)
            ySmall = 100
            for l in range(self.x.size):
                if (-0.5e-9 < self.x[l] - self.peakPos[k] < 0) and self.y[l] <= ySmall:
                    ysmall = self.y[l]
                    minx = self.x[l]

            for l in range(self.x.size):
                if (minx+0.25e-9 <= self.x[l] < (0.8e-9+self.peakPos[k])):
                    fitDataX = np.append(fitDataX, self.x[l])
                    fitDataY = np.append(fitDataY, self.y[l])

            # First round of fitting using highest point as estimated peak
            fitDataX = fitDataX - self.peakPos[k]
            self.popt, self.pcov = curve_fit(self.__piecewise1, fitDataX, fitDataY, 
                                             sigma=np.zeros(fitDataX.size)+8.5e-4, 
                                             absolute_sigma=True,
                                             p0=[1,1e-10])
            self.xFit = np.linspace(minx-self.peakPos[k]+0.25e-9, 0.8e-9, 100)
            self.yFit = self.__piecewise1(self.xFit, self.popt[0], self.popt[1])
            xmax = 0
            for l in range(self.xFit.size):
                if self.yFit[l] == max(self.yFit):
                    xmax = self.xFit[l] + self.peakPos[k]
            
            # second round of fitting

            fitDataX = np.zeros(0)
            fitDataY = np.zeros(0)
            
            for l in range(self.x.size):
                if (-0.5e-9 <= (self.x[l]-xmax) < 0.8e-9):
                    fitDataX = np.append(fitDataX, self.x[l])
                    fitDataY = np.append(fitDataY, self.y[l])
            fitDataX = fitDataX - xmax

            self.popt, self.pcov = curve_fit(self.__piecewise1, fitDataX, fitDataY, sigma=np.zeros(fitDataX.size)+8.5e-4,
                                             absolute_sigma=True, p0=[1, 1e-10])
            self.xFit = np.linspace(-2e-9, 2.5e-9, 100)
            self.yFit = self.__piecewise1(self.xFit, self.popt[0], self.popt[1])
            true = 1
            for l in range(self.xFit.size):
                if self.yFit[l] == max(self.yFit) and true == 1:
                    self.xFitPeaks = np.append(self.xFitPeaks, self.xFit[l])
                    true = 0
            self.yFitPeaks = np.append(self.yFitPeaks, max(self.yFit))

            self.fittedX = np.append(self.fittedX, self.xFit+xmax)
            self.fittedY = np.append(self.fittedY, self.yFit)
            
            self.curveX.append(self.xFit+xmax)
            self.curveY.append(self.yFit)
            #plt.plot(self.x, self.y, 'ro', markersize=2)
            #plt.plot(self.xFit+xmax, self.yFit, 'g-') 
    # corrects any potential mismatch between highest amplitude fitted pulse and 
    # highest voltage from raw voltage data
    def __finaleReference(self): # not sure what to name this method
        self.newPeak = self.__peaks(self.fittedX, self.fittedY, self.thresh, 5)
        self.peakPositionsRefer = self.newPeak[1]
        self.peakValuesRefer = self.newPeak[2]

        for i in range(self.peakValuesRefer.size):
            if self.peakValuesRefer[i] == max(self.peakValuesRefer):
                self.maxPeak = self.peakPositionsRefer[i]

        self.peakPositions1 = self.peakPositionsRefer[0::]

        if self.peakPositions1[0]>self.basePos[0]:
            self.peakValues1 = self.newPeak[2] - self.baseVal
        else:
            raise Exception("Mismatched base and peak position")
    # corrects any potential mismatch between highest amplitude fitted pulse and 
    # highest voltage from raw voltage data
    # also calculates the ratio plot
    def __finaleTransmission(self):
        self.newPeak = self.__peaks(self.fittedX, self.fittedY, self.thresh, 6)
        self.peakPositions2 = self.newPeak[1]
        
        if self.peakPositions2[0] > self.basePos2[0]:
            self.peakValues2 = self.newPeak[2] - self.baseVal2
        else:
            raise Exception("Mismatched base and peak position")
        
        self.reflectivity=(self.peakValues2)/self.peakValues1

    # main function responsible for running the pulse fitting algorithm in order
    # also generates output plots
    # commented out plt methods were used in testing. Leaving them in in case this file is being used to analyze a .csv and not
    # data directly from an oscilloscope. Do not uncomment the plots when program is being used to analyze data stright from an
    # oscilloscope. It will cause errors as matplotlib needs to be run on the main thread and not any other threads. 
    def runAnalysis(self, data1Array=None, data2Array=None): 
        # selects which get data method depending on input data source .csv vs direct from usb
        if data1Array is not None:
            self.getData1FromArray(data1Array[0], data1Array[1])
        else:
            self.__getData1()

        # creates raw data plot for output
        figureRef = Figure(figsize=(5,4), dpi=100)
        ax1 = figureRef.add_subplot(111)
        ax1.plot(self.x, self.y, linewidth=1)
        ax1.plot(self.x, self.y, 'o', markersize=1)
        ax1.set_title("Reference Traces")
        self.outputPlots[0] = figureRef

        # plt.plot(self.x, self.y,linewidth=1) # this is just used in testing
        # plt.plot(self.x,self.y,'o', markersize=1)
        # plt.title("this is the first plot")
        # plt.show()

        # starting pulse fitting algorithm here
        self.__grabPeakValues(reference=True)
        self.__setEveryFourth()
        self.__calculateNewx()
        self.__pulseAverageShape()
        self.__arrangeSamplesLowHigh()

        # fig,ax=plt.subplots()
        # ax.plot(self.xTot, self.yTot,'o')
        # plt.show() 
        self.__fitPulse() # fitting superimposed pulse

        # if photodiode parameters not set then the algorithm continues on with calculated values
        if all(v is None for v in (self.amp,self.c1, self.d1,self.lam,self.d2)):
            self.setPhotodiodeParams(amp=self.popt[0],b=self.popt[1], c1=self.popt[2], d1=self.popt[4],lam=self.popt[3], d2=self.popt[5])
        
        # plt.plot(self.xFit,self.yFit,'ro', markersize=1)
        # xplot = np.linspace(-2*10**(-9), 2.8*10**(-9), 100)
        # yplot = self.rightleft2(xplot, self.amp,self.b,self.c1,self.lam,self.d1,self.d2)
        # plt.plot(xplot,yplot)
        # plt.show()

        # prepares first part fitted pulses (only dots/oscilloscope samples)
        figureRefFitted = Figure(figsize=(5,4), dpi=100)
        ax2 = figureRefFitted.add_subplot(111)
        ax2.plot(self.x[self.dataStart:self.dataEnd], self.y[self.dataStart:self.dataEnd], "ro", markersize = 1)

        self.__updateH()
        self.peakNum, self.peakPos, self.peakVal = self.__peaks(self.x[self.dataStart:self.dataEnd], self.y[self.dataStart:self.dataEnd], self.thresh, updown=7)
        self.baseNum, self.basePos, self.baseVal = self.__bases(self.x[self.dataStart:self.dataEnd], self.y[self.dataStart:self.dataEnd], updown=7)
        self.__sliceXY()
        self.__fitAllPulse()

        # prepares second part fitted pulses (the fitted functions)
        for i in range(50): # for plotting each of the 50 fitted curves
            ax2.plot(self.curveX[i], self.curveY[i], 'g-')
        ax2.set_title("Fitted Transmission Pulses")
        self.outputPlots[1] = figureRefFitted

        self.__finaleReference() 
        
        # part 2 for transmission photodiode
        # selects get data method depending on data input source .csv vs USB
        if data2Array is not None:
            self.getData2FromArray(data2Array[0], data2Array[1])
        else:
            self.__getData2()

        # creates raw data plot for output
        figureTransmission = Figure(figsize=(5,4), dpi=100)
        ax3 = figureTransmission.add_subplot(111)
        ax3.plot(self.x, self.y, linewidth=1)
        ax3.plot(self.x, self.y, 'o', markersize=1)
        ax3.set_title("Signal Traces")
        self.outputPlots[2] = figureTransmission

        # plt.plot(self.x, self.y,linewidth=1)
        # plt.plot(self.x,self.y,'o', markersize=1)
        # plt.title("this is the first plot")
        # plt.show() 

        # begins pulse fitting algorithm for transmission signal photodiode
        self.__grabPeakValues(reference=False)
        self.__setEveryFourth()
        self.__calculateNewx()
        self.__pulseAverageShape()
        self.__arrangeSamplesLowHigh()

        # fig,ax=plt.subplots()
        # ax.plot(self.xTot, self.yTot, 'o')
        # plt.show()
        
        self.__fitPulse() # fits superimposed pulse

        # if user does not set transmission photodiode parameters then algorithm continues on with calculated values
        if all(v is None for v in (self.amp,self.c1, self.d1,self.lam,self.d2)):
            self.setPhotodiodeParams(amp=self.popt[0],b=self.popt[1], c1=self.popt[2], d1=self.popt[4],lam=self.popt[3], d2=self.popt[5])

        # prepares first part fitted pulses (only dots/oscilloscope samples)
        figureTransmissionFitted = Figure(figsize=(5,4), dpi=100)
        ax4 = figureTransmissionFitted.add_subplot(111)
        ax4.plot(self.x, self.y, "ro", markersize=1)

        # fits each of the 50 pulses for transmitted signal
        self.__updateH()
        self.peakNum, self.peakPos, self.peakVal = self.__peaks(self.x, self.y, self.thresh, updown=6)
        self.baseNum2, self.basePos2, self.baseVal2 = self.__bases(self.x, self.y, updown=6)
        self.__fitAllPulse()
        self.__finaleTransmission()


        # prepares second part fitted pulses (the fitted functions/curves)
        for i in range(50):  # for plotting each of the 50 fitted curves
            ax4.plot(self.curveX[i], self.curveY[i], "g-")
        ax4.set_title("Fitted Signal Pulses")
        self.outputPlots[3] = figureTransmissionFitted
        figureRatio = Figure(figsize=(5,4), dpi=100)
        ax5 = figureRatio.add_subplot(111)
        ax5.plot(self.peakPositions1, self.reflectivity)
        ax5.set_title("Ratio")

        # prepares ratio plot
        self.outputPlots[4] = figureRatio

        # output for saving ratio plot as .csv file via GUI class save button
        self.ratioPlotToCSV[0] = self.peakPositions1
        self.ratioPlotToCSV[1] = self.reflectivity
        # plt.plot(self.peakPositions1, self.reflectivity)
        # plt.show() 
        

# just testing
if __name__ == "__main__":
    analysis = Analysis('C1trace0', 'C3trace0',thresh=0.01, dataStart=343)
    analysis.setFileName1Tail(38)
    analysis.setFileName2Tail(38)
    analysis.runAnalysis()