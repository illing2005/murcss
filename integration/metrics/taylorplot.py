'''
Created on 13.12.2012

@author: Sebastian Illing

Copyright (C) 2014  Sebastian Illing
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>
'''
import numpy as NP
import matplotlib.pyplot as plt 
import abc
import os
from cdo import *
cdo = Cdo()

from filehandler import FileHandler

class TaylorPlotAbstract(object):
    '''
    classdocs
    '''
    __metaclass__ = abc.ABCMeta
    _tmpFiles = []
    _modelValues = {}
    colorsMarker = (['red', 'green', 'blue', 'yellow', 'black', 'grey'], ['s', 'd', '*', 'D'])

    def __init__(self, refContour = False, negativeCorr = False, *args, **kwargs):
        '''
        Constructor
        '''
        self.refStd = 1
        
        self.refContour = refContour
        
        self._modelValues = {}
        
        if(negativeCorr):
            self.width = NP.pi 
            self.labelDirection = "bottom"
            self.labelPosition = "right"
        else:
            self.width = NP.pi / 2
            self.labelDirection = "left"
            self.labelPosition = "center"
        
    def generateGrid(self, fig=None, srange=(0,1.5), rect=111):
        from matplotlib.projections import PolarAxes
        import mpl_toolkits.axisartist.floating_axes as FA
        import mpl_toolkits.axisartist.grid_finder as GF
        
        tr = PolarAxes.PolarTransform()
        
        # Correlation labels
        rlocs = NP.concatenate((NP.arange(10)/10.,[0.95,0.99],-NP.arange(10)/10., [-0.95,-0.99]))
        tlocs = NP.arccos(rlocs)        # Conversion to polar angles
        gl1 = GF.FixedLocator(tlocs)    # Positions
        tf1 = GF.DictFormatter(dict(zip(tlocs, map(str,rlocs))))
        
        # Standard deviation axis extent (in units of reference stddev)
        self.smin = srange[0]*self.refStd
        self.smax = srange[1]*self.refStd

        steps = NP.round((self.smax - self.smin) / 10,1)
        xticks = NP.arange(self.smin, self.smax, steps)
        tf2 = GF.DictFormatter(dict(zip(xticks, map(str, xticks))))
        gl2= GF.FixedLocator(xticks)  
        
        ghelper = FA.GridHelperCurveLinear(tr,
                                           extremes=(0, self.width, 
                                           self.smin, self.smax),
                                           grid_locator1=gl1,
                                           tick_formatter1=tf1,
                                           tick_formatter2=tf2,
                                           grid_locator2=gl2)

        if fig is None:
            fig = plt.figure(figsize=(10, 10))
        
        self.figure = fig
        
        ax = FA.FloatingSubplot(fig, rect, grid_helper=ghelper)
        fig.add_subplot(ax)
        
        # Adjust axes
        ax.axis["top"].set_axis_direction("bottom")  # "Angle axis"
        ax.axis["top"].toggle(ticklabels=True, label=True)
        ax.axis["top"].major_ticklabels.set_axis_direction("top")
        ax.axis["top"].label.set_axis_direction("top")
        ax.axis["top"].label.set_text("Correlation")

        ax.axis["left"].set_axis_direction("bottom") # "X axis"
        ax.axis["left"].label.set_text("Standard deviation")
        ax.axis["left"].label.set_ha(self.labelPosition)

        ax.axis["right"].set_axis_direction("top")   # "Y axis"
        ax.axis["right"].toggle(ticklabels=True)
        ax.axis["right"].major_ticklabels.set_axis_direction(self.labelDirection)
        

        ax.axis["bottom"].set_visible(False)         # Useless
       
        self._ax = ax                   # Graphical axes
        self.ax = ax.get_aux_axes(tr)   # Polar coordinates

        self._add_contours(colors = '0.5')
                    
        plt.title("Taylor diagram") # Figure title
        
        return fig
        
    def showPlot(self):
        plt.show()

    def _add_contours(self, levels=5, **kwargs):
        """Add constant centered RMS difference contours, defined by
        *levels*."""

        rs,ts = NP.meshgrid(NP.linspace(self.smin,self.smax),
                            NP.linspace(0,self.width))

        # Compute centered RMS difference
        rms = NP.sqrt(self.refStd**2 + rs**2 - 2*self.refStd*rs*NP.cos(ts))
       
        contours = self.ax.contour(ts, rs, rms, levels, **kwargs)
        plt.clabel(contours, inline=2, fontsize=10)
        
        return contours
    
    def _addReferencePoint(self):
        #Add reference point and stddev contour
        label='Reference'
        l, = self.ax.plot([0], self.refStd, 'k*', ls='', ms=10, label=label)
        
        #Collect sample points for latter use (e.g. legend)
        self.samplePoints = {'Reference': l}
        #Add contour at self.reStd
        if(self.refContour):
            t = NP.linspace(0, self.width)
            r = NP.zeros_like(t) + self.refStd
            self.ax.plot(t,r, 'k--', label='_')
    
    @abc.abstractmethod
    def _addSample(self, *args, **kwargs):
        """
        Abstract method to add sample values to the plot 
        """
        i=0
        j=0
        for mName in sorted(self._modelValues.iterkeys()):
            l, = self.ax.plot(NP.arccos(self._modelValues[mName]['corr']), self._modelValues[mName]['std'] / self._modelValues[mName]['ref'], 
                              marker= self.colorsMarker[1][j], color= self.colorsMarker[0][i], label = mName, *args, **kwargs) # (theta,radius)
            self.samplePoints.update({mName:l})
            
            i += 1
            if(i == len(self.colorsMarker[0])):
                j += 1
                i = 0
                if(j == len(self.colorsMarker[1])):
                    j = 0
            
    def _addLegend(self):
        """
        Abstract method to add a legend
        """
        # Add a figure legend and title
        #print self.samplePoints
        self._ax.legend([ p for k,p in self.samplePoints.items()],
                   [ k for k,p in self.samplePoints.items() ],
                   numpoints=1, prop=dict(size='small'), loc='upper right')
        
    @abc.abstractmethod
    def _computeValues(self):
        """
        Abstract method to calculate Std and Correlation
        """
        
    def _savePlot(self):
        """
        Abstract method to save the taylor plots
        """
        plt.savefig('taylor.eps', dpi=500, format="eps")
        
        
    @abc.abstractmethod
    def constructPlot(self):
        """
        Abstract method to construct the plot 
        """
        self._addSample()
        self._addLegend()
        if(self.config_dict['displayGraph']):
            self.showPlot()

        

class TaylorPlotMurCSS(TaylorPlotAbstract):
    
    def __init__(self, *args, **kwargs):
        super(TaylorPlotMurCSS, self).__init__(*args, **kwargs)
        
    def _addSample(self, *args, **kwargs):
        """
        Implementation
        """
        super(TaylorPlotMurCSS, self)._addSample()
        
        
    def _computeValues(self, resultList):
        """
        Implementation
        """
        import re
        for res in resultList:
            #get Correlation
            regex = re.compile(".*(correlation).*")
            corr_list = [m.group(0) for l in res for m in [regex.search(l)] if m]
            #std ratio
            regex = re.compile(".*(std_ratio).*")
            std_list = [m.group(0) for l in res for m in [regex.search(l)] if m]
            
            #add to model values
            for i in range(0,len(std_list)):
                name = std_list[i].split('/')[-1]
                name = name.split('_')
                name = name[0] + name[1] + name[2]
                tmp_dict = dict(corr=FileHandler.openNetCDFFile(corr_list[i], 'var'), 
                                std=1/FileHandler.openNetCDFFile(std_list[i], 'var'),
                                ref=1)
                self._modelValues[name] = tmp_dict
                
                
    def _savePlot(self, outputFolder):        
            
        plt.savefig(outputFolder+'taylor.eps', dpi=500, format="eps")    
        
            
                    

    def constructPlot(self, result, outputFolder):
        """
        Abstract method to construct the plot 
        """
        self.generateGrid()
        self._addReferencePoint()
        #self._computeValues(result)
        #self._addSample()
        
        self._addLegend()
        self._savePlot(outputFolder)
        
if __name__ == '__main__':

    tmp = TaylorPlotMurCSS(negativeCorr=True)
    figure = tmp.generateGrid(rect = 111)
    #tmp.generateGrid(rect = 222, fig = figure)
    plt.show()      
    import matplotlib
    print matplotlib.__version__  
