'''
Created on 11.03.2013

:author: Sebastian Illing
:contact: sebastian.illing@met.fu-berlin.de

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

import abc
import os
from cdo import *
#from integration.metrics.filehandler import FileHandler
cdo = Cdo()
import time

import numpy as N
import sys

import murcss_config
if murcss_config.file_system == 'miklip':
    from findFiles import FindFiles
    from findFilesSeason import FindFilesSeason
else:
    from findFilesCustom import FindFilesCustom as FindFiles
    FindFilesSeason = FindFiles

from filehandler import FileHandler
from plotter import Plotter
from nondeamonpool import MyPool
from tool_abstract import ToolAbstract#, unwrap_self_f

class MetricAbstract(ToolAbstract):
    '''
    Abstract class for calculating any metric of decadal runs with CDO
    Every reusable function which uses CDO commands should go in here.
    '''
    
    def __init__(self, 
                 tmpDir = './tmp', 
                 output = '/', 
                 output_plots='/', 
                 baseDir='/', 
                 result_grid=None,
                 observation='HadCrut',
                 leadtimes_mode='yearly',
                 decadals='1960,1965,1970,1975,1980,1985,1990,1995,2000',
                 months=None
                ):
        '''
        Constructor. Sets variables and creates Folders
         
        :param tmpDir: 
        :param output: 
        :param output_plots: self.outputPlots
        :param baseDir: Path to search for source filesself.outputPlots
        :param result_grid: grid to remap to
        :param observation: observation name or file used for analysis
        '''    
        setattr(self, 'obsExp', observation) #TODO: Wrong named kwarg
        self.observation = ''    
        if(self.obsExp == 'HadCrut' and self.variable == 'tas'):  
            self.observation = baseDir + '/../src/obs2/tas_obs_HadCrut3.nc'
        elif os.path.isfile(self.obsExp):
            self.observation = self.obsExp
            self.obsExp = self.extractFilename(self.obsExp)
            self.obsExp = self.obsExp.split('.')[0]
        
        
        #Initialize MaxLeadYear
        self.getYearRange()        
        outputFolder = self.constructName(self.outputFolderStruct)
        super(MetricAbstract, self).__init__(output_folder = self.checkPath(output)+outputFolder, 
                                             output_plots = self.checkPath(output_plots)+outputFolder)
        self.outputDir = self.checkPath(self.checkPath(output)+outputFolder)
        self.outputPlots = self.checkPath(self.checkPath(output_plots)+outputFolder)
        self.baseDir = baseDir
        self.tmpDir = self.checkPath(tmpDir)
        if not os.path.isdir(self.checkPath(tmpDir)):
            os.makedirs(self.checkPath(tmpDir))
        
        if result_grid is None: 
            self.gridFile = self.baseDir + '/../src/grids/griddes_HadCrut3_5x5.txt'
        else:
            self.gridFile = result_grid         
        
        
        if leadtimes_mode == 'monthly':
            if months is not None:
                month_list = months.split(',')
            else:
                month_list = ['']
            init_list = list()
            for month in month_list:
                for year in decadals:
                    init_list.append(int(str(year)+month))
            self.decadals = init_list
        else:
            self.decadals=decadals
        #set special options for "decadal" or "seasonal"
        self.leadtimes_mode = leadtimes_mode
        if self.leadtimes_mode == 'yearly':
            self.analysisOptions = {'cdoSelect' : 'selyear',
                                    'cdomean' : 'yearmean'}
        else:
            self.analysisOptions = {'cdoSelect' : 'selmon',
                                    'cdomean' : 'monmean'}
                
        #findFiles instance
        if leadtimes_mode == 'yearly':
            self.findFiles = FindFiles(tmpDir=self.tmpDir, level=self.level)
        elif leadtimes_mode == 'monthly':
            self.findFiles = FindFilesSeason(tmpDir=self.tmpDir, level=self.level)
        
    @abc.abstractmethod    
    def analyze(self):
        '''
        Abstract Method to calculate a metric
        Must be implemented in child classes
        '''
        pass
        
    def _remapFiles(self, ensList, flag=''):
        '''
        Multiprocessing for remapping files
        remap observations and model files
        
        :param ensList: list of files
        :return: list of remapped files
        '''
        
        ensCount = len(ensList)
        poolArgs = self.getPoolArgs(ensCount,self,ensList,flag,'remapFile')
        return self.multiProcess(poolArgs)

    def remapFile(self, *arg, **kwargs): 
        '''
        Remaps a single file width CDO's remapcon
        
        :param inputfile 
        '''
        iFile= arg  
        iFile=iFile[0]
        
        try:
            flag = arg[1]
            if flag == 'output*':
                flag = ''
        except:
            flag = ''
        flag += self.getRandomStr()
        fName = iFile.split('/')[-1]
        return cdo.remapcon(self.gridFile, input=iFile, output=self.tmpDir+fName+flag+'.remapped')
    
    def getEnsembleMean(self, fileList):
        '''
        :todo: Multiprocessing
        Should do ensmean in multiprocessing!!! At the moment only single approach
        
        :param fileList: dict of filelists
        :return: dict of ensemble means
        '''
        ensMeanList = self.yeardictToList(fileList)
            
        count = len(ensMeanList)
        poolArgs = self.getPoolArgs(count,self,ensMeanList,'_getEnsembleMean')
        return self.listToYeardict(self.multiProcess(poolArgs))

    
    def _getEnsembleMean(self, fileList, flag='ENSMEAN'):
        '''
        Calculates ensemble Mean
        
        :param fileList: list of files
        :return: ensemble mean file 
        '''
        if len(fileList) == 1:
            return fileList[0]
        fileStr = ' '.join(fileList)
        return cdo.ensmean(input=fileStr, output=self.tmpDir+self.extractFilename(fileList[0])+'_'+flag) 
    
    def tempSmoothing(self, ensList, startyear, endyear):
        '''
        Multiprocessing approach for temporal smoothing
        Starts a process for every file in ensList
        
        :param ensList: List of files
        :param startyear: int
        :param endyear: int
        :return: list of temporal smoothed files
        '''
        if(isinstance(ensList,list)):
            ensCount = len(ensList)
        else:
            ensCount = 1
        poolArgs = self.getPoolArgs(ensCount,self,ensList,startyear,endyear,'_tempSmoothing')
        return  self.multiProcess(poolArgs)
    
    def _tempSmoothing(self, fileName, startyear, endyear):
        '''
        Temporal smoothing. I.e. timmean over year 2-9
        
        :param fileName: filepath
        :param startyear: first year to select (int)
        :param endyear: last year to select (int)
        :return: temp smoothed file
        '''  
#        pid= os.getpid()
#        filestartyear = int(cdo.showyear(input=fileName)[0].split(' ')[0])
#        selyearstr = ''
#        fn = fileName.split('/')[-1]
#        fileNameOut = self.tmpDir+fn
#        for i in range(startyear, endyear+1):
#            selyearstr = ','.join([selyearstr, str(filestartyear+i-1)])
#        selyear = cdo.selyear(selyearstr, input=fileName, output=fileNameOut+'_selYear_'+str(startyear)+'-'+str(endyear)+str(pid))
#        return cdo.timmean(input=selyear, output=selyear+'Timmean')
        
        #perform yearmean or monthmean
        #try:
        
        
        pid= os.getpid()
        fn = fileName.split('/')[-1]
        fileNameOut = self.tmpDir+fn
        
        
        #self.analysisOptions['cdomean'] = 'seasmean'
        #method = getattr(cdo, self.analysisOptions['cdomean'])
        #meanFile =  method(input = fileName, output=self.tmpDir+fn+str(pid)+self.analysisOptions['cdomean'])
        
#        #test monthly analysis
#        self.analysisOptions['cdomean'] = 'monmean'
#        method = getattr(cdo, self.analysisOptions['cdomean'])
#        meanFile =  method(input = fileName, output=self.tmpDir+fn+str(pid)+self.analysisOptions['cdomean'])        
#        selstr = ','.join(map(str,range((startyear-1)*12+1,endyear*12+1)))
#        seldate = cdo.seltimestep(selstr,input=meanFile, output=fileNameOut+'_selYear_'+str(startyear)+'-'+str(endyear)+str(pid)+self.getRandomStr())
#        tmp = cdo.ymonmean(input=seldate, output=seldate+'Timmean')
#        return tmp
#        
        #TEST SEASONANALYSIS
#        self.analysisOptions['cdomean'] = 'seasmean'
#        method = getattr(cdo, self.analysisOptions['cdomean'])
#        meanFile =  method(input = fileName, output=self.tmpDir+fn+str(pid)+self.analysisOptions['cdomean'])
#        meanFile = cdo.selseas('JJA',input=meanFile,output=meanFile+'JJA')
#        selstr = ','.join(map(str,range(startyear,endyear+1)))
#        seldate = cdo.seltimestep(selstr,input=meanFile, output=fileNameOut+'_selYear_'+str(startyear)+'-'+str(endyear)+str(pid)+self.getRandomStr())
#        tmp = cdo.timmean(input=seldate, output=seldate+'Timmean')
#        return tmp
    
    
        method = getattr(cdo, self.analysisOptions['cdomean'])
        meanFile =  method(input = fileName, output=self.tmpDir+fn+str(pid)+self.analysisOptions['cdomean'])  
        selstr = ','.join(map(str,range(startyear,endyear+1)))
        seldate = cdo.seltimestep(selstr,input=meanFile, output=fileNameOut+'_selYear_'+str(startyear)+'-'+str(endyear)+str(pid)+self.getRandomStr())
        tmp = cdo.timmean(input=seldate, output=seldate+'Timmean')
        return tmp
        
    def removeSeasonalCycle(self, fn):
        '''
        Subtracts the seasonal cycle of monthly data
        '''
        print fn
        seasonalCycle = cdo.ymonavg(input=fn, output=fn+'_meanCycle')
        print seasonalCycle
        tmp = cdo.ymonsub(input=' '.join([fn,seasonalCycle]),output=fn+'_noseason')
        print tmp
        return tmp
    
    
    
    def getAnomalies(self, ensMeanList, crossMean):
        '''
        Subtracts cross-validated mean from ensembles
        
        :param ensMEanList: dict with filenames
        :param crossMean: dict width cross-val means
        :return: dict width cross-validated anomalies
        '''
        anomalies = dict()
        for year in self.decadals:
            anomalies[year] = cdo.sub(input=' '.join([ensMeanList[year], crossMean[year]]), output=ensMeanList[year]+'_ANOMALIE')
            #
            #anomalies[year] = cdo.sub(input=' '.join([ensMeanList[year], cdo.timmean(input=crossMean[year],output=crossMean[year]+'timmean')]), output=ensMeanList[year]+'_ANOMALIE')        
        return anomalies 
    
    def _getAnomalies(self, ensMeanList, crossMean):
        
        if type(ensMeanList) == list:
            tmpList = list()
            for i,item in enumerate(ensMeanList):
                tmpList.append(self._getAnomalies(item, crossMean))
            return tmpList
        else:
            return cdo.sub(input=' '.join([ensMeanList, crossMean]), output=ensMeanList+'_ANOMALIE')
    
    def getCrossValMean(self, ensMeanList, flag):
        '''
        Calculates cross-validated means (averages through forecast times, excluding the forecast in question)
        
        :param ensMeanList: Dict of filnames
        :param flag: String value added to new filenames
        :return: dict with cross-validated means
        '''
        poolMeanList = list()
        
        yearList = list()
        flagList = list()
        for year in self.decadals:
            yearList.append(year)
            meanList = list()
            flagList.append(str(year)+'_'+flag+self.getRandomStr())
            for key,val in ensMeanList.iteritems():
                if( key != year):
                    meanList.append(val)
        
            poolMeanList.append(meanList)
        
        count = len(poolMeanList)
        poolArgs = self.getPoolArgs(count,self,poolMeanList,flagList,'_getEnsembleMean')
        return self.listToYeardict(self.multiProcess(poolArgs))
            
    
    def getVariance(self, list, flag):
        '''
        Calculate "variance"
        
        :deprecated: Not used in goddard metrics --> should not be used 
        :param list: filelist
        :param flag: string for filenames
        :return: kind of standard deviation of anomalies
        '''
        squareList = dict()
        squareStr = ''
        for year in self.decadals:
            squareList[year] = cdo.sqr(input=list[year], output=list[year]+'_square')
            squareStr += ' '+squareList[year]
            
        squareSum = cdo.ensmean(input=squareStr, output=self.tmpDir+flag+'_goddardVar')
        return cdo.sqrt(input=squareSum, output=self.tmpDir+flag+'_goddardStd')
    
    def createConstantFile(self, gridFile, flag=''):
        '''
        Creates a constant field with value 1
        
        :param gridFile: grid description file (txt)
        :return: filename
        '''
        if self.fieldmean:
            return cdo.const('1,r1x1', output=self.tmpDir+flag+'constantField.nc', options = '-f nc')
        elif self.zonalmean:
            tmp_const = cdo.setrtoc('-1e99,1e99,1',input=self.observationRemapped[self.decadals[0]],output=self.tmpDir+flag+'constantField_tmp.nc',
                                    options = '-f nc')  
            return cdo.timmean(input=tmp_const,output=self.tmpDir+flag+'constantField.nc')      
        else:
            const = cdo.const('1,'+gridFile, output=self.tmpDir+flag+'constantField.nc', options = '-f nc')
            if self.lonlatbox is not None:
                return self._sellonlatbox(const)
            else:
                return const
   
    def getEnsembleStd(self, ensList):
        '''
        Calculates the mean ensemble STD for a given dictionary
        Keys in dict have to be the starting years 
        
        :param ensList: dict with ensembles
        :return: file with mean ensemble std
        '''
        tmpList = list()
        for year in self.decadals:
            tmpList.append(ensList[year])
        ensCount = len(tmpList)   
        poolArgs = self.getPoolArgs(ensCount,self,tmpList,'_getEnsembleVar')
        result = self.listToYeardict(self.multiProcess(poolArgs))
        return cdo.sqrt(input=cdo.ensmean(input=' '.join(result.values()), output=result.values()[0]+'ensmean_afgst'), output=result.values()[0]+'AVGSTD')
            
    def _getEnsembleVar(self, fileList):
        '''
        Calculates the ensemble Variance of a given file list (ensembles)
        
        :note: CDO norally divides by n here this is changed to n-1
        :param fileList: List of files
        :return: file with ensemble variance
        '''
        fileStr = ' '.join(fileList)
        norm = len(fileList)/(len(fileList)-1.)
        ensvar = cdo.ensvar(input=fileStr, output=fileList[0]+'_ENSVAR'+self.getRandomStr())
        #print len(fileList)
        
        return cdo.mulc(norm, input=ensvar, output=ensvar+'_BETTER')
    
    
    def printValues(self, variable, msg):
        '''
        Method for debugging. Prints our values of a list, dict or file.
        
        :note: Printing is only activated for TESTDATA
        
        :todo: Maybe extend for numpy arrays
        
        :param variable: variable to print
        :param msg: name of variable
        '''
        if self.obsExp != 'TESTDATA':
            return
        
        typeName = type(variable).__name__
        if typeName == 'str':
            self.printFileValue(variable, msg)
        elif typeName == 'list':
            for item in variable:
                self.printValues(item, msg)
        elif typeName == 'dict':
            for year in self.decadals:
                self.printValues(variable[year], msg+'_'+str(year))
        else:
            print 'Cant recognize filetype.'         
            
    def printFileValue(self,filename, msg):
        '''
        Method for debugging. Prints a specific value of a given netCDF File.
        
        :param filename:
        :param msg: additional message to print out. Like the name of the variable
        '''
        variable = FileHandler.openNetCDFFile(filename, mode='var')
        print np.shape(variable)
        if len(np.shape(variable)) == 2:
            print 'Value of %s is: %s' % (msg,variable[16,37])
        else:
            print 'Value of %s is: %s' % (msg,variable[0,16,37])
        
    
    def calcMissingValueMask(self, observations):
        '''
        Calculates a missing value mask using the observations. 
        At the moment all gridpoints are masked at missing value where at least one value is missing
        
        :TODO: find a less strict solution. Maybe 10% available?
        DONE!
        
        :param observations: dict with all observations
        :return netcdf file with missing values mask 1/0
        '''

#        for year in self.decadals:
#            
#            tmp_array = FileHandler.openNetCDFFile(observations[year], 'var')
#            tmp_array = np.expand_dims(tmp_array, axis=2)
#            try:
#                master = np.concatenate((master,tmp_array),axis=2)
#            except NameError:
#                master = tmp_array  
#        missing_value_mask = np.zeros((master.shape[0],master.shape[1]))
#        #This loop is probably not the best idea
#        #TODO: Change to numpy function 
#        for x in xrange(master.shape[0]):
#            for y in xrange(master.shape[1]):
#                missing_count=0
#                for i in xrange(len(self.decadals)):
#                    
#                    if master[x,y,i] > 0.9e20:
#                        missing_count+=1
#                if missing_count > round(0.1*len(self.decadals)):
#                    missing_value_mask[x,y] = 1e20
#                else:
#                    missing_value_mask[x,y] = 1
#                #print missing_count
#        misval_new = FileHandler.saveToNetCDF(missing_value_mask, observations.values()[0], 'missing_value_new')


        for year in self.decadals:
            
            tmp_array = FileHandler.openNetCDFFile(observations[year], 'var')
            #tmp_array = np.expand_dims(tmp_array, axis=2)
            tmp_array_1d = tmp_array.reshape((1,tmp_array.size))
            try:
                master = np.concatenate((master,tmp_array_1d))
            except NameError:
                orig_shape = tmp_array.shape
                master = tmp_array_1d  
        #sum along axis
        sum_array = np.sum(master,axis=0)
        #missing_value_mask = np.zeros((1,tmp_array.size))
        missing_value_mask = np.where(sum_array<0.1*(len(self.decadals)*1e20),sum_array,1e20)
        missing_value_mask = np.where(sum_array>=0.1*(len(self.decadals)*1e20),missing_value_mask,1)
                #print missing_count
        missing_value_mask = np.reshape(missing_value_mask, orig_shape)
        misval_new = FileHandler.saveToNetCDF(missing_value_mask, observations.values()[0], 'missing_value_new2')
        misval_new = cdo.timmean(input=misval_new, output=misval_new+'timmean')
#        tmpList = list()
#        for year in self.decadals:    
#            #Mark every gridpoint as missing, if it is missing once
#            tmpList.append(cdo.timavg(input=observations[year], output=self.tmpDir+'misvalmask_tmp'+str(year)))  
#            #Mark only as missing if it is missing the whole year
#            print observations[year]
#            #tmpList.append(cdo.yearmean(input=observations[year], output=self.tmpDir+'misvalmask_tmp'+str(year))) 
#        misval_tmp = cdo.mergetime(input=' '.join(tmpList), output=self.tmpDir+'misvalmask_tmp_all') 
#        misval = cdo.timavg(input=misval_tmp, output=self.tmpDir+'misvalmask_all_years')
#        misval_old = cdo.setrtoc('-1000,1000,1', input=misval, output=self.tmpDir+'missingValueMask.nc')                    

        return misval_new        
                
    
    def sellonlatbox(self, fileList):
        '''
        Method to select a lon-lat-box using CDO
        
        :param fileList: list with files of single file 
        :return: new filelist 
        '''
        if type(fileList) == list:        
            #Mult Processor            
            count = len(fileList)
            poolArgs = self.getPoolArgs(count,self,fileList,'_sellonlatbox')
            resultList = self.multiProcess(poolArgs)     
            return resultList 
        else:
            return self._sellonlatbox(fileList)
    
    def _sellonlatbox(self, file):
        '''
        Single process for selecting lon-lat-box with cdo
        
        :param: file
        :return: new file
        '''
        return cdo.sellonlatbox(self.lonlatbox, input=file, output=self.tmpDir+self.extractFilename(file)+'_sellonlatbox')
    
    def fieldMean(self, fileList):
        '''
        Multiprocess Method to calc fieldmean
        
        :param fileList:
        :return: fieldmean list
        '''
        #Mult Processor            
        count = len(fileList)
        poolArgs = self.getPoolArgs(count,self,fileList,'_fieldMean')
        resultList = self.multiProcess(poolArgs) 
        return resultList
    
    def _fieldMean(self, fn):
        '''
        Single process calculating fieldmean using CDO
        
        :param fn:
        :return: fieldmean fn
        '''
        return cdo.fldmean(input=fn, output=self.tmpDir+self.extractFilename(fn)+'_fldmean')
    
    def multiProcessCdo(self, fileList, cdo_command, arguments=''):
        '''
        Multi Process Wrapper for "single" cdo commands like "fldmeam" with 1 input and 1 output file
        '''
        count = len(fileList)
        poolArgs = self.getPoolArgs(count,self,fileList,cdo_command, arguments, '_singleProcessCdo')
        return self.multiProcess(poolArgs)
        
    def _singleProcessCdo(self, fn, cdo_command, arguments=''):
        '''
        Single Process Wrapper for "single" cdo commands like "fldmeam" with 1 input and 1 output file
        '''
        cdo_command_method = getattr(cdo, cdo_command) 
        return cdo_command_method(arguments, input=fn, output=self.tmpDir+self.extractFilename(fn)+'_'+cdo_command)
        
    def _plotField(self, fileName, vmin, vmax):
        '''
        @deprecated: use the static Plotter class instead
        Plot any field variable
        
        :param fileName: filepath
        :param vmin: min value for colorbar
        :param vmax: max value for colorbar
        ''' 
        Plotter.plotField(fileName, vmin, vmax, colormap='', output_folder=self.outputPlots, lonlatbox=self.lonlatbox)
        Plotter.saveFig(self.outputPlots, fileName.split(self.outputDir)[-1])
       
    def applyMissingValueMask(self, start, end, outputDirList):
        '''
        Multiply results with a missing value mask
        
        :todo: At the moment all files in the outputdir are multiplied. This should be changed
        :param start: startyear
        :param end: endyear
        :return: list of changed files
        '''
        fileList = list()
        for outputDir in outputDirList:
            for files in os.listdir(outputDir):
                if(os.path.isfile(outputDir+files) and (outputDir+files).split('_')[-1] != 'masked'):
                    file_parts = str(files).split('_')
                    if(file_parts[0] == start and file_parts[1] == end):
                        if self.maskMissingValues:
                            fileList.append(cdo.mul(input = ' '.join([outputDir+files, self.misvalMask]),
                                                    output = outputDir+files+'_masked'))
                        else:
                            fileList.append(outputDir+files) 
        return fileList 
    
    
    def detrendTimeSeries(self, fn, keepMean=True):
        '''
        Subtracts trend of a timeseries. 
        If keepMean is True the mean is kept 
        '''
#        print fn
#        #split file in mothly files
#        tmp = cdo.copy(input=fn,output=fn+'copy',options='-f nc')
#        tmp = cdo.splitmon(input=tmp,output=fn+'monthly',options='-f nc')
#        
#        months = ['01','02','03','04','05','06','07','08','09','10','11','12']
#        
#        detrended_list = list()
#        for mon in months:
#            fn_mon = tmp+mon+'.nc'
#            trend_a = fn_mon + 'trend_a'
#            trend_b = fn_mon + 'trend_b'
#            cdo.trend(input=fn_mon, output=' '.join([trend_a,trend_b]))
#            
#            if keepMean:
#                trend_a = cdo.mulc('0.00',input=trend_a,output=trend_a+'_keepMean')
#                
#            detrended_tmp = cdo.subtrend(input=' '.join([fn_mon,trend_a,trend_b]), output=fn_mon+'_detrended')
#            detrended_list.append(detrended_tmp)
#
#        detrended = cdo.mergetime(' '.join(detrended_list), output=fn+'_detrended')    
        detrended = cdo.detrend(input=fn,output=fn+'_detrended')
        return detrended
        
    def getMissingMaskForFieldMean(self):
        '''
        If we calculate field mean we have to take the missing values of the observations into account.
        In this method the missing field mask of all observations is calculated. 
        DON'T mix it up with calcMissingValueMask
        '''
        yearly_obs = dict()
        for year in self.decadals:
            yearly_obs[year] = cdo.yearmean(input=self.observationRemapped[year], output=self.tmpDir+self.extractFilename(self.observationRemapped[year]+'_ymean'))
        mismask = cdo.mergetime(input=' '.join(yearly_obs.values()), output=self.tmpDir+'mergedobs')
        mismask = cdo.timavg(input=mismask, output=mismask+'_tim_avg')
        mismask = cdo.setrtoc('-1e19,1e19,1', input=mismask, output=mismask+'_setrtoc'  )
        return mismask
    
    def applyMissingMaskForFieldMean(self, fileList, missmask):
        '''
        Apply missing value mask to a list of files
        '''
        count = len(fileList)
        poolArgs = self.getPoolArgs(count,self,fileList,missmask,'_applyMissingMaskForFieldMean')
        resultList = self.multiProcess(poolArgs) 
        return resultList
    
    def _applyMissingMaskForFieldMean(self, data, missmask):
        '''
        Apply missing value mask to a file
        '''
        return cdo.mul(input=' '.join([data,missmask]), output=data+'_masked' )
        
    def getLevelIntersection(self,fn1,fn2,fn3=None):
        '''
        Takes 2 3d netcdf files and calculates common levels
        
        :return: list of common levels
        '''
        def roundList(il):
            t = il[0].split(' ')
            t = map(float,t)
            t = map(round,t)
            t = map(int,t)
            return map(str,t)
        level_f1 = roundList(cdo.showlevel(input=fn1))
        level_f2 = roundList(cdo.showlevel(input=fn2))
        result = [val for val in level_f1 if val in level_f2]
        
        if fn3 is not None:
            level_f3 = roundList(cdo.showlevel(input=fn3))
            result = [val for val in result if val in level_f3]             
        return result
        
        
        
        
        
        
        
        
           
        
         
        
        
    
        