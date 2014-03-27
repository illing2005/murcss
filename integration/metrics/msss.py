'''
Created on 19.02.2013

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

from cdo import *
cdo = Cdo()
from time import time

from metricAbstract import MetricAbstract
from taylorplot import TaylorPlotMurCSS

class MsssError(Exception): pass
    

    
class Msss(MetricAbstract):
    '''
    Class to calculate the MSSS. Call method "analyze" to get results
    '''
    #Here the output filenames and folder structure get defined
    outputFolderStruct = 'TIME_variable_project1_product1_institute1_model1_experiment1_project2_product2_institute2_model2_experiment2_obsExp_LEVEL_SELLONLAT_YEARRANGE'.split('_')
    insideFolderNameFlag = 'project_product_model_experiment'.split('_')
    fileNameFlag = 'LEADTIMES_variable_project_product_model_experiment_YEARRANGE'.split('_')
    fileNameFlagVs = 'LEADTIMES_variable_project1_product1_model1_experiment1_vs_project2_product2_model2_experiment2_YEARRANGE'.split('_')

    def __init__(    self, 
                     output='/tmp/msss/output/', 
                     output_plots='/tmp/msss/output/', 
                     decadals='1960,1965,1970,1975,1980,1985,1990,1995,2000',
                     variable='tas', 
                     
                     #CMOR parameter for input1
                     project1='baseline1', 
                     product1='output',
                     institute1='mpi-m', 
                     model1 = 'mpi-esm-lr', 
                     experiment1='decs4e', 
                     ensemblemembers1='*', 
                     
                     #CMOR parameter for input2
                     project2='baseline0', 
                     product2='output1', 
                     institute2='mpi-m', 
                     model2='mpi-esm-lr', 
                     experiment2='decadal',
                     ensemblemembers2='*', 
                     
                     leadtimes='1,2-9', 
                     observation = 'HadCrut', 
                     maskMissingValues = True,
                     result_grid=None,
                     cache='/tmp/msss/cache/', 
                     baseDir = '..', 
                     timeFreq = 'mon',

                     #for bootstrapping
                     bootstrap=None, 
                     bootstrap2=None, 
                     obsRemapped = None, 
                     
                     analysis_type = 'decadal',
                     level=None, 
                     lonlatbox=None, 
                     fieldmean=False, #not finally implemented
                     
                     colormap='RedBlu', #not used
                     observation_type = 'REANALYSIS', #not used ???
                     reanpath='', #not used ???
                     months = None,
                 ):
        '''
        Constructor. Set parameters and calls parent constructor (metricAbstract)
        
        :param output: Path where the results are saved
        :param output_plots: Path where the plots are saved
        :param decadals: a list with decadal experiments, i.e. [1960,1960,...,1995]
        :param variable: CMOR-Name of variable like 'tas' (Near Surface Temperature)
        :param project1: CMOR-Name of the project like CMIP5 or BASELINE1
        :param product1: CMOR-Parameter
        :param institute1: CMOR-Name of the institute like MPI-M
        :param model1: CMOR-Name of the used model like MPI-ESM-LR
        :param experiment1: Experiment name like decs4e or historical
        :param project2: CMOR-Name of the project like CMIP5 or BASELINE1
        :param product2: CMOR-Parameter
        :param institute2:CMOR-Name of the institute like MPI-M
        :param model2: CMOR-Name of the used model like MPI-ESM-LR
        :param experiment2: Experiment name like decs4e or historical
        :param ensemblemembers: Comma separated string of esemble members like 'r1i1p1,r2i1p1,...'
        :param leadtimes: Leadtimes to analyze, ie. 1,2-9
        :param observation: Observation or Reanalysis Experiment, HadCrut or ERA-Int
        :param maskMissingValues: Boolean 
        :param result_grid: Griddescription of resultgrid like r72x36
        :param cache: Path for cachedir during analysis
        :param baseDir: Path of Class to find default files
        :param timeFreq: Timefrequency of the files, ie 'mon' or 'day'
        :param bootstrap:
        :param bootstrap2:
        :param obsRemapped:
        :param level: For 3D-Files --> Select a single level
        :param lonlatbox: If you want to select a specific lonlat-Box
        :param fieldmean: Boolean - If you want to calculate global mean
        :param colormap: Matplotlib Colormap used for ploting
        :param observation_type: NOT USED ANYMORE!!!
        :param reanpath: NOT USED ANYMORE!!!
        '''            
        #Add kwargs to attributes
        listOfAttributes = ['decadals','variable',
                            'project1','product1','institute1','model1','experiment2',
                            'project2','product2','institute2','model2','experiment1',
                            'leadtimes','maskMissingValues','level','lonlatbox','fieldmean',
                            'bootstrap','bootstrap2','obsRemapped',
                            'observation_type','colormap','timeFreq']
        for attr in listOfAttributes:
            setattr(self, attr, locals()[attr])
               
        
        if ensemblemembers1 != None and ensemblemembers1 != '*':
            self.ensemblemembers1 = ensemblemembers1.split(',')
        else:
            self.ensemblemembers1 = '*'
        if ensemblemembers2 != None and ensemblemembers2 != '*':
            self.ensemblemembers2 = ensemblemembers2.split(',')
        else:
            self.ensemblemembers2 = '*'

        super(Msss,self).__init__(tmpDir = self.checkPath(cache), output=output, output_plots=output_plots,
                                  baseDir=baseDir, result_grid=result_grid, observation=observation,
                                  analysis_type=analysis_type, decadals=decadals, months=months )


    def prepareInput(self):
        self.findFiles.observation = self.observation  #WORKAROUND
             
        self.observationDict,self.input1Remapped,self.observationRemapped,self.input2Remapped = [{} for dummy in range(4)]        
        if self.bootstrap is None:
        
            print 'Searching Files' 
            countYears = len(self.decadals)
            poolArgs = self.getPoolArgs(countYears, self.findFiles, list(self.decadals), self.project2, self.model2, self.variable,
                                        'mon', self.product2, self.ensemblemembers2, self.institute2, self.experiment2,
                                        self.maxLeadtime, 'getFiles')
            self.input2Dict = self.listToYeardict(self.multiProcess(poolArgs)) 
            #print self.input2Dict
            poolArgs = self.getPoolArgs(countYears, self.findFiles, list(self.decadals), self.project1,self.model1,self.variable,
                                        'mon',self.product1,self.ensemblemembers1,self.institute1,self.experiment1,
                                        self.maxLeadtime,'getFiles')
            self.input1Dict = self.listToYeardict(self.multiProcess(poolArgs))
            for year in self.decadals:
                self.observationDict[year] = self.findFiles.getReanalysis(year, self.observation_type, self.obsExp, self.variable, maxLeadtime=self.maxLeadtime)     
            print 'Remapping Files'
            for year in self.decadals:
                self.input1Remapped[year] = self._remapFiles(self.input1Dict[year], flag=self.product1)
                self.input2Remapped[year] = self._remapFiles(self.input2Dict[year], flag=self.product2)
                self.observationRemapped[year] = self.remapFile(self.observationDict[year])
            self.ensList = self.input1Dict
            self.histList = self.input2Dict   
        else:
            for year in self.decadals:
                self.input2Remapped[year] = self.bootstrap2[year]
                self.input1Remapped[year] = self.bootstrap[year]
                if self.obsRemapped == None:
                    self.observationDict[year] = self.findFiles.getReanalysis(year, self.observation_type, self.obsExp, self.variable)
                    self.observationRemapped[year] = self.remapFile(self.observationDict[year])
                else:
                    self.observationRemapped[year] = self.obsRemapped[year]

    def analyze(self):
        '''
        Main method to calculate the MSSS. Following steps are performed:
        
        1. Searching for file using solr_search
        2. Remapping Files to coarser grid
        3. Calculating Ensemble Mean
        4. Calculating cross-validated mean
        5. Calculating anomalies for both models
        6. Analyzing different timeranges 
        '''  
        
        #check if prepare output was called or variables are set
        try:
            self.input1Remapped
            self.input2Remapped
            self.observationRemapped
        except AttributeError:
            raise MsssError, 'Please run "prepareInput() first or provide input data'
            
        #select lon/lat box
        if self.lonlatbox is not None:
            print 'Selecting lon-lat-box %s' %(self.lonlatbox)
            for year in self.decadals:
                self.input1Remapped[year] = self.sellonlatbox(self.input1Remapped[year])
                self.input2Remapped[year] = self.sellonlatbox(self.input2Remapped[year])
                self.observationRemapped[year] = self._sellonlatbox(self.observationRemapped[year])         

        if self.fieldmean:
            print 'Calculating field mean'
            for year in self.decadals:
                self.input1Remapped[year] = self.fieldMean(self.input1Remapped[year])
                self.input2Remapped[year] = self.fieldMean(self.input2Remapped[year])
                self.observationRemapped[year] = self._fieldmean(self.observationRemapped[year]) 
            

        if self.bootstrap is None:
            print 'Calculating ensemble mean'
            input1Ensemblemean = self.getEnsembleMean(self.input1Remapped)
            input2Ensemblemean = self.getEnsembleMean(self.input2Remapped)
        else:
            input1Ensemblemean = self.input1Remapped
            input2Ensemblemean = self.input2Remapped
            
        print 'Calculating crossvalidated mean'
        input1CrossvalMean = self.getCrossValMean(input1Ensemblemean, 'model')
        input2CrossvalMean = self.getCrossValMean(input2Ensemblemean, 'hist')
        if self.obsExp != 'HadCrut':
            if self.obsRemapped == None:
                obsCrossValMean = self.getCrossValMean(self.observationRemapped, 'obs')
                self.observationRemapped = self.getAnomalies(self.observationRemapped, obsCrossValMean)
        self.obsRemapped = self.observationRemapped
            
        print 'Calculating Anomalies'
        input1Anomalies = self.getAnomalies(input1Ensemblemean, input1CrossvalMean)
        input2Anomalies = self.getAnomalies(input2Ensemblemean, input2CrossvalMean)
        self.constantField = self.createConstantFile(self.gridFile)
        
        self.misvalMask = self.calcMissingValueMask(self.observationRemapped)
        
        rangeCount = len(self.getYearRange())
        firstYearList = list()
        lastYearList = list()
        
        for yr in self.getYearRange():
            firstYearList.append(yr[0])
            lastYearList.append(yr[1])

        #escape attributes 
        self.escapeAttributes()

        poolArgs = self.getPoolArgs(rangeCount, self, input1Anomalies,self.observationRemapped,input2Anomalies,
                                    firstYearList,lastYearList,'analyzeYearRange')
        resultList = self.multiProcess(poolArgs)
        
        if self.fieldmean:
            #construct a taylor plot
            taylor = TaylorPlotMurCSS(negativeCorr=True)
            taylor.constructPlot(resultList, self.outputPlots)
        else:
            #plot fields
            if self.bootstrap is None:
                fileList = list()
                for rList in resultList:  
                    fileList = fileList + rList          
                filesToPlot = len(fileList)
                #print fileList
                min = np.array([-1]*filesToPlot)
                max = np.array([1]*filesToPlot)
                std_ind = [i for i,s in enumerate(fileList) if 'std_ratio' in s]
                min[std_ind] = 0.5
                max[std_ind] = 2
                poolArgs = self.getPoolArgs(filesToPlot,self,fileList,list(min),list(max),'_plotField')
                resultList = self.multiProcess(poolArgs)
                  
        
           
    def analyzeYearRange(self, input1Anomalies, observationRemapped, input2Anomalies, startYear, endYear):
        '''
        Calculate the MSSS for given timerange. The following steps are perfomed in this method:
        
        1. temporal smoothing
        2. calculation of correlation, bias, and msss
        3. calculation of model1 vs. model2
        4. plotting of the 9 fields 
        
        :todo: maybe this could be done in multiprocessing --> so different timeranges at once
        
        :param input1Anomalies: dict with anomalies of model1
        :param observationRemapped: dict with anomalies of observations
        :param input2Anomalies: dict with anomalies of model2
        :param startYear: 
        :param endYear: 
        
        '''
        print 'Analyzing year %s to %s' %(startYear, endYear)

        base = '%s-%s'%(startYear, endYear)
        #get file names
        flag1 = self.constructName(self.insideFolderNameFlag,exp='1')
        flag2 = self.constructName(self.insideFolderNameFlag,exp='2')
        fnFlag1 = self.constructName(self.fileNameFlag, exp='1', startYear=str(startYear), endYear=str(endYear))
        fnFlag2 = self.constructName(self.fileNameFlag, exp='2', startYear=str(startYear), endYear=str(endYear))
        fnFlagVs = self.constructName(self.fileNameFlagVs, exp='', startYear=str(startYear), endYear=str(endYear))

        #make outputfolders
        if self.bootstrap is None:
            model1Output = self.makeFolder(self.outputDir+base+'/'+flag1+'/msss')
            model2Output = self.makeFolder(self.outputDir+base+'/'+flag2+'/msss')
            vsOutput = self.makeFolder(self.outputDir+base+'/%s_vs_%s/msss'%(flag1,flag2))
            #make plot folders
            model1Plot = self.makeFolder(self.outputPlots+base+'/'+flag1+'/msss')
            model2Plot = self.makeFolder(self.outputPlots+base+'/'+flag2+'/msss')
            vsPlot = self.makeFolder(self.outputPlots+base+'/%s_vs_%s/msss'%(flag1,flag2))
        else:
            model1Output = self.outputDir
            model2Output = self.outputDir
            vsOutput = self.outputDir
                 
        yearCount = len(self.decadals)
        poolArgs = self.getPoolArgs(yearCount,self,self.yeardictToList(input1Anomalies),startYear,endYear,'_tempSmoothing')
        ensListSelYearTimmean = self.listToYeardict(self.multiProcess(poolArgs))
        poolArgs = self.getPoolArgs(yearCount,self,self.yeardictToList(observationRemapped),startYear,endYear,'_tempSmoothing')
        obsSelYearTimmean = self.listToYeardict(self.multiProcess(poolArgs))
        poolArgs = self.getPoolArgs(yearCount,self,self.yeardictToList(input2Anomalies),startYear,endYear,'_tempSmoothing')
        histSelYearTimmean = self.listToYeardict(self.multiProcess(poolArgs))    
    
        startYear = str(startYear)
        endYear = str(endYear)
        correlation = self.getCorrelationAndVariancesAndMSSS(ensListSelYearTimmean,obsSelYearTimmean, fnFlag1, model1Output)
        correlationHist = self.getCorrelationAndVariancesAndMSSS(histSelYearTimmean,obsSelYearTimmean, fnFlag2, model2Output)
        
        msssInitVsUninit = self.getInitVsUninit(correlation, correlationHist, fnFlagVs, vsOutput)
        
        rmss = self.getRMSS(correlation['msss'], fnFlag1, model1Output)
        rmssHist = self.getRMSS(correlationHist['msss'], fnFlag2, model2Output)
        rmssInitvsUninit = self.getRMSS(msssInitVsUninit, fnFlagVs, vsOutput)  
        
        return self.applyMissingValueMask(startYear, endYear, [model1Output,model2Output,vsOutput])

    def getInitVsUninit(self, initResults, uninitResults, flag, outputDir):
        '''
        Calculates comparing values between two models/modelversions i.e. initialized vs. uninitialized.
        Following variables are computed:
        
        - correlation: corr1 - corr2 
        - conditional bias: abs(bias1) - abs(bias2)
        - msss: (msss1-msss2)/(1-msss2)
        
        :param initResults: dict with "corr", "bias", and "msss" of first model
        :param uninitResults: dict with "corr", "bias", and "msss" of second model
        :param tempRange: start- and endyear --> used for filenames 
        :return: msss filename
        '''
        constantField = self.createConstantFile(self.gridFile, flag=flag)    #Create Constant Field for MSSS
        
        corr = cdo.sub(input=' '.join([initResults['corr'], uninitResults['corr']]), 
                       output = outputDir+flag+'_correlation.nc')
        bias = cdo.sub(input='-abs ' + uninitResults['bias'] + ' -abs ' + initResults['bias'], 
                       output=outputDir+flag+'_conditional_bias.nc')
        #MSSS vs Plot        
        msss = cdo.div(input='-sub ' + initResults['msss'] + ' ' + uninitResults['msss'] + 
                             ' -sub ' + constantField + ' ' + uninitResults['msss'], 
                       output=outputDir+flag+'_msss.nc')
        
        #conditional bias slope 
        slopeBiasInit = cdo.abs(input=initResults['biasSlopeMinusOne'], 
                                                     output=self.tmpDir+self.extractFilename(initResults['biasSlopeMinusOne'])+'_abs')
        slopeBiasUninit = cdo.abs(input=uninitResults['biasSlopeMinusOne'], 
                                                     output=self.tmpDir+self.extractFilename(uninitResults['biasSlopeMinusOne'])+'_abs')
        slopeBias = cdo.sub(input=' '.join([slopeBiasInit, slopeBiasUninit]),
                             output=outputDir+flag+'_biasslope-1.nc')
        
        #conditional bias skill score
        biasSKillScore = cdo.sub(input=' '.join([self.constantField,
                                                 cdo.div(input=' '.join([slopeBiasInit, slopeBiasUninit]), 
                                                         output=self.tmpDir+flag+'_biasskillDIV.nc')]),
                                 output=outputDir+flag+'_biasSkillScore.nc')
        
        return msss
           
    def getCorrelationAndVariancesAndMSSS(self, hindcast, observation, flag, outputDir):
        '''
        Calculates the following quantities:
        
        - correlation bewtween observation and hindcast
        - temporal variance of hindcast
        - temporal variance of observation
        - conditional bias 
        - msss
        - biasSlope so/sh * c
        - std ratio so/sh
        
        :param hindcast: dict of hindcasts
        :param observation: dict of observations
        :param flag: additional string for filenames
        :return: dict with "msss", "bias", "corr", "biasSlope, "stdRatio", "biasSlopeMinusOne" 
        '''
        hindcastStr = ''
        observationStr = ''
        for year in self.decadals:
            hindcastStr += ' ' + hindcast[year]
            observationStr += ' ' + observation[year]
        
        mergeHindcast = cdo.mergetime(input=hindcastStr, output=hindcast[self.decadals[0]]+flag+self.getRandomStr()+'_merged.nc')
        mergeObservations = cdo.mergetime(input=observationStr, output=observation[self.decadals[0]]+flag+self.getRandomStr()+'_merged.nc')
        varianceHindcast = cdo.timstd(input=mergeHindcast, output=hindcast[self.decadals[0]]+flag+'_variance')
        varianceObservation = cdo.timstd(input=mergeObservations, output=observation[self.decadals[0]]+flag+'_variance')
        corr = cdo.timcor(input=' '.join([mergeHindcast, mergeObservations]), output=outputDir+flag+'_correlation.nc')
        
        #standard deviation ratio
        stdRatio = cdo.div(input=' '.join([varianceObservation, varianceHindcast]), 
                           output=outputDir+flag+'_std_ratio.nc')
        
        #slope of conditional bias
        biasSlope = cdo.mul(input=' '.join([stdRatio, corr]),
                            output=outputDir+flag+'_biasslope.nc')
        
        #slope -1
        biasSlopeMinusOne = cdo.sub(input=' '.join([biasSlope, self.constantField]),
                                    output=outputDir+flag+'_biasslope-1.nc')
        bias = cdo.sub(input=corr+' '+ '-div ' + varianceHindcast + ' ' + varianceObservation, 
                       output=outputDir+flag+'_conditional_bias.nc')
        msss = cdo.sub(input='-sqr ' + corr + ' -sqr ' + bias, 
                       output=outputDir+flag+'_msss.nc')
        
        return dict(msss= msss, bias= bias, corr= corr, biasSlope=biasSlope, stdRatio=stdRatio, biasSlopeMinusOne=biasSlopeMinusOne)
        
    def getMSE(self, anomalies, obsAnom):
        '''
        Method calculates the Mean Squared Error. 
        
        :note: method ist not used at the moment
        :param anomalies: dictionary with anomalies
        :param obsAnom: dict with observation anomalies
        :return: file with mse field
        '''        
        MSEList = list()
        for year in self.decadals:
            subTmp = cdo.sub(input=' '.join([anomalies[year], obsAnom[year]]))
            MSEList.append(cdo.sqr(input=subTmp)) 
        return cdo.ensmean(input=' '.join(MSEList), output=self.tmpDir+'mse.nc')
    
    def getRMSS(self, MSSS, flag, outputDir):
        '''
        Calculates the RMSSS using the MSSS
        
        :param MSSS: netcdf file
        :param flag: name flag
        :param outputDir: output path 
        :return: file with RMSSS field 
        '''      
        tmp = cdo.sqrt(input="-sub " + self.constantField + " " + MSSS,                 
                       output=self.tmpDir+flag+'rmss_sqrt.nc')
        return cdo.sub(input=' '.join([self.constantField,tmp]), output=outputDir+flag+'_rmss.nc')
