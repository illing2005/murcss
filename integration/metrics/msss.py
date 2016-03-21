"""
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
"""

from cdo import *


from metricAbstract import MetricAbstract
from plotter import Plotter
import murcss_config

cdo = Cdo()


class MsssError(Exception):
    pass

    
class Msss(MetricAbstract):
    """
    Class to calculate the MSSS. Call method "analyze" to get results
    
    :param output: Path where the results are saved
    :param output_plots: Path where the plots are saved
    :param decadals: a list with decadal experiments, i.e. [1960,1960,...,1995]
    :param variable: CMOR-Name of variable like 'tas' (Near Surface Temperature)
    :param project1: CMOR-Name of the project like CMIP5 or BASELINE1
    :param product1: CMOR-Parameter
    :param institute1: CMOR-Name of the institute like MPI-M
    :param model1: CMOR-Name of the used model like MPI-ESM-LR
    :param experiment1: Experiment name like decadal or historical
    :param ensemblemembers1: Comma separated string of esemble members like 'r1i1p1,r2i1p1,...'
    :param project2: CMOR-Name of the project like CMIP5 or BASELINE1
    :param product2: CMOR-Parameter
    :param institute2: CMOR-Name of the institute like MPI-M
    :param model2: CMOR-Name of the used model like MPI-ESM-LR
    :param experiment2: Experiment name like decs4e or historical
    :param ensemblemembers2: Comma separated string of esemble members like 'r1i1p1,r2i1p1,...'
    :param leadtimes: Leadtimes to analyze, ie. 1,2-9
    :param observation: Observation or Reanalysis Experiment, HadCrut or ERA-Int
    :param observation_ensemble: If your observation has multiple ensemble members. Not necesary if you use an observation file    
    :param maskMissingValues: Boolean
    :param result_grid: Griddescription of resultgrid like r72x36
    :param cache: Path for cachedir during analysis
    :param baseDir: Path of Class to find default files
    :param timeFreq: Timefrequency of the files, ie 'mon' or 'day'
    :param basic_output: Boolean - If true only basic out is produced (Correlation, MSESS, and Conditional Bias). 
    :param bootstrap:
    :param bootstrap2:
    :param obsRemapped:
    :param level: For 3D-Files --> Select a single level
    :param lonlatbox: If you want to select a specific lonlat-Box
    :param fieldmean: Boolean - If you want to calculate field means
    :param zonalmean: Boolean - If you want to calculate zonal means
    :param colormap: Matplotlib Colormap used for ploting
    :param observation_type: NOT USED ANYMORE!!!
    :param reanpath: NOT USED ANYMORE!!!
    :param months: Usefull if your analysing monthly and your experiment name has a month included
    """
    # Here the output file names and folder structure get defined
    outputFolderStruct = 'TIME_variable_project1_product1_institute1_model1_experiment1_project2_product2_institute2_model2_experiment2_obsExp_LEVEL_SELLONLAT_YEARRANGE'.split('_')
    insideFolderNameFlag = 'project_product_model_experiment'.split('_')
    fileNameFlag = 'LEADTIMES_variable_project_product_model_experiment_YEARRANGE'.split('_')
    fileNameFlagVs = 'LEADTIMES_variable_project1_product1_model1_experiment1_vs_project2_product2_model2_experiment2_YEARRANGE'.split('_')

    def __init__(self,
                 output='/tmp/msss/output/',
                 output_plots='/tmp/msss/output/',
                 decadals='1960,1965,1970,1975,1980,1985,1990,1995,2000',
                 variable='tas',
                     
                 # CMOR parameter for input1
                 project1='baseline1',
                 product1='output',
                 institute1='mpi-m',
                 model1='mpi-esm-lr',
                 experiment1='decs4e',
                 ensemblemembers1='*',

                 # CMOR parameter for input2
                 project2='baseline0',
                 product2='output1',
                 institute2='mpi-m',
                 model2='mpi-esm-lr',
                 experiment2='decadal',
                 ensemblemembers2='*',
                     
                 leadtimes='1,2-9',
                 observation='HadCrut',
                 observation_ensemble='*',
                 maskMissingValues=True,
                 result_grid=None,
                 cache='/tmp/msss/cache/',
                 baseDir='..',
                 timeFreq='mon',

                 basic_output = True,
                 # for bootstrapping
                 bootstrap=None,
                 bootstrap2=None,
                 obsRemapped=None,
                     
                 leadtimes_mode='yearly',
                 level=None,
                 lonlatbox=None,
                 fieldmean=False,
                 zonalmean=False,
                     
                 colormap='RedBlu', #not used
                 observation_type = 'REANALYSIS', #not used ???
                 reanpath='', #not used ???
                 months=None):
        """
        Constructor. Set parameters and calls parent constructor (metricAbstract)
        """
        # Add kwargs to attributes
        listOfAttributes = ['decadals','variable',
                            'project1', 'product1', 'institute1', 'model1', 'experiment2',
                            'project2', 'product2', 'institute2', 'model2', 'experiment1',
                            'leadtimes', 'maskMissingValues', 'level', 'lonlatbox', 'fieldmean', 'zonalmean',
                            'bootstrap', 'bootstrap2', 'obsRemapped', 'observation_ensemble',
                            'observation_type', 'colormap', 'timeFreq', 'basic_output']
        for attr in listOfAttributes:
            setattr(self, attr, locals()[attr])
               
        if ensemblemembers1 is not None and ensemblemembers1 != '*':
            self.ensemblemembers1 = ensemblemembers1.split(',')
        else:
            self.ensemblemembers1 = '*'
        if ensemblemembers2 is not None and ensemblemembers2 != '*':
            self.ensemblemembers2 = ensemblemembers2.split(',')
        else:
            self.ensemblemembers2 = '*'

        super(Msss, self).__init__(tmpDir=self.checkPath(cache), output=output, output_plots=output_plots,
                                   baseDir=baseDir, result_grid=result_grid, observation=observation,
                                   leadtimes_mode=leadtimes_mode, decadals=decadals, months=months)

    def prepareInput(self):
        self.findFiles.observation = self.observation  # WORKAROUND
             
        self.observationDict, self.input1Remapped, self.observationRemapped, self.input2Remapped = [{} for dummy in range(4)]
        if self.bootstrap is None:
        
            print 'Searching Files' 
            countYears = len(self.decadals)

            poolArgs = self.getPoolArgs(countYears, self.findFiles, list(self.decadals), self.project1, self.model1,
                                        self.variable, 'mon', self.product1, self.ensemblemembers1, self.institute1,
                                        self.experiment1, self.maxLeadtime, self.minLeadtime, 'getFiles')
            self.input1Dict = self.listToYeardict(self.multiProcess(poolArgs))
            poolArgs = self.getPoolArgs(countYears, self.findFiles, list(self.decadals), self.project2, self.model2,
                                        self.variable, 'mon', self.product2, self.ensemblemembers2, self.institute2,
                                        self.experiment2, self.maxLeadtime, self.minLeadtime, 'getFiles')
            self.input2Dict = self.listToYeardict(self.multiProcess(poolArgs)) 
            for year in self.decadals:
                self.observationDict[year] = self.findFiles.getReanalysis(year, self.observation_type, self.obsExp,
                                                                          self.variable, maxLeadtime=self.maxLeadtime,
                                                                          observation_ensemble=self.observation_ensemble,
                                                                          minLeadtime=self.minLeadtime)
            print 'Remapping Files'
            obslist_tmp = list()
            for year in self.decadals:
                self.input1Remapped[year] = self._remapFiles(self.input1Dict[year], flag=self.product1)
                self.input2Remapped[year] = self._remapFiles(self.input2Dict[year], flag=self.product2)
                obslist_tmp.append(self.observationDict[year])
            self.observationRemapped = self.listToYeardict(self._remapFiles(obslist_tmp))
            
            self.ensList = self.input1Dict
            self.histList = self.input2Dict   
        else:
            for year in self.decadals:
                self.input2Remapped[year] = self.bootstrap2[year]
                self.input1Remapped[year] = self.bootstrap[year]
                if self.obsRemapped is None:
                    self.observationDict[year] = self.findFiles.getReanalysis(year, self.observation_type, self.obsExp,
                                                                              self.variable, maxLeadtime=self.maxLeadtime,
                                                                              observation_ensemble=self.observation_ensemble,
                                                                              minLeadtime=self.minLeadtime)
                    self.observationRemapped[year] = self.remapFile(self.observationDict[year])
                else:
                    self.observationRemapped[year] = self.obsRemapped[year]

        # select lon/lat box
        if self.lonlatbox is not None:
            print 'Selecting lon-lat-box %s' % self.lonlatbox
            for year in self.decadals:
                self.input1Remapped[year] = self.sellonlatbox(self.input1Remapped[year])
                self.input2Remapped[year] = self.sellonlatbox(self.input2Remapped[year])
                self.observationRemapped[year] = self._sellonlatbox(self.observationRemapped[year])         

        if self.fieldmean:
            print 'Calculating field mean'
            # Calculate missing value fields
            self.obsmissmask = self.getMissingMaskForFieldMean()     
            for year in self.decadals:
                # Apply Missing value Mask to all fields
                self.observationRemapped[year] = self._applyMissingMaskForFieldMean(self.observationRemapped[year],
                                                                                    self.obsmissmask)
                self.input1Remapped[year] = self.applyMissingMaskForFieldMean(self.input1Remapped[year],
                                                                              self.obsmissmask)
                self.input2Remapped[year] = self.applyMissingMaskForFieldMean(self.input2Remapped[year],
                                                                              self.obsmissmask)

                self.input1Remapped[year] = self.fieldMean(self.input1Remapped[year])
                self.input2Remapped[year] = self.fieldMean(self.input2Remapped[year])
                self.observationRemapped[year] = self._fieldMean(self.observationRemapped[year]) 

        if self.zonalmean:
            print 'Calculating zonal mean'
            level_intersection = self.getLevelIntersection(self.observationRemapped[self.decadals[0]],
                                                           self.input1Remapped[self.decadals[0]][0],
                                                           self.input2Remapped[self.decadals[0]][0])
            print 'Common levels: %s' % (','.join(level_intersection))
            for year in self.decadals:
                self.observationRemapped[year] = cdo.zonmean(input=self.observationRemapped[year],
                                                             output=self.observationRemapped[year]+'_zonmean')
                self.observationRemapped[year] = cdo.sellevel(','.join(level_intersection),
                                                              input=self.observationRemapped[year],
                                                              output=self.observationRemapped[year]+'sellevel')

                self.input1Remapped[year] = self.multiProcessCdo(self.input1Remapped[year], 'zonmean')
                self.input1Remapped[year] = self.multiProcessCdo(self.input1Remapped[year], 'sellevel',
                                                                 ','.join(level_intersection))
                
                self.input2Remapped[year] = self.multiProcessCdo(self.input2Remapped[year], 'zonmean')
                self.input2Remapped[year] = self.multiProcessCdo(self.input2Remapped[year], 'sellevel',
                                                                 ','.join(level_intersection))

    def analyze(self):
        """
        Main method to calculate the MSSS. Following steps are performed:
        
        1. Searching for file using solr_search
        2. Remapping Files to coarser grid
        3. Calculating Ensemble Mean
        4. Calculating cross-validated mean
        5. Calculating anomalies for both models
        6. Analyzing different timeranges 
        """
        
        # check if prepare output was called or variables are set
        try:
            self.input1Remapped
            self.input2Remapped
            self.observationRemapped
        except AttributeError:
            raise MsssError, 'Please run "prepareInput() first or provide input data'

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

        if self.obsRemapped is None or self.bootstrap:
            obsCrossValMean = self.getCrossValMean(self.observationRemapped, 'obs')
            self.observationRemappedAnom = self.getAnomalies(self.observationRemapped, obsCrossValMean)
        self.obsRemapped = self.observationRemappedAnom
            
        print 'Calculating Anomalies'
        input1Anomalies = self.getAnomalies(input1Ensemblemean, input1CrossvalMean)
        input2Anomalies = self.getAnomalies(input2Ensemblemean, input2CrossvalMean)
        self.constantField = self.createConstantFile(self.gridFile)
        
        if self.maskMissingValues:
            if not self.fieldmean:
                self.misvalMask = self.calcMissingValueMask(self.observationRemappedAnom)
            else:
                self.misvalMask = self.constantField
            
        rangeCount = len(self.getYearRange())
        firstYearList = list()
        lastYearList = list()
        
        for yr in self.getYearRange():
            firstYearList.append(yr[0])
            lastYearList.append(yr[1])

        # escape attributes
        self.escapeAttributes()

        if murcss_config.subsub_procs:

            poolArgs = self.getPoolArgs(rangeCount, self, input1Anomalies, self.observationRemappedAnom, input2Anomalies,
                                        firstYearList, lastYearList, 'analyzeYearRange')
            resultList = self.multiProcess(poolArgs)
        else:
            resultList = list()
            for i in range(0, rangeCount):
                resultList.append(self.analyzeYearRange(input1Anomalies, self.observationRemappedAnom, input2Anomalies,
                                                        firstYearList[i], lastYearList[i]))
 
        if self.fieldmean:
            # construct a taylor plot
            if self.bootstrap is None:
                outputPlots = self.outputPlots+'accuracy/'
                self.makeFolder(outputPlots)
                #taylor = TaylorPlotMurCSS(negativeCorr=False)
                #taylor.constructPlot(resultList, outputPlots)
                print 'Plotting fldmean'
                
                flag1 = self.constructName(self.fileNameFlag, exp='1', startYear='1', endYear='1')
                flag1 = flag1[4:]
                flag2 = self.constructName(self.fileNameFlag, exp='2', startYear='1', endYear='1')
                flag2 = flag2[4:]
                fnFlagVs = self.constructName(self.fileNameFlagVs, exp='', startYear='1', endYear='1')
                fnFlagVs = fnFlagVs[4:]
                plot_list = [('correlation', 'Anomaly Correlation', [0, 1]),
                             ('msss', 'Mean Squared Error Skill Score', [-1, 1])]
                Plotter.plotLeadtimeseries(resultList, [flag1, flag2, fnFlagVs], plot_list,
                                           ['input1', 'input2', ''])
                Plotter.saveFig(outputPlots, 'accuracy_leadtimeseries_all')
                
                plot_list = [('correlation', 'Anomaly Correlation', [0, 1]),
                             ('msss', 'Mean Squared Error Skill Score', [-1, 1])]
                Plotter.plotLeadtimeseries(resultList, [flag1], plot_list)
                Plotter.saveFig(outputPlots, flag1+'accuracy_leadtimeseries_input1')
                
                plot_list = [('correlation', 'Anomaly Correlation', [0, 1]),
                             ('msss', 'Mean Squared Error Skill Score', [-1, 1])]
                Plotter.plotLeadtimeseries(resultList, [flag2], plot_list)
                Plotter.saveFig(outputPlots, flag2+'accuracy_leadtimeseries_input2')
                
                plot_list = [('correlation', 'Anomaly Correlation', [0, 1]),
                             ('msss', 'Mean Squared Error Skill Score', [-1, 1])]
                Plotter.plotLeadtimeseries(resultList, [fnFlagVs], plot_list)
                Plotter.saveFig(outputPlots, fnFlagVs+'accuracy_leadtimeseries_versus')
        
        elif self.zonalmean:
            fileList = list()
            for rList in resultList:  
                fileList = fileList + rList          
            filesToPlot = len(fileList)
            for fn in fileList:
                Plotter.plotVerticalProfile(fn, -1, 1)
                Plotter.saveFig(self.outputPlots, fn.split(self.outputDir)[-1])
        else:
            # plot fields
            if self.bootstrap is None:
                fileList = list()
                for rList in resultList:  
                    fileList = fileList + rList          
                filesToPlot = len(fileList)

                min = np.array([-1]*filesToPlot)
                max = np.array([1]*filesToPlot)
                std_ind = [i for i, s in enumerate(fileList) if 'std_ratio' in s]
                min[std_ind] = 0.5
                max[std_ind] = 2
                poolArgs = self.getPoolArgs(filesToPlot, self, fileList, list(min), list(max), '_plotField')
                resultList = self.multiProcess(poolArgs)
                
                # Plotting without multiprocessing
#                for i,f in enumerate(fileList):     
#                    if not 'biasslope' in f:
#                        if not 'std_ratio' in f:
#                            self._plotField(f, min[i], max[i])

    def analyzeYearRange(self, input1Anomalies, observationRemapped, input2Anomalies, startYear, endYear):
        """
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
        """
        if startYear == endYear:
            print 'Analyzing leadtime %s' % startYear
        else:
            print 'Analyzing leadtime %s to %s' % (startYear, endYear)

        base = '%s-%s' % (startYear, endYear)
        # get file names
        flag1 = self.constructName(self.insideFolderNameFlag, exp='1', extra='input1')
        flag2 = self.constructName(self.insideFolderNameFlag, exp='2', extra='input2')
        fnFlag1 = self.constructName(self.fileNameFlag, exp='1', startYear=str(startYear), endYear=str(endYear))
        fnFlag2 = self.constructName(self.fileNameFlag, exp='2', startYear=str(startYear), endYear=str(endYear))
        fnFlagVs = self.constructName(self.fileNameFlagVs, exp='', startYear=str(startYear), endYear=str(endYear))

        # make outputfolders
        if self.bootstrap is None:
            model1Output = self.makeFolder(self.outputDir+base+'/'+flag1+'/accuracy')
            model2Output = self.makeFolder(self.outputDir+base+'/'+flag2+'/accuracy')
            vsOutput = self.makeFolder(self.outputDir+base+'/%s_vs_%s/accuracy' % (flag1, flag2))
            # make plot folders
            model1Plot = self.makeFolder(self.outputPlots+base+'/'+flag1+'/accuracy')
            model2Plot = self.makeFolder(self.outputPlots+base+'/'+flag2+'/accuracy')
            vsPlot = self.makeFolder(self.outputPlots+base+'/%s_vs_%s/accuracy' % (flag1, flag2))
        else:
            model1Output = self.outputDir
            model2Output = self.outputDir
            vsOutput = self.outputDir
                 
        yearCount = len(self.decadals)
        poolArgs = self.getPoolArgs(yearCount, self, self.yeardictToList(input1Anomalies),
                                    startYear, endYear, '_tempSmoothing')
        ensListSelYearTimmean = self.listToYeardict(self.multiProcess(poolArgs))
        poolArgs = self.getPoolArgs(yearCount, self, self.yeardictToList(observationRemapped),
                                    startYear, endYear, '_tempSmoothing')
        obsSelYearTimmean = self.listToYeardict(self.multiProcess(poolArgs))
        poolArgs = self.getPoolArgs(yearCount, self, self.yeardictToList(input2Anomalies),
                                    startYear, endYear, '_tempSmoothing')
        histSelYearTimmean = self.listToYeardict(self.multiProcess(poolArgs))
    
        startYear = str(startYear)
        endYear = str(endYear)
        correlation = self.getCorrelationAndVariancesAndMSSS(ensListSelYearTimmean, obsSelYearTimmean,
                                                             fnFlag1, model1Output)
        correlationHist = self.getCorrelationAndVariancesAndMSSS(histSelYearTimmean, obsSelYearTimmean,
                                                                 fnFlag2, model2Output)
        
        msssInitVsUninit = self.getInitVsUninit(correlation, correlationHist, fnFlagVs, vsOutput)
        if not self.basic_output: 
            rmss = self.getRMSS(correlation['msss'], fnFlag1, model1Output)
            rmssHist = self.getRMSS(correlationHist['msss'], fnFlag2, model2Output)
            rmssInitvsUninit = self.getRMSS(msssInitVsUninit, fnFlagVs, vsOutput)  
        
        return self.applyMissingValueMask(startYear, endYear, [model1Output, model2Output, vsOutput])

    def getInitVsUninit(self, initResults, uninitResults, flag, outputDir):
        """
        Calculates comparing values between two models/modelversions i.e. initialized vs. uninitialized.
        Following variables are computed:
        
        - correlation: corr1 - corr2 
        - conditional bias: abs(bias1) - abs(bias2)
        - msss: (msss1-msss2)/(1-msss2)
        
        :param initResults: dict with "corr", "bias", and "msss" of first model
        :param uninitResults: dict with "corr", "bias", and "msss" of second model
        :param tempRange: start- and endyear --> used for filenames 
        :return: msss filename
        """
        constantField = self.createConstantFile(self.gridFile, flag=flag)   # Create Constant Field for MSSS
        
        corr = cdo.sub(input=' '.join([initResults['corr'], uninitResults['corr']]), 
                       output=outputDir+flag+'_correlation.nc')
        bias = cdo.sub(input='-abs ' + uninitResults['bias'] + ' -abs ' + initResults['bias'], 
                       output=outputDir+flag+'_conditional_bias.nc')
        # MSSS vs Plot
        msss = cdo.div(input='-sub ' + initResults['msss'] + ' ' + uninitResults['msss'] + 
                             ' -sub ' + constantField + ' ' + uninitResults['msss'], 
                       output=outputDir+flag+'_msss.nc')
        if self.basic_output:
            return msss
        else:
            # conditional bias slope
            slopeBiasInit = cdo.abs(input=initResults['biasSlopeMinusOne'], 
                                    output=self.tmpDir+self.extractFilename(initResults['biasSlopeMinusOne'])+'_abs')
            slopeBiasUninit = cdo.abs(input=uninitResults['biasSlopeMinusOne'], 
                                      output=self.tmpDir+self.extractFilename(uninitResults['biasSlopeMinusOne'])+'_abs')
            slopeBias = cdo.sub(input=' '.join([slopeBiasInit, slopeBiasUninit]),
                                output=outputDir+flag+'_biasslope-1.nc')
        
            # conditional bias skill score
            biasSKillScore = cdo.sub(input=' '.join([self.constantField,
                                                     cdo.div(input=' '.join([slopeBiasInit, slopeBiasUninit]),
                                                             output=self.tmpDir+flag+'_biasskillDIV.nc')]),
                                     output=outputDir+flag+'_biasSkillScore.nc')
            return msss
           
    def getCorrelationAndVariancesAndMSSS(self, hindcast, observation, flag, outputDir):
        """
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
        """
        hindcastStr = ''
        observationStr = ''
        for year in self.decadals:
            hindcastStr += ' ' + hindcast[year]
            observationStr += ' ' + observation[year]
        mergeHindcast = cdo.mergetime(input=hindcastStr,
                                      output=hindcast[self.decadals[0]]+self.getRandomStr()+'_merged.nc')
        mergeObservations = cdo.mergetime(input=observationStr,
                                          output=observation[self.decadals[0]]+self.getRandomStr()+'_merged.nc')

        #print mergeHindcast
        #remove seasonal cycle
        #mergeHindcast = self.removeSeasonalCycle(mergeHindcast)
        #mergeObservations = self.removeSeasonalCycle(mergeObservations)
        #print mergeHindcast
        
        #detrend data
        #mergeHindcast = self.detrendTimeSeries(mergeHindcast,keepMean=False)
        #mergeObservations = self.detrendTimeSeries(mergeObservations,keepMean=False)
        
        varianceHindcast = cdo.timstd(input=mergeHindcast,
                                      output=hindcast[self.decadals[0]]+self.getRandomStr()+'_variance')
        varianceObservation = cdo.timstd(input=mergeObservations,
                                         output=observation[self.decadals[0]]+self.getRandomStr()+'_variance')
        corr = cdo.timcor(input=' '.join([mergeHindcast, mergeObservations]), output=outputDir+flag+'_correlation.nc')
        
        bias = cdo.sub(input=corr+' '+ '-div ' + varianceHindcast + ' ' + varianceObservation,    
                       output=outputDir+flag+'_conditional_bias.nc')
        msss = cdo.sub(input='-sqr ' + corr + ' -sqr ' + bias,
                       output=outputDir+flag+'_msss.nc')

        if self.basic_output:
            return dict(msss= msss, bias= bias, corr= corr)
        else:
            # standard deviation ratio
            stdRatio = cdo.div(input=' '.join([varianceObservation, varianceHindcast]),
                               output=outputDir+flag+'_std_ratio.nc')
        
            # slope of conditional bias
            biasSlope = cdo.mul(input=' '.join([stdRatio, corr]),
                                output=outputDir+flag+'_biasslope.nc')
        
            # slope -1
            biasSlopeMinusOne = cdo.sub(input=' '.join([biasSlope, self.constantField]),
                                        output=outputDir+flag+'_biasslope-1.nc')
        
            return dict(msss=msss, bias=bias, corr=corr, biasSlope=biasSlope, stdRatio=stdRatio,
                        biasSlopeMinusOne=biasSlopeMinusOne)
        
    def getMSE(self, anomalies, obsAnom):
        """
        Method calculates the Mean Squared Error. 
        
        :note: method ist not used at the moment
        :param anomalies: dictionary with anomalies
        :param obsAnom: dict with observation anomalies
        :return: file with mse field
        """
        MSEList = list()
        for year in self.decadals:
            subTmp = cdo.sub(input=' '.join([anomalies[year], obsAnom[year]]))
            MSEList.append(cdo.sqr(input=subTmp)) 
        return cdo.ensmean(input=' '.join(MSEList), output=self.tmpDir+'mse.nc')
    
    def getRMSS(self, MSSS, flag, outputDir):
        """
        Calculates the RMSSS using the MSSS
        
        :param MSSS: netcdf file
        :param flag: name flag
        :param outputDir: output path 
        :return: file with RMSSS field 
        """
        tmp = cdo.sqrt(input="-sub " + self.constantField + " " + MSSS,                 
                       output=self.tmpDir+flag+'rmss_sqrt.nc')
        return cdo.sub(input=' '.join([self.constantField, tmp]), output=outputDir+flag+'_rmss.nc')
