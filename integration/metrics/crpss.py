'''
Created on 07.03.2013

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
import scipy.stats as stats
import numpy as np
from random import choice
import os
from cdo import *
cdo = Cdo()

from significance2 import Significance
from plotter import Plotter
from filehandler import FileHandler
from metricAbstract import MetricAbstract

from time import time

class CrpssError(Exception): pass
class NotEnoughEnsemblemembersFound(CrpssError): pass
class VariableName(CrpssError): pass

class Crpss(MetricAbstract):
    '''
    Class to calculate the CRPSS. Call method "analyze" to get results
    
    :param output: Path where the results are saved
    :param output_plot: Path where the plots are saved
    :param decadals: a list with decadal experiments, i.e. [1960,1960,...,1995]
    :param variable: CMOR-Name of variable like 'tas' (Near Surface Temperature)
    :param project: CMOR-Name of the project like CMIP5 or BASELINE1
    :param product1: CMOR-Parameter 
    :param institute1: CMOR-Name of the institute like MPI-M
    :param model: CMOR-Name of the used model like MPI-ESM-LR
    :param experiment: Experiment name like decs4e or historical
    
    :param leadtimes: Leadtimes to analyze, ie. 1,2-9
    :param observation: Observation or Reanalysis Experiment, HadCrut or ERA-Int or path to observation file        
    
    :param ensemblemembers: Comma separated string of esemble members like 'r1i1p1,r2i1p1,...'
      
    :param maskMissingValues: Boolean value. Whether you want to mask missing values or not
    :param result_grid: Griddescription of resultgrid like r72x36
    :param cache: Path for cachedir during analysis    
    :param baseDir: Path of Class to find default files
    :param timeFreq: Timefrequency of the files, ie 'mon' or 'day'
    :param basic_output: Boolean - If true only basic out is produced (CRPSS_ES and LESS).
    :param bootstrapSwitch: If you want calculating Significance using Bootstrap methods
    :param bootstrap_number: Number of Bootstrap runs    
    :param leatimes_mode: 'monthly' or 'yearly' calculation
    :param level: For 3D-Files --> Select a single level
    :param lonlatbox: If you want to select a specific lonlat-Box   
    :param fieldmean: Boolean - If you want to calculate field means
    :param zonalmean: Boolean - If you want to calculate zonal means
    
    :param input_part: String - Will be added to output path and cache file
    
    :param colormap: Matplotlib Colormap used for ploting
    :param observation_type: NOT USED ANYMORE!!!
    :param reanpath: NOT USED ANYMORE!!!
    
    :param months: Usefull if your analysing monthly and your experiment name has a month included
    '''

    outputFolderStruct = 'TIME_variable_project_model_product1_obsExp_YEARRANGE'.split('_')
    insideFolderNameFlag = 'project_product_model_experiment'.split('_')
    fileNameFlag = 'LEADTIMES_variable_project_product_model_experiment_YEARRANGE'.split('_')

    def __init__(   self, 
                    output='/tmp/crpss/output', 
                    output_plots='/tmp/crpss/plots', 
                    
                    decadals='1960,1965,1970,1975,1980,1985,1990,1995,2000',
                    variable='tas', 
                    
                    #CMOR params for input data
                    project='baseline0', 
                    product1='output', 
                    institute1='mpi-m', 
                    model = 'mpi-esm-lr',
                    experiment='decadal',
                    
                    leadtimes='1,2-9', 
                    observation = 'HadCrut',
                    observation_ensemble = '*',
                    ensemblemembers='*', 
                     
                    maskMissingValues = True,
                    result_grid=False, 
                    cache='/tmp/crpss/cache', 
                    baseDir = './', 
                    timeFreq = 'mon',  
                                        
                    basic_output = True,
		            #for bootstrapping
                    bootstrapSwitch=False,
                    bootstrap_number = 2,
                    
                    leadtimes_mode = 'yearly',
                    level=None, 
                    lonlatbox='', 
                    fieldmean=False, 
                    zonalmean=False,
                    
                    input_part = 'input1',
                    
                    #not really used
                    colormap='goddard',
                    observation_type = 'REANALYSIS', 
                    reanpath='', 
                    months = None,
                    ):
        '''
        Constructor. Set parameters and calls parent constructor (metricAbstarct)
        '''
        #Add kwargs to attributes
        listOfAttributes = ['decadals','variable',
                            'project','product1','institute1','model','experiment', #input data
                            'leadtimes','maskMissingValues','level','lonlatbox','fieldmean', 'zonalmean', 'observation_ensemble', 
                            'bootstrapSwitch','bootstrap_number', #for bootstrapping
                            'observation_type','colormap','timeFreq','basic_output','input_part']  #not used 
        for attr in listOfAttributes:
            setattr(self, attr, locals()[attr])

        if ensemblemembers != None and ensemblemembers != '*':
            self.ensemblemembers = ensemblemembers.split(',')
        else:
            self.ensemblemembers = '*'

        super(Crpss,self).__init__(tmpDir = self.checkPath(cache), output=output, output_plots=output_plots,
                                   baseDir=baseDir, result_grid=result_grid,observation=observation,
                                   leadtimes_mode=leadtimes_mode, decadals=decadals, months=months)


        
    def prepareInput(self):
        '''
        This method is to search and prepare the input files. 
        Integrated methods:
        1. Searching files
        2. Remapping
        3. Selecting lonlatbox
        4. Calculating fieldmean
        5. Calculating zonalmean 
        
        The prepared files are stored in "self.inputRemapped" and "self.observationRemapped"
        '''
        #TODO: find better way...
        self.findFiles.observation = self.observation  #WORKAROUND
        
        self.inputDict = dict()
        self.observationDict = dict()
        print 'Searching Files' 
        countYears= len(self.decadals)
        poolArgs = self.getPoolArgs(countYears,self.findFiles,list(self.decadals),self.project,self.model,self.variable,
                                    self.timeFreq,self.product1,self.ensemblemembers,self.institute1,
                                    self.experiment,self.maxLeadtime,self.minLeadtime,'getFiles')    
        self.inputDict = self.listToYeardict(self.multiProcess(poolArgs))
        for year in self.decadals:
            self.observationDict[year] = self.findFiles.getReanalysis(year, self.observation_type, self.obsExp, self.variable, maxLeadtime=self.maxLeadtime, observation_ensemble=self.observation_ensemble, minLeadtime=self.minLeadtime)

        self.inputRemapped = dict()
        self.observationRemapped = dict()
        obslist_tmp = list()
        print "Remapping Files"
        for year in self.decadals:
            self.inputRemapped[year] = self._remapFiles(self.inputDict[year])
            #self.observationRemapped[year] = self.remapFile(self.observationDict[year])  
            obslist_tmp.append(self.observationDict[year])    
        self.observationRemapped = self.listToYeardict(self._remapFiles(obslist_tmp))
        
        #More than 1 ensemblemember found?
        self.checkEnsembleError(self.inputDict)  

        if self.lonlatbox is not None:
            print 'Selecting lon-lat-box %s' %(self.lonlatbox)
            for year in self.decadals:
                self.inputRemapped[year] = self.sellonlatbox(self.inputRemapped[year])
                self.observationRemapped[year] = self._sellonlatbox(self.observationRemapped[year])  
        
        if self.fieldmean:
            print 'Calculating field mean'
            missmask = self.getMissingMaskForFieldMean()
            
            for year in self.decadals:
                #Apply Missing value Mask to all fields
                self.inputRemapped[year] = self.applyMissingMaskForFieldMean(self.inputRemapped[year], missmask)
                self.observationRemapped[year] = self._applyMissingMaskForFieldMean(self.observationRemapped[year], missmask)
                
                self.inputRemapped[year] = self.fieldMean(self.inputRemapped[year])
                self.observationRemapped[year] = self._fieldMean(self.observationRemapped[year]) 
                
        
        #self.zonalmean = True
        if self.zonalmean:
            print 'Calculating zonal mean'
            level_intersection = self.getLevelIntersection(self.observationRemapped[self.decadals[0]], self.inputRemapped[self.decadals[0]][0])
            print 'Common levels: %s' % (','.join(level_intersection))
            for year in self.decadals:
                self.observationRemapped[year] = cdo.zonmean(input=self.observationRemapped[year], output=self.observationRemapped[year]+'_zonmean')
                self.observationRemapped[year] = cdo.sellevel(','.join(level_intersection), input=self.observationRemapped[year], output=self.observationRemapped[year]+'sellevel')
                
                self.inputRemapped = self.multiProcessCdo(self.inputRemapped[year], 'zonmean')
                self.inputRemapped[year] = self.multiProcessCdo(self.inputRemapped[year], 'sellevel', ','.join(level_intersection))
#                
#                tmpList = list()
#                for fn in self.inputRemapped[year]:
#                    tmp_file = cdo.zonmean(input=fn,output=fn+'_zonmean')
#                    tmpList.append(cdo.sellevel(','.join(level_intersection), input=tmp_file, output=tmp_file+'sellevel'))
#                self.inputRemapped[year] = tmpList

        
    def analyze(self):
        '''
        Main function to calculate the CRPSS after Goddard et al. (2012)
        It also calculates the ESS and LESS after Kadow et al. (2015)
        The following steps are performed:
        
        1. Temporal Smoothing
        2. Calculating ensemble means of hindcasts
        3. Calculating cross-validated means of ensemblemembers
        4. Calculating anomalies of hindcast and observations
        5. Removing conditional bias
        6. Calculating mean ensemble variance and reference STD
        7. Calculating climatological standard deviation
        8. Calculating CRPS for mean ensemble variance and for reference STD
        9. Calculating CRPSS
        10. Calculating Ensemble Spread Score (ESS) and Logarithmic ESS
        11. Plotting        
        '''
        self.constant = self.createConstantFile(self.gridFile)
        #check if prepare output was called or variables are set
        try:
            self.inputRemapped
            self.observationRemapped
        except AttributeError:
            raise CrpssError, 'Please run "prepareIntput() first or provide input data'

        #get missing value mask
        if not self.fieldmean and self.maskMissingValues:
            self.misvalMask = self.calcMissingValueMask(self.observationRemapped)
        else:
            self.misvalMask = self.constant
        #escape attributes 
        self.escapeAttributes()

        
        ##########################
        ## START YEARRANGE LOOP
        ##########################
        temporalSmoothing = self.getYearRange()
        for yearRange in temporalSmoothing:
            startyear = yearRange[0]
            endyear = yearRange[1]
            if startyear == endyear:
                print '--------------\nAnalyzing leadtime %s\n--------------' %(startyear)
            else:
                print '--------------\nAnalyzing leadtime %s to %s\n--------------' %(startyear, endyear)
            #print '--------------\nAnalyzing year %s to %s\n--------------' % (startyear, endyear)
            rangeStr = '%s-%s' % (startyear, endyear)  
            tempSmoothedEns = dict()
            tempSmoothedObs = dict()
            
            tmp_outdir = self.outputDir
            tmp_outplot = self.outputPlots 
            flag1 = self.constructName(self.insideFolderNameFlag,exp='1',extra=self.input_part)
            fnFlag1 = self.constructName(self.fileNameFlag, exp='1', startYear=str(startyear), endYear=str(endyear))
            self.outputDir = self.makeFolder(self.outputDir+'%s/%s/ensemble_spread'%(rangeStr,flag1))
            self.outputPlots = self.makeFolder(self.outputPlots+'%s/%s/ensemble_spread'%(rangeStr,flag1))
            
            print "Temporal Smoothing"
            for year in self.decadals:
                tempSmoothedEns[year] = self.tempSmoothing(self.inputRemapped[year], startyear, endyear)
                tempSmoothedObs[year] = self._tempSmoothing(self.observationRemapped[year], startyear, endyear)
            
            print 'Calculating Ensemble Mean'
            ensMeanList = self.getEnsembleMean(tempSmoothedEns)
            
            print 'Calculating cross-validated mean'
            crossValMean = self.getCrossValMean(ensMeanList, 'model'+rangeStr)
            if(self.obsExp != 'HadCrut'):
                obsCrossValMean = self.getCrossValMean(tempSmoothedObs, 'obs')
                tempSmoothedObs = self.getAnomalies(tempSmoothedObs, obsCrossValMean)
            
            print 'Calculating anomalies'
            anomalies = self.getAnomalies(ensMeanList, crossValMean)  
            #self.printValues(anomalies, 'EnsAnoms')      
            #Calculate Anomalies for Ensemble
            ensList,crossList = list(),list()
            for year in self.decadals:
                ensList.append(tempSmoothedEns[year])
                crossList.append(crossValMean[year])

            poolArgs = self.getPoolArgs(len(ensList),self,ensList,crossList,'_getAnomalies')
            ensAnoms = self.listToYeardict(self.multiProcess(poolArgs))
            #self.printValues(ensAnoms, 'EnsAnoms')
            tmp2 = cdo.ensmean(input=' '.join(tempSmoothedObs.values()))
            #Calculate obs anomalies again because we use only a few points
            for year in self.decadals:
                tempSmoothedObs[year] = cdo.sub(input=' '.join([tempSmoothedObs[year], tmp2]), output=tempSmoothedObs[year]+'meanAgain'+str(startyear)+str(endyear))
            
            print 'Removing Conditional Bias'
            ensAnoms = self.removeConditionalBiasNew(anomalies, tempSmoothedObs, ensAnoms)
            anomalies = self.getEnsembleMean(ensAnoms)

            print 'Calculating ensemble variance'
            ensembleVariance = self.getEnsembleStd(ensAnoms)  
#            ensembleVariance = self.getEnsembleStd(tempSmoothedEns) 
#            ensembleVariance = self.getEnsembleStd(self.inputRemapped) 
            #print ensembleVariance
            print 'Calculating MSE'
            ensembleVarianceRef = self.getReferenceStd(anomalies, tempSmoothedObs)
            self.printValues(ensembleVarianceRef, 'MSE')
            
            print 'Calculating climatological standard deviation'
            norm = np.sqrt(len(tempSmoothedObs.values())/(len(tempSmoothedObs.values())-1.))
            climStd = cdo.mulc(norm, input=cdo.ensstd(input=' '.join(tempSmoothedObs.values())),
                                     output=self.tmpDir+rangeStr+'obsSTD_NORM.nc')
            ensListSelYearTimmean,obsSelYearTimmean,ensvarSelYear,ensvarRefSelYear,climStdSelYear = dict(),dict(),dict(),dict(),dict()
            #get obs mean
            observationMean = cdo.ensmean(input=' '.join(tempSmoothedObs.values()), 
                                          output=self.tmpDir+rangeStr+'observation_mean')
            obsMeanList = dict()
            #rename a few variables 
            for year in self.decadals:
                ensListSelYearTimmean[year] = anomalies[year]
                obsSelYearTimmean[year] = tempSmoothedObs[year]
                obsMeanList[year] = observationMean
                ensvarSelYear[year] = ensembleVariance
                ensvarRefSelYear[year] = ensembleVarianceRef
                climStdSelYear[year] = climStd
            biasCorrected = ensListSelYearTimmean
            
            print 'Calculating CRPS'
            crpsEns = self.getCrps(biasCorrected, obsSelYearTimmean, ensvarSelYear, 'ens')
            crpsRef = self.getCrps(biasCorrected, obsSelYearTimmean, ensvarRefSelYear, 'ref')
            crpsClim = self.getCrps(obsMeanList, obsSelYearTimmean, climStdSelYear, 'clim')
            
            print 'Calculating CRPSS'
            #print ensvarSelYear
            #print ensvarRefSelYear
            ensemblespreadscore = cdo.div(input=' '.join([ensembleVariance,ensembleVarianceRef]), output=self.tmpDir+fnFlag1+'_ensspread_vs_referror_tmp.nc')
            ensemblespreadscore = cdo.div(input=' '.join([ensembleVariance,ensembleVarianceRef]), output=self.outputDir+fnFlag1+'_ensspread_vs_referror.nc')
            ensemblespreadscore_ln_tmp = cdo.ln(input=cdo.div(input=' '.join([ensembleVariance,ensembleVarianceRef]), output=self.tmpDir+fnFlag1+'_ensspread_vs_referror.nc'), output = self.tmpDir+fnFlag1+'_ensspread_vs_referror_ln.nc',options='-b 64')
            #We do this because ncview can't handly (-)Infinity values
	    ensemblespreadscore_ln = cdo.setvals('-Infinity,1e20',input=ensemblespreadscore_ln_tmp, output=self.outputDir+fnFlag1+'_ensspread_vs_referror_ln.nc')

	    crpss =  self.getCrpss(crpsEns, crpsRef, '_'.join([fnFlag1,'ens-vs-ref']))
            
            if not self.basic_output:
                ensemblespreadscore = cdo.div(input=' '.join([ensembleVariance,ensembleVarianceRef]), output=self.outputDir+fnFlag1+'_ensspread_vs_referror.nc')
                crpssEC =  self.getCrpss(crpsEns, crpsClim, '_'.join([fnFlag1,'ens-vs-clim']))
                crpssRC =  self.getCrpss(crpsRef, crpsClim, '_'.join([fnFlag1,'ref-vs-clim']))
            
            
            
            #apply missing value mask
            filesToPlot = self.applyMissingValueMask(str(startyear), str(endyear), [self.outputDir])
            #print filesToPlot
            if self.maskMissingValues:
                if self.basic_output:
                    ensemblespreadscore = cdo.mul(input=' '.join([ensemblespreadscore,self.misvalMask]),output=ensemblespreadscore+'miss_val')
                    crpss = [s for s in filesToPlot if 'ens-vs-ref_crpss.nc' in s][0]  #filesToPlot[2]
                    ensemblespreadscore_ln = [s for s in filesToPlot if 'ensspread_vs_referror_ln' in s][0] #filesToPlot[1]
                else:
                    ensemblespreadscore = filesToPlot[0]
                    ensemblespreadscore_ln = [s for s in filesToPlot if 'ensspread_vs_referror_ln' in s][0] #filesToPlot[1]
                    crpss = [s for s in filesToPlot if 'ens-vs-ref_crpss.nc' in s][0] #filesToPlot[2]
                    crpssEC = filesToPlot[3]
                    crpssRC = filesToPlot[4]
            if self.fieldmean:
                try:
                    fieldmean_files.append(filesToPlot)
                except:
                    fieldmean_files = [filesToPlot]
                if self.bootstrapSwitch:
		    b_ens_ref = self.bootstrap(crpsEns, crpsRef, 'ens-vs-ref_'+rangeStr, crpss, self.bootstrap_number, plot_range=[-0.5,0.], colorbar='goddard')
                    if not self.basic_output:
                        b_ens_clim = self.bootstrap(crpsEns, crpsClim, 'ens-vs-clim_'+rangeStr, crpssEC, self.bootstrap_number)
                        b_ref_clim = self.bootstrap(crpsRef, crpsClim, 'ref-vs-clim_'+rangeStr, crpssRC, self.bootstrap_number)
                    self.bootstrapEnsembleSpreadscore(ensAnoms, anomalies, tempSmoothedObs, ensemblespreadscore, self.bootstrap_number, ensemblespreadscore_ln)
            elif self.zonalmean:
                if self.bootstrapSwitch:
                    print 'Bootstrapping'
                    b_ens_ref = self.bootstrap(crpsEns, crpsRef, 'ens-vs-ref_'+rangeStr, crpss, self.bootstrap_number, plot_range=[-0.5,0.], colorbar='goddard')
                    if not self.basic_output:
                        b_ens_clim = self.bootstrap(crpsEns, crpsClim, 'ens-vs-clim_'+rangeStr, crpssEC, self.bootstrap_number)
                        b_ref_clim = self.bootstrap(crpsRef, crpsClim, 'ref-vs-clim_'+rangeStr, crpssRC, self.bootstrap_number)
                    self.bootstrapEnsembleSpreadscore(ensAnoms, anomalies, tempSmoothedObs, ensemblespreadscore, self.bootstrap_number, ensemblespreadscore_ln)
                else:
                    Plotter.plotVerticalProfile(crpss, -0.5, 0., 'goddard', lonlatbox=self.lonlatbox)
                    Plotter.saveFig(self.outputPlots, self.extractFilename(crpss))
                    Plotter.plotVerticalProfile(ensemblespreadscore_ln, -1, 1, 'RdBu_r', lonlatbox=self.lonlatbox)
                    Plotter.saveFig(self.outputPlots, self.extractFilename(ensemblespreadscore_ln))
                    if not self.basic_output:
                        Plotter.plotVerticalProfile(crpssEC, -1, 1, 'RdBu_r', lonlatbox=self.lonlatbox)
                        Plotter.saveFig(self.outputPlots, self.extractFilename(crpssEC))
                        Plotter.plotVerticalProfile(crpssRC, -1, 1, 'RdBu_r', lonlatbox=self.lonlatbox)
                        Plotter.saveFig(self.outputPlots, self.extractFilename(crpssRC))
                        Plotter.plotVerticalProfile(ensemblespreadscore, 0, 2, 'RdBu_r', lonlatbox=self.lonlatbox)
                        Plotter.saveFig(self.outputPlots, self.extractFilename(ensemblespreadscore))
            else:
                if self.bootstrapSwitch:
                    print 'Bootstrapping'
                    b_ens_ref = self.bootstrap(crpsEns, crpsRef, 'ens-vs-ref_'+rangeStr, crpss, self.bootstrap_number, plot_range=[-0.5,0.], colorbar='goddard')
                    if not self.basic_output:
                        b_ens_clim = self.bootstrap(crpsEns, crpsClim, 'ens-vs-clim_'+rangeStr, crpssEC, self.bootstrap_number)
                        b_ref_clim = self.bootstrap(crpsRef, crpsClim, 'ref-vs-clim_'+rangeStr, crpssRC, self.bootstrap_number)
                    self.bootstrapEnsembleSpreadscore(ensAnoms, anomalies, tempSmoothedObs, ensemblespreadscore, self.bootstrap_number, ensemblespreadscore_ln)
    #                Plotter.plotField(ensemblespreadscore, 0, 2, 'RdBu_r', lonlatbox=self.lonlatbox)
    #                Plotter.saveFig(self.outputPlots, self.extractFilename(ensemblespreadscore))
                else:
                    Plotter.plotField(crpss, -0.5, 0., 'goddard', lonlatbox=self.lonlatbox)
                    Plotter.saveFig(self.outputPlots, self.extractFilename(crpss))
                    Plotter.plotField(ensemblespreadscore_ln, -1, 1, 'RdBu_r', lonlatbox=self.lonlatbox)
                    Plotter.saveFig(self.outputPlots, self.extractFilename(ensemblespreadscore_ln))
                    if not self.basic_output:
                        Plotter.plotField(crpssEC, -1, 1, 'RdBu_r', lonlatbox=self.lonlatbox)
                        Plotter.saveFig(self.outputPlots, self.extractFilename(crpssEC))
                        Plotter.plotField(crpssRC, -1, 1, 'RdBu_r', lonlatbox=self.lonlatbox)
                        Plotter.saveFig(self.outputPlots, self.extractFilename(crpssRC))
                        Plotter.plotField(ensemblespreadscore, 0, 2, 'RdBu_r', lonlatbox=self.lonlatbox)
                        Plotter.saveFig(self.outputPlots, self.extractFilename(ensemblespreadscore))

                
            self.outputDir = tmp_outdir
            self.outputPlots = tmp_outplot
    
        if self.fieldmean:
            outputPlots = self.outputPlots+'ensemble_spread/'
            self.makeFolder(outputPlots)
            flag1 = self.constructName(self.fileNameFlag, exp='1', startYear='1', endYear='1')
            flag1 = flag1[4:]
            plot_list = [('ensspread_vs_referror_ln','Logarithmic Ensemble Spread Score',[-1,1]),('ens-vs-ref_crpss','Continous Ranked Probabillity Score',[-0.5,0])]
            if self.bootstrapSwitch:
                new_list = list()
                for item in fieldmean_files:
                    new_list += item
                Plotter.plotLeadtimeseriesSign(new_list,[flag1],plot_list)
            else:
                Plotter.plotLeadtimeseries(fieldmean_files,[flag1],plot_list)
            Plotter.saveFig(outputPlots, flag1+'ensemble_spread_leadtimeseries')
    
    def removeConditionalBiasNew(self, hindcast, observations, ensembleMembers):
        '''
        Remove the conditional bias from the ensemble members. Following Murphy 1988 and Goddard et al. 2012 
        
        :param hindcast: dict of hindcast ensemble means
        :param observations: dict of observation anomalies
        :param ensembleMembers: dict of hindcast ensemble members
        :return: dict of bias corrected ensemble members
        '''
        hindcastStr = ''
        observationsStr = ''
        for year in self.decadals:
            tmp = cdo.yearmean(input=hindcast[year])
            hindcastStr += ' ' + tmp
            observationsStr += ' ' + observations[year]
        mergeHindcast = cdo.mergetime(input=hindcastStr, output=hindcast[self.decadals[0]]+'_merged.nc')
        mergeObservations = cdo.mergetime(input=observationsStr, output=observations[self.decadals[0]]+'_merged.nc')
        varHindcast = cdo.timstd(input=mergeHindcast, output=hindcast[self.decadals[0]]+'_variance')
        self.printValues(varHindcast, 'Model STD')
        varObservations = cdo.timstd(input=mergeObservations, output=observations[self.decadals[0]]+'_variance')
        self.printValues(varObservations, 'Observation STD')
        corr = cdo.timcor(input=' '.join([mergeHindcast, mergeObservations]), output=self.tmpDir+'_correlation.nc')
        self.printValues(corr, 'Correlation')
        #######################################
        ###This is were Goddard et al. made the mistake
        ###they flipped the quotient varObservations/varHindcast and calculated varHindcast/varObservations
        ###Right: varObservations/varHindcast
        ###Wrong: varHindcast/varObservations
        #######################################
        bias = cdo.mul(input=corr+' '+ cdo.div(input=varObservations+' '+varHindcast), output=self.tmpDir+'conditional_biasNew.nc')
        self.printValues(bias, 'Conditional Bias')     
        hindList = list()
        ensembleMembersList = list()
        for year in self.decadals:
            hindList.append(hindcast[year])
            ensembleMembersList.append(ensembleMembers[year])

        poolArgs = self.getPoolArgs(len(hindList),self,hindList,bias,ensembleMembersList,'_removeConditionalBias')
        biasCorrected = self.listToYeardict(self.multiProcess(poolArgs))
        return biasCorrected
    
    def _removeConditionalBias(self, hindcast, bias, ensembleMembers):
        '''
        This is the single process version of removeConditionalBias
        '''
        shift = cdo.mul(input=' '.join([hindcast, cdo.subc('1', input=bias, output=hindcast+'_bias_shift')]),
                        output=hindcast+'_shift')
        newList = list()
        for f in ensembleMembers:
            newList.append(cdo.add(input=' '.join([f, shift]), output=f+'_BIASCORRECTED2'))
        return newList  
        
            
    def getCrps(self, hindcast, observation, variances, tag):
        '''
        Calculates crpss for dictionarys after Gneiting and Raferty (2007)
        
        :param hindcast: dict of anomaly conditional bias corrected hindcast
        :param observation: dict of anomaly observation
        :param variances: dict of Variances you want to test
        :param tag: name tag
        :return: dict of crps        
        '''
        hindList,obsList,varList,tagList = list(),list(),list(),list()
        for year in self.decadals:
            hindList.append(hindcast[year])
            obsList.append(observation[year])
            varList.append(variances[year])
            tagList.append(tag+str(year))

        poolArgs = self.getPoolArgs(len(hindList),self,hindList,obsList,varList,tagList,'_calcCrps')
        crps = self.listToYeardict(self.multiProcess(poolArgs))
         
        return crps

    def _calcCrps(self, hindcast, observation, variance, tag):
        '''
        Method calculates the CRPS after the formula by Gneiting and Raferty (2007).
        
        :param hindcast: anomaly conditional bias corrected hindcast
        :param observation: anomaly observation
        :param variance: Variance or STD you want to test
        :param tag: name tag
        :return: crps file
        '''
        np.seterr('ignore')
        try:
            x = cdo.div(input=' '.join([cdo.sub(input=' '.join([observation, hindcast]), output=hindcast+'_min_hindcast' + tag, options = '-f nc'), variance]),
                        output=hindcast+tag+'_div_')#, returnMaArray=cdo.showname(input=hindcast)[0])
	    x = FileHandler.openNetCDFFile(x,mode='var')
        except KeyError:
            raise VariableName, 'Could not find variablename "%s" in observation files.' % (cdo.showname(input=hindcast)[0])
            
        xCopy = np.ma.masked_greater(x, 0.8e20)
        varVar = FileHandler.openNetCDFFile(variance, mode='var')
        
        x = x.squeeze()
        xCopy = xCopy.squeeze()
        crps = - varVar * (1/np.sqrt(np.pi) - 2. * stats.norm.pdf(x) - x * (2. * stats.norm.cdf(x) - 1.))
        crps = np.ma.array(crps, mask=xCopy.mask, fill_value=1e20)
        crps = crps.filled(1e20)
        crpsFile = FileHandler.saveToNetCDF(crps, hindcast, 'crps'+tag) 
        return crpsFile
    
    def getCrpss(self, crps1, crps2, tag):
        '''
        Calculates the crpss out of two different crps dicts. 
        
        :param crps1: dict
        :param crps2: dict
        :param tag: name tag
        :return: crpss filepath
        '''
        crpsSum1 = cdo.enssum(input=' '.join(crps1.values()), output=self.tmpDir+'crpsSum1'+tag)
        crpsSum2 = cdo.enssum(input=' '.join(crps2.values()), output=self.tmpDir+'crpsSum2'+tag)
        crpsQuot = cdo.div(input=' '.join([crpsSum1, crpsSum2]), output=crpsSum1+'div')
        crpsQuotMinus = cdo.mulc('-1',input=crpsQuot,output=crpsQuot+'minus')
        crpss = cdo.addc('1',input=crpsQuotMinus,output=self.outputDir + tag + '_crpss.nc')
#        crpss = cdo.sub(input=' '.join([self.constant, cdo.div(input=' '.join([crpsSum1, crpsSum2]), output=crpsSum1+'div')]),
#                        output = self.outputDir + tag + '_crpss.nc')
        return crpss
              
    def getReferenceStd(self, hindcast, observation):
        '''
        Calculate Reference STD after Goddard et al. Calculates the Root Mean Squared Error between hindcast and observation 
        and standardizes by n-2
        
        :param hindcast: dict of anomalies
        :param observation: dict of observation anomalies
        :return: reference STD
        '''
        varSum = list()
        for year in self.decadals:
            varSum.append(cdo.sqr(input=cdo.sub(input=' '.join([hindcast[year], observation[year]]), output=hindcast[year]+self.getRandomStr()),
                                  output=hindcast[year]+'-obs_'+str(year)+'_SQR'))
        rms = cdo.enssum(input=' '.join(varSum),output=self.tmpDir+'enssum_rms'+self.getRandomStr())
        self.printValues(rms, 'RMS')
        return cdo.sqrt(input=cdo.divc((len(varSum)-2), input=rms, output=hindcast[year]+self.getRandomStr()+'divc'), output=hindcast[year]+self.getRandomStr()+'refStd')
   
   
    def bootstrap(self, crps, crpsRef, tag, crpss, bootstrap_number=500, plot_range=[-1,1], colorbar='RdBu_r'):
        '''
        Bootstrap function to calculate significance crosses for crpss and other metrics. 
        It also plots the fields with significance crosses
        '''
        tmpDir = self.outputDir
        self.outputDir = self.tmpDir + 'bootstrap_' + tag + '/'
        if not os.path.isdir(self.outputDir):
            os.mkdir(self.outputDir)

        poolArgs = self.getPoolArgs(bootstrap_number,self,crps,crpsRef,tag,crpss,range(0,bootstrap_number),'_multiBootstrap')
        bootstrapList = self.multiProcess(poolArgs)
            
        self.outputDir = tmpDir      
        
        significance = Significance(self.tmpDir, self.outputPlots)
        if self.fieldmean:
            (sig_lon, sig_lat) = significance.checkSignificanceFldmean(bootstrapList, crpss)
        else:
            (sig_lon, sig_lat) = significance.checkSignificance(bootstrapList, crpss)
            if self.zonalmean:
                m = Plotter.plotVerticalProfile(crpss, plot_range[0], plot_range[1], colorbar, lonlatbox=self.lonlatbox)
                Plotter.addCrossesXY(m, crpss+'_significance_mask')
    
            else:
                m = Plotter.plotField(crpss, plot_range[0], plot_range[1], colorbar, lonlatbox=self.lonlatbox)
                Plotter.addCrosses(m, sig_lon, sig_lat)
            Plotter.saveFig(self.outputPlots, crpss.split('/')[-1])
            
        return bootstrapList   
    
    def _multiBootstrap(self, crps, crpsRef, tag, crpss, i):
        '''
        Single process version used in crpss.bootstrap
        '''
        tmp_crps = dict()
        tmp_crps_ref = dict()
        for year in self.decadals:
            choice_year = choice(self.decadals)
            tmp_crps[year] = crps[choice_year]
            tmp_crps_ref[year] = crpsRef[choice_year]
                
        return self.getCrpss(tmp_crps, tmp_crps_ref, tag+'bootstrap_number_'+str(i)) 
    
    def _multiEnsspreadBootstrap(self, ensAnoms, anomalies, tempSmoothedObs, i):
        '''
        TODO: Not working!!!
        '''
        tmpDir = self.tmpDir + 'bootstrap_ensspread/'
        tmp_ensAnoms = dict()
        tmp_anoms = dict()
        tmp_obs = dict()
        for year in self.decadals:
            choice_year = choice(self.decadals)
            tmp_ensAnoms[year] = ensAnoms[choice_year]
            tmp_anoms[year] = anomalies[year]
            tmp_obs[year] = tempSmoothedObs[year]
        
        tmp_ensspread = self.getEnsembleStd(tmp_ensAnoms)
        tmp_rmse = self.getReferenceStd(tmp_anoms, tmp_obs)
        tmp_spreadscore = cdo.div(input=' '.join([tmp_ensspread, tmp_rmse]),
                                  output=tmpDir+str(i)+'spreadscore_bootstrap')   
        
        return tmp_spreadscore 
    
    def bootstrapEnsembleSpreadscore(self, ensAnoms, anomalies, tempSmoothedObs, spreadscore, bootstrap_number, spreadscore_ln):
        '''
        Bootstrapping method for ESS and LESS
        '''
        print 'Bootstrapping spreadscore'
        tmpDir = self.tmpDir + 'bootstrap_ensspread/'
        if not os.path.isdir(tmpDir):
            os.mkdir(tmpDir)
        
        bootstrapList = []       
        for i in range(0,bootstrap_number):
            tmp_ensAnoms = dict()
            tmp_anoms = dict()
            tmp_obs = dict()
            for year in self.decadals:
                choice_year = choice(self.decadals)
                tmp_ensAnoms[year] = ensAnoms[choice_year]
                tmp_anoms[year] = anomalies[year]
                tmp_obs[year] = tempSmoothedObs[year]
            
            tmp_ensspread = self.getEnsembleStd(tmp_ensAnoms)
            tmp_rmse = self.getReferenceStd(tmp_anoms, tmp_obs)
            tmp_spreadscore = cdo.div(input=' '.join([tmp_ensspread, tmp_rmse]),
                                      output=tmpDir+str(i)+'spreadscore_bootstrap')        
            bootstrapList.append(tmp_spreadscore)
   
        significance = Significance(self.tmpDir, self.outputPlots)
        if self.fieldmean:
            (sig_lon, sig_lat) = significance.checkSignificanceFldmean(bootstrapList, spreadscore, check_value=1)
            cdo.ln(input=spreadscore+'_bootstrap_max_val', output=spreadscore_ln+'_bootstrap_max_val')
            cdo.ln(input=spreadscore+'_bootstrap_min_val', output=spreadscore_ln+'_bootstrap_min_val')
        else:
            (sig_lon, sig_lat) = significance.checkSignificance(bootstrapList, spreadscore, check_value=1)
            if not self.basic_output:
                if self.zonalmean:
                    m = Plotter.plotVerticalProfile(spreadscore, 0, 2, 'RdBu_r', lonlatbox=self.lonlatbox)
                    Plotter.addCrossesXY(m, spreadscore+'_significance_mask')
                else:
                    m = Plotter.plotField(spreadscore, 0, 2, 'RdBu_r', lonlatbox=self.lonlatbox)
                    Plotter.addCrosses(m, sig_lon, sig_lat)
                Plotter.saveFig(self.outputPlots, spreadscore.split('/')[-1])
            
            if self.zonalmean:
                m = Plotter.plotVerticalProfile(spreadscore_ln, -1, 1, 'RdBu_r', lonlatbox=self.lonlatbox)
                Plotter.addCrossesXY(m, spreadscore+'_significance_mask')
            else:
                m = Plotter.plotField(spreadscore_ln, -1, 1, 'RdBu_r', lonlatbox=self.lonlatbox)
                Plotter.addCrosses(m, sig_lon, sig_lat)
            Plotter.saveFig(self.outputPlots, spreadscore_ln.split('/')[-1])
         
        return bootstrapList  
        
    
    def checkEnsembleError(self, ensList):
        '''
        Checks if enough ensemblemembers are selected for analysis
        Crpss needs at least 2 ensemblemembers. And is only reasonable for more
        '''
        for year in self.decadals:
            if len(ensList[year]) < 2:
                raise NotEnoughEnsemblemembersFound, 'Not enough ensemblemembers found for %s %s %s %s %s' %(self.variable, self.project, self.model, self.project, year)
            
        
         
       
