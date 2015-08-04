'''
Created on 16.08.2013

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

from msss import Msss
from significance2 import Significance
from plotter import Plotter

from cdo import *
cdo = Cdo()
import os
from random import choice
import shutil
    
    
class MsssBootstrap(Msss):
    '''
    Class with special method for bootstrapping
    '''
    def _prepareBootstrap(self, year):
        '''
        Searches the files for specific year, remaps them, and calculates ensemblemean
        
        :param year
        :return ensmeanProject1, ensmeanProject2
        '''
        if(hasattr(self, 'ensList')):
            ensListProject1 = self.ensList[year]
        else:
            ensListProject1 = self.findFiles.getFiles(year, self.project1, self.model1, self.variable, 
                                                    ensemblemembers=self.ensemblemembers1, product=self.product1, 
                                                    institute=self.institute1, exp_prefix=self.experiment1,
                                                    maxLeadtime=self.maxLeadtime,minLeadtime=self.minLeadtime)
        tmpList = list()
        for ensfile in ensListProject1:
            #yearmean
            #tmpList.append(cdo.yearmean(input=ensfile, output=self.tmpDir+self.extractFilename(ensfile)+str(year)+self.project1))
            tmpList.append(ensfile)
            #remap
            tmpList2 = list()
            for f in tmpList:
                tmpList2.append(self.remapFile(f))
        tmpOut = self.tmpDir+self.project1+'_'+self.product1+'_'+self.model1+self.experiment1+'/'+str(year)+self.getRandomStr()+'/'
        self.makeFolder(tmpOut)
        
        ensmeanProject1 = cdo.ensmean(input=' '.join(tmpList2), output=tmpOut+str(year)+self.experiment1+self.project1+'_'+self.model1+'_'+self.product1+self.getRandomStr()+'.nc')
        
        if(hasattr(self, 'histList')):
            ensListProject2 = self.histList[year]
        else:
            ensListProject2 = self.findFiles.getFiles(year, self.project2, self.model2, self.variable,
                                                    ensemblemembers=self.ensemblemembers2, product=self.product2, 
                                                    institute=self.institute2, exp_prefix=self.experiment2,
                                                    maxLeadtime=self.maxLeadtime,minLeadtime=self.minLeadtime)
        tmpList = list()
        for ensfile in ensListProject2:
            #yearmean
            #tmpList.append(cdo.yearmean(input=ensfile, output=self.tmpDir+self.extractFilename(ensfile)+str(year)+self.project2))
            tmpList.append(ensfile)
            #remap
            tmpList2 = list()
            for f in tmpList:
                tmpList2.append(self.remapFile(f))
        tmpOut = self.tmpDir+self.project2+'_'+self.product2+'_'+self.model2+self.experiment2+'/'+str(year)+self.getRandomStr()+'/'
        #os.makedirs(tmpOut)
        self.makeFolder(tmpOut)
        ensmeanProject2 = cdo.ensmean(input=' '.join(tmpList2), output=tmpOut+str(year)+self.experiment2+self.project2+'_'+self.model2+'_'+self.product2+self.getRandomStr()+'.nc')
        
        return (ensmeanProject1, ensmeanProject2,)


    def prepareBootstrap(self):
        '''
        Multiprocess approach
        Uses _prepareBootstrap to prepare the Bootstrap files for input1 and input2
        
        :return: ensMeanProject1Dict, ensMeanProject2Dict
        '''
        #print self.findFiles
        procs = len(self.decadals)
        poolArgs =  zip([self]*procs, self.decadals, ['_prepareBootstrap'] * procs)
        result = self.multiProcess(poolArgs)
        ensMeanProject1Dict = dict()
        ensMeanProject2Dict = dict()
        for i in range(0, len(result)):    
            ensMeanProject1Dict[self.decadals[i]] = result[i][0]
            ensMeanProject2Dict[self.decadals[i]] = result[i][1]
        
        if self.lonlatbox is not None:
            print 'Selecting lon-lat-box %s' %(self.lonlatbox)
            for year in self.decadals:
                ensMeanProject1Dict[year] = self.sellonlatbox(ensMeanProject1Dict[year])
                ensMeanProject2Dict[year] = self.sellonlatbox(ensMeanProject2Dict[year]) 

        if self.fieldmean:
            print 'Calculating field mean'
            #Calculate missing value fields     
            for year in self.decadals:
                #Apply Missing value Mask to all fields                                          output=self.observationRemapped[year]+'_masked' )
                #self.observationRemapped[year] = self._applyMissingMaskForFieldMean(self.observationRemapped[year], missmask)
                ensMeanProject1Dict[year] = self._applyMissingMaskForFieldMean(ensMeanProject1Dict[year], self.obsmissmask)
                ensMeanProject2Dict[year] = self._applyMissingMaskForFieldMean(ensMeanProject2Dict[year], self.obsmissmask)
                
                ensMeanProject1Dict[year] = self._fieldMean(ensMeanProject1Dict[year])
                ensMeanProject2Dict[year] = self._fieldMean(ensMeanProject2Dict[year])
        
        if self.zonalmean:
            print 'Calculating zonal mean'
            level_intersection = self.getLevelIntersection(self.observationRemapped[self.decadals[0]], self.input1Remapped[self.decadals[0]][0])
            for year in self.decadals:
                tmpList = list()
                #for fn in ensMeanProject1Dict[year]:
                tmp_file = cdo.zonmean(input=ensMeanProject1Dict[year],output=ensMeanProject1Dict[year]+'_zonmean')
                #tmpList.append()
                ensMeanProject1Dict[year] = cdo.sellevel(','.join(level_intersection), input=tmp_file, output=tmp_file+'sellevel')
                tmpList = list()
                #for fn in ensMeanProject2Dict[year]:
                tmp_file = cdo.zonmean(input=ensMeanProject2Dict[year],output=ensMeanProject2Dict[year]+'_zonmean')
                #tmpList.append()
                ensMeanProject2Dict[year] = cdo.sellevel(','.join(level_intersection), input=tmp_file, output=tmp_file+'sellevel')
        
        self.bootstrapPoolProject1 = ensMeanProject1Dict
        self.bootstrapPoolProject2 = ensMeanProject2Dict
    
        return (ensMeanProject1Dict, ensMeanProject2Dict,)
    
    def bootstrapGoddard(self, outputFolder='./tmp/'):
        '''
        Actual method for selecting files of a bootstrap run.
        Select input data from bootstrap pool by choice from bootstrap pool
        together with corresponding obs data
        --> constraint is to select within a 5 Year range only (because of a trend in the data) 
        
        :param outputFolder: temp folder 
        '''
        if(not os.path.isdir(outputFolder)):
            os.makedirs(outputFolder)
        
        self.tmpDir = outputFolder
        
        newPoolList1 = dict()
        newPoolList2 = dict()
        for year in self.decadals:
            f = self.bootstrapPoolProject1[year]
            new_f = self.tmpDir + self.extractFilename(f)
            shutil.copyfile(f, new_f)
            newPoolList1[year] = new_f
            
            f = self.bootstrapPoolProject2[year]
            new_f = self.tmpDir + self.extractFilename(f)
            shutil.copyfile(f, new_f)
            newPoolList2[year] = new_f
        self.bootstrapPoolProject1 = newPoolList1
        self.bootstrapPoolProject2 = newPoolList2        
        
        bootstrapResultM1 = dict()
        bootstrapResultM2 = dict()
        bootstrapObservations = dict()
        
        for year in self.decadals:
            try:
                tmp = bootstrapResultM1[year]
            except:
                yearToSelect = choice(self.decadals)

                
                bootstrapResultM1[year] = self.bootstrapPoolProject1[yearToSelect]
                bootstrapResultM2[year] = self.bootstrapPoolProject2[yearToSelect]
                bootstrapObservations[year] = self.obsRemapped[yearToSelect]
                
                for i in range(1,6):
                    if yearToSelect+i < max(self.decadals) and year+i <= max(self.decadals):
                        
                        if yearToSelect+i in self.decadals: 
                            bootstrapResultM1[year+i] = self.bootstrapPoolProject1[yearToSelect+i]
                            bootstrapResultM2[year+i] = self.bootstrapPoolProject2[yearToSelect+i]
                            bootstrapObservations[year+i] = self.obsRemapped[yearToSelect+i]
                            

        return [bootstrapResultM1, bootstrapResultM2, bootstrapObservations]
   
    def _calcSignificance(self, bootstrap_folders, output_folder, plot_folder, file_to_check):
        '''
        Caluclate the significance for the intput field (file_to_check)
        The field is also plotted with significance crosses
        
        :param bootstrap_folder: path of bootstraped data
        :param output_folder: path
        :param plot_folder: path
        :param file_to_check: fn
        '''
        fn = file_to_check
        file_to_check = self.extractFilename(file_to_check)
        b_array_list = list()
        for folder in bootstrap_folders:
            
            b_array_list.append(folder+file_to_check)
            
        
        significance = Significance(self.tmpDir, self.outputPlots)
        if self.fieldmean:
            (min_val, max_val) = significance.checkSignificanceFldmean(b_array_list, fn)
        else:
            check_value = 1 if 'std_ratio' in fn else 0
            (sig_lon, sig_lat) = significance.checkSignificance(b_array_list, fn, check_value=check_value)
        
            if self.zonalmean:
                m = Plotter.plotVerticalProfile(fn, -1, 1, colormap='RedBlu', lonlatbox=self.lonlatbox)
                Plotter.addCrossesXY(m, fn+'_significance_mask')            
            else:
                min_val = -1
                max_val = 1
                if 'std_ratio' in fn:
                    min_val = 0.5
                    max_val = 2.
                m = Plotter.plotField(fn, min_val, max_val, colormap='RedBlu', lonlatbox=self.lonlatbox)
                Plotter.addCrosses(m, sig_lon, sig_lat)
            Plotter.saveFig(plot_folder, fn.split(output_folder)[-1])
     
    def calcSignificance(self, bootstrap_folders, output_folder, plot_folder):
        '''
        Multiprocessing of "_calcSignificance". Start 1 Process for every result file
        
        :param bootstrap_folder: path of bootstraped data
        :param output_folder: path
        :param plot_folder: path
        '''
        files_to_check = self.findFiles.getAllFilesInSubfolders(output_folder)
        if self.maskMissingValues:
            files_to_check = [fn for fn in files_to_check if fn.find('masked') != -1]
        
        procs = len(files_to_check)
        poolArgs =  zip([self]*procs, [bootstrap_folders]*procs, [output_folder]*procs, 
                        [plot_folder]*procs, files_to_check, ['_calcSignificance'] * procs)
        #self.multiProcess(poolArgs)
        
        #SINGLE PROCESS FOR DEBUGGING
        for fn in files_to_check:
            self._calcSignificance(bootstrap_folders, output_folder, plot_folder, fn)
        
        
        if self.fieldmean:
            outputPlots = plot_folder+'accuracy/'
            self.makeFolder(outputPlots)
            #taylor = TaylorPlotMurCSS(negativeCorr=False)
            #taylor.constructPlot(resultList, outputPlots)
            print 'Plotting fldmean'
            
            flag1 = self.constructName(self.fileNameFlag, exp='1', startYear='1', endYear='1')
            flag1 = flag1[4:]
            flag2 = self.constructName(self.fileNameFlag, exp='2', startYear='1', endYear='1')
            flag2= flag2[4:]
            fnFlagVs = self.constructName(self.fileNameFlagVs, exp='', startYear='1', endYear='1')
            fnFlagVs = fnFlagVs[4:]
            plot_list = [('correlation','Anomaly Correlation',[0,1]),('msss','Mean Squared Error Skill Score',[-1,1])]
            Plotter.plotLeadtimeseriesSign(files_to_check,[flag1,flag2,fnFlagVs],plot_list)
            Plotter.saveFig(outputPlots, 'accuracy_leadtimeseries_all')
            
            plot_list = [('correlation','Anomaly Correlation',[0,1]),('msss','Mean Squared Error Skill Score',[-1,1])]
            Plotter.plotLeadtimeseriesSign(files_to_check,[flag1],plot_list)
            Plotter.saveFig(outputPlots, flag1+'accuracy_leadtimeseries_input1')
            
            plot_list = [('correlation','Anomaly Correlation',[0,1]),('msss','Mean Squared Error Skill Score',[-1,1])]
            Plotter.plotLeadtimeseriesSign(files_to_check,[flag2],plot_list)
            Plotter.saveFig(outputPlots, flag2+'accuracy_leadtimeseries_input2')
            
            plot_list = [('correlation','Anomaly Correlation',[0,1]),('msss','Mean Squared Error Skill Score',[-1,1])]
            Plotter.plotLeadtimeseriesSign(files_to_check,[fnFlagVs],plot_list)
            Plotter.saveFig(outputPlots, fnFlagVs+'accuracy_leadtimeseries_versus')
         
    
def main(config_dict, baseDir):
    '''
    Main fuinction for the whole bootstrap process
    
    1. Normal MSSS calculation
    2. Prepare bootstrap data
    3. Calculate bootstraps
    4. Calcualte Significance and plot files
    
    :param config_dict: dictionary with all params needed for Msss calcualtion
    :param baseDir: dir of class
    '''
    try: config_dict.pop('baseDir') 
    except: pass
    bootstrap_number = config_dict.pop('bootstrap_number')
    config_dict.pop('significance')
    Msss = MsssBootstrap(baseDir=baseDir,**config_dict)
    try:
        Msss.prepareInput()
        Msss.analyze()
        remappedObservations = Msss.obsRemapped
        
        print 'Prepare Bootstrap data'
        (bootstrapList1, bootstrapList2) = Msss.prepareBootstrap()
        
        bootstrap_folders = list()
        for i in range(1,bootstrap_number+1):
            print '###############################################'
            print 'Bootstrap number '+str(i)
            print '###############################################'
            Msss.obsRemapped = remappedObservations
            bootstrap_config_dict = config_dict.copy()
            outputFolder = '/'.join([config_dict['cache'], 'bootstrap', 'msss', 'number'+str(i)]) + '/'
            cacheFolder = outputFolder + 'cache/'
            
            bootstrap_config_dict['output']= outputFolder
            bootstrap_config_dict['output_plots']= outputFolder+'plots/'
            bootstrap_config_dict['cache']= cacheFolder
            bootstrap_config_dict['bootstrap']= True
            bootstrap_config_dict['obsRemapped'] = remappedObservations
            bootstrap = MsssBootstrap(baseDir=baseDir,**bootstrap_config_dict)
            (bootstrap.input1Remapped, bootstrap.input2Remapped, bootstrap.observationRemapped) = Msss.bootstrapGoddard(cacheFolder)
            
            bootstrap.analyze()
            bootstrap_folders.append(bootstrap.outputDir)
    
        print '###############################################'
        print 'Calculating Significance'
        print '###############################################'
        bootstrap.calcSignificance(bootstrap_folders, Msss.outputDir, Msss.outputPlots)
        print "Plots produced in %s" % Msss.outputPlots   
        return Msss 
    finally: 
        Msss.tmpDir = config_dict['cache']
        #Msss.deleteCache()
        
        
