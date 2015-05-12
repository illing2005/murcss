'''
Created on 12.03.2013

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
import multiprocessing
import abc
import os
from cdo import *
cdo = Cdo()
from string import lowercase, translate, maketrans
import shutil

#from evaluation_system.model.file import *
from evaluation_system.model.solr import SolrFindFiles

from tool_abstract import ToolAbstract, unwrap_self_f
from findFilesAbstract import FindFilesAbstract

class FileError(Exception): pass
class NoFilesFoundError(FileError): pass
class UnexpectedFileFormat(FileError): pass
class NotEnoughYearsInFile(FileError): pass
class WrongDrsStruct(FileError): pass
class LevelNotFound(FileError): pass
class EnsembleMemberError(FileError): pass

class FindFiles(FindFilesAbstract):
    '''
    Wrapper class to use solr_search with "python friendly" output --> lists or dicts
    '''
               
    def getFiles(self,year,fileType, model, variable, time_frequency='mon', product='*', ensemblemembers='*', institute='*', exp_prefix='d*', maxleadtime=10, minLeadtime=1):
        '''
        Method to get model files with solr_search.
        
        :param year: decadal starting year
        :param fileType: baseline1, cmip5, historical or...
        :param model: model name i.e. MPI-ESM-LR
        :param variable: CMOR variable
        :param time_frequency: monthly, yearly, daily and so on
        
        :return: list with all ensemblemembers members found
        '''   
        #TODO: BUGFIX for minLeadyear
        minLeadtime=1                     
        output = list() 
        decStr = exp_prefix+str(year)
        project = fileType.lower()    
        tmpList = list()
        for fn in SolrFindFiles.search(experiment=decStr, latest_version=True, product=product, institute=institute,
                                      variable=variable, time_frequency=time_frequency, model=model, project=project):
            if(str(fn).split('.')[-1] == 'nc'):
                tmpList.append(str(fn))
        try:
            test = tmpList[0]
        except:
            import time
            time.sleep(5) # delays for 5 seconds
            for fn in SolrFindFiles.search(experiment=decStr, latest_version=True, product=product, institute=institute,
                                      variable=variable, time_frequency=time_frequency, model=model, project=project):
                if(str(fn).split('.')[-1] == 'nc'):
                    tmpList.append(str(fn))
            try:
                test = tmpList[0]
            except:
                if exp_prefix.find('*') != -1:
                    raise NoFilesFoundError, "Couldn't find files for %s in %s %s %s experiment: %s" % (variable, fileType, model, product, year)
                #OK we can't find files, now try one last time using only the exp_prefix, i.e. "historical"
                decStr = exp_prefix
                for fn in SolrFindFiles.search(experiment=exp_prefix, latest_version=True, product=product, institute=institute,
                                      variable=variable, time_frequency=time_frequency, model=model, project=project):
                    if(str(fn).split('.')[-1] == 'nc'):
                        tmpList.append(str(fn))
                try:
                    test = tmpList[0]
                except:
                    #OK, there are no Files...
                    raise NoFilesFoundError, "Couldn't find files for %s in %s %s %s experiment: %s" % (variable, fileType, model, product, year)  
        

              
        #Check if we have time-splitted files
        time_values = SolrFindFiles.facets(facets='time', experiment=decStr, latest_version=True, product=product, institute=institute,
                                           variable=variable, time_frequency=time_frequency, model=model, project=project)
        if len(time_values['time'])>1:
            tmpList = self.mergeSplittedFiles(tmpList)        
              
        #select only wanted ensemblemembers
        if type(ensemblemembers) == list and ensemblemembers[0] != '*':
            ensList = list()
            for ens in ensemblemembers:
                onlyfiles =  [f for f in tmpList if f.lower().find(ens) != -1]
                if len(onlyfiles) > 0:
                    ensList.append(onlyfiles[0])
		else:
		    raise EnsembleMemberError, "Ensemble member %s not found for  %s %s %s for starting year %s" % (ens,fileType, model, product, year)
            tmpList = ensList
        
        for fn in tmpList:
            years = cdo.showyear(input=str(fn))[0]
            yearList = years.split(' ')
            #print years    
            #print fn       
            if str(year+minLeadtime) not in yearList or str(year+maxleadtime) not in yearList:
                raise NotEnoughYearsInFile, "1Not enough years in %s %s %s for starting year %s" % (fileType, model, product, year)
            
            #if(len(years.split(' ')) > maxleadtime):
            selStr = ','.join(map(str,range(year+minLeadtime,year+1+maxleadtime)))
            fileName = str(fn).split('/')[-1]
            output.append(cdo.selyear(selStr, input=str(fn), output=self.tmpDir+fileName+self.getRandomStr()+'_'+str(year+minLeadtime)+'-'+str(year+maxleadtime), options='-f nc'))
            #else:    
            #    output.append(str(fn))
                
            if len(cdo.showyear(input=output[-1])[0].split(' ')) < maxleadtime-minLeadtime: 
                raise NotEnoughYearsInFile, "2Not enough years in %s %s %s for starting year %s" % (fileType, model, product, year)
                
        if(not output or not isinstance(output, list)):
            raise NoFilesFoundError, "Couldn't find files for %s in %s %s %s for starting year %s" % (variable, fileType, model, product, year)

        #check for curvilinear grid
        if(not hasattr(self,'curvilinearGrid') or self.curvilinearGrid == True):
            output = self.checkGrid(output, model)

        #user wants to select levels
        if self.level is not None:
            return self.selectLevel(output)
        else:
            return output
    
    def getReanalysis(self,year,fileType, experiment, variable, filePath='', time_frequency='mon', maxLeadtime=10, observation_ensemble='*', minLeadtime=1):
        '''
        Wrapper method to find reanalysis file with solr_search.
        
        :param year: startyear
        :param fileType: reanalysis or observation
        :param experiment: i.e. NCEP, HadCrut or MERRA
        :param variable: CMOR Variable
        :param time_frequency: monthly, yearly, daily and so on
        :return: "decadal" file with observations  
        '''
        #TODO: BUGFIX for minLeadyear
        minLeadtime=1
        reanFiles = list()
        if((experiment == 'HadCrut') and (variable == 'tas')):
            return self.getObsFiles(variable, year, maxLeadtime=maxLeadtime)
        
        
        
        #to use your own reanalysis data
        if os.path.isfile(self.observation):
            return self.getObsFiles(variable, year, maxLeadtime=maxLeadtime, minLeadtime=minLeadtime)
        
        
        
        if(not hasattr(self,'mergedReanFile')):
            #Observation or reanalysis?
            facet = SolrFindFiles.facets(facets='data_type', experiment=experiment, variable=variable, 
                                         time_frequency=time_frequency)
	    try:
                if 'reanalysis' in facet['data_type']:
                    searchList = SolrFindFiles.search(data_type=['reanalysis','observations'], experiment=experiment, variable=variable, 
                                         time_frequency=time_frequency, ensemble=observation_ensemble)
                else:
                    searchList = SolrFindFiles.search(data_type=['reanalysis','observations'], experiment=experiment, variable=variable, 
                                         time_frequency=time_frequency, data_structure='grid')
            except IndexError:
                raise NoFilesFoundError, "Couldn't find files for %s in %s" % (variable, experiment)
	    for fn in searchList:
		yearTmp = cdo.showyear(input=str(fn))[0]      
                fname = str(fn).split('/')[-1]
                #reanFiles.append(cdo.yearmean(input=str(fn), output=self.tmpDir+fname+'_YEARMEAN'))
                reanFiles.append(str(fn))
                #print reanFiles
                #if more than one year in File we break the loop and expect it to be a observationsfile
                if(len(yearTmp.split(' ')) > 1 ):
                    break
            if(len(reanFiles) == 0):
                raise NoFilesFoundError, "Couldn't find files for %s in %s " % (variable, experiment)    
            mergedFile = cdo.mergetime(input=' '.join(reanFiles), output=self.tmpDir+'mergedREAN_YEARMEAN')
            tmpMean = cdo.timmean(input=mergedFile)
            self.mergedReanFile = cdo.sub(input=' '.join([mergedFile, tmpMean]), output=self.tmpDir+'reananomalies.nc')
            #self.mergedReanFile = cdo.detrend(input=self.tmpDir+'reananomalies.nc', output=self.tmpDir+'reananomalies.nc_notrend')
            #print self.mergedReanFile
            if self.level is not None:
                self.mergedReanFile = self._selectLevel(self.mergedReanFile)
            
            #print self.mergedReanFile
        
        if(not hasattr(self,'mergedReanFile')):
            raise NoFilesFoundError, "Couldn't find files for %s in %s" % (variable, experiment)
            
        years = cdo.showyear(input=self.mergedReanFile)[0]
        if((years.find(str(year+minLeadtime)) != -1) and (years.find(str(year+maxLeadtime)) != -1)):
            #create tmp decadal file
            fileStr = ','.join(map(str,range(year+minLeadtime,year+maxLeadtime+1)))
            tmp= cdo.selyear(fileStr, input=self.mergedReanFile, output=self.tmpDir+'reanalysis_'+experiment+str(year+1)+'-'+str(year+maxLeadtime)+'.nc',options='-f nc')
            return tmp
        else:
            raise NotEnoughYearsInFile, "%s-%s are not part of %s reanalysis" % (year+minLeadtime, year+maxLeadtime, experiment)            
     

    def getObsFiles(self, variable, year, maxLeadtime=10, minLeadtime=1):
        '''
        Get the observation files from an specified folder
        
        :param variable:
        :param year: start year of decadal
        :return tmp file with maxLeadtime years of observation 
        '''    
        if not os.path.isfile(self.observation):
            raise NoFilesFoundError, '%s does not exist.' % (self.observation)
        
        variable_file = cdo.showname(input=self.observation)[0]
        if variable != variable_file:
            print 'WARNING: Variable in observation file is not %s. \n Variable %s will be renamed.' % (variable, variable_file)
            self.observation = cdo.chvar(variable_file+','+variable, input=self.observation, output=self.tmpDir+self.getFilename(self.observation))
        
        years = cdo.showyear(input=self.observation)[0]
        if(years.find(str(year+minLeadtime)) != -1) and (years.find(str(year+maxLeadtime)) != -1):
            #create tmp decadal file
            fileStr = ','.join(map(str,range(year+minLeadtime,year+maxLeadtime+1)))
            tmpFile =  cdo.selyear(fileStr, input=self.observation, 
                                   output=self.tmpDir+self.getFilename(self.observation)+'_'+str(year+minLeadtime)+'-'+str(year+maxLeadtime),options='-f nc')
	    if self.level is not None:
                return self._selectLevel(tmpFile)
            else:
                return tmpFile    
        else:
            if years.find(str(year+minLeadtime)) == -1:
                raise FileError, 'Can\'t find data for year %s in observational data! \n%s' % (year+minLeadtime, self.observation)
            if years.find(str(year+maxLeadtime)) == -1:
                raise FileError, 'Can\'t find data for year %s in observational data! \n%s' % (year+maxLeadtime, self.observation)   
