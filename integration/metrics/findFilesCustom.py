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
import glob
#from evaluation_system.model.file import *
#from evaluation_system.model.solr import SolrFindFiles

from tool_abstract import ToolAbstract, unwrap_self_f
import murcss_config

class FileError(Exception): pass
class NoFilesFoundError(FileError): pass
class UnexpectedFileFormat(FileError): pass
class NotEnoughYearsInFile(FileError): pass
class WrongDrsStruct(FileError): pass
class LevelNotFound(FileError): pass


class FindFilesCustom(ToolAbstract):
    '''
    Wrapper class to use solr_search with "python friendly" output --> lists or dicts
    '''
    def __init__(self, tmpDir = '/', observation='', level=None, output='/'):
        '''
        Constructor
        
        :param tmpDir: cache folder
        :param observation: folder of "special" observation data 
        '''
        self.tmpDir = self.checkPath(tmpDir)
        self.output = self.checkPath(output)
        self.observation = observation
        self.level = level
        
        super(FindFilesCustom,self).__init__(output_tmp=tmpDir, output_dir=output)
               
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
        bl = murcss_config.DRS_STRUCTURE
        search_dict = bl['defaults'].copy()
        if exp_prefix == 'uninitialized':
            experiment = 'uninitialized'
        else:
            experiment = exp_prefix+str(year)
        search_dict.update({'variable':variable,
                       'project':fileType,
                       'model':model,
                       'time_frequency':time_frequency,
                       'product':product,
                       'institute':institute,
                       'experiment':experiment,  
                       'ensemble':ensemblemembers,                     
                       })
        #search_dict.update(locals())

        local_path = bl['root_dir']
        for key in bl['parts_dir']:
            if key in search_dict:
                local_path = os.path.join(local_path, search_dict[key])
                del search_dict[key]    #remove it so we can see if all keys matched
            else:
                local_path = os.path.join(local_path, "*")

        #NOTE: We might have defaults for the datasets that are not appearing in the directories.
        if set(search_dict) - set(bl['defaults']):
            #ok, there are typos or non existing constraints in the search.
            #which are not in the defaults. Those are "strange" to the selected structure.
            mesg = "Unknown parameter(s) %s " % (','.join(search_dict))
            raise WrongDrsStruct, mesg

        files = list()
        for path in glob.iglob(local_path):
            files.append(path)
        
        if len(files)==0:
            raise NoFilesFoundError, "Couldn't find files for %s in %s %s %s for starting year %s" % (variable, fileType, model, product, year)
       
        output = list() 
        for fn in files:
            
            years = cdo.showyear(input=str(fn))[0]
            yearList = years.split(' ')    
            if str(year+minLeadtime) not in yearList or str(year+maxleadtime) not in yearList:
                raise NotEnoughYearsInFile, "Not enough years in %s %s %s for starting year %s" % (fileType, model, product, year)
            
            if(len(years.split(' ')) > maxleadtime):
                selStr = ','.join(map(str,range(year+minLeadtime,year+1+maxleadtime)))
                fileName = str(fn).split('/')[-1]
                output.append(cdo.selyear(selStr, input=str(fn), output=self.tmpDir+fileName+'_'+str(year+minLeadtime)+'-'+str(year+maxleadtime)))
            else:    
                output.append(str(fn))            
            
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
        #to use your own reanalysis data
        if os.path.isfile(self.observation):
            return self.getObsFiles(variable, year, maxLeadtime=maxLeadtime, minLeadtime=minLeadtime)
 
    def getObsFiles(self, variable, year, maxLeadtime=10, minLeadtime=1):
        '''
        Get the observation files from an specified folder
        
        :param variable:
        :param year: start year of decadal
        :return tmp file with maxLeadtime years of observation 
        '''    
        if not os.path.isfile(self.observation):
            raise NoFilesFoundError, '%s does not exist.' % (self.observation)
        
        years = cdo.showyear(input=self.observation)[0]
        if(years.find(str(year+minLeadtime)) != -1) and (years.find(str(year+maxLeadtime)) != -1):
            #create tmp decadal file
            fileStr = ','.join(map(str,range(year+minLeadtime,year+maxLeadtime+1)))
            tmpFile =  cdo.selyear(fileStr, input=self.observation, 
                                   output=self.tmpDir+self.getFilename(self.observation)+'_'+str(year+minLeadtime)+'-'+str(year+maxLeadtime))
            if self.level is not None:
                return self._selectLevel(tmpFile)
            else:
                return tmpFile    
        else:
            if years.find(str(year+minLeadtime)) == -1:
                raise FileError, 'Can\'t find data for year %s in observational data! \n%s' % (year+minLeadtime, self.observation)
            if years.find(str(year+maxLeadtime)) == -1:
                raise FileError, 'Can\'t find data for year %s in observational data! \n%s' % (year+maxLeadtime, self.observation)

    def getFilename(self, fn):
        '''
        Helper to extract a filename out of a path
        :deprecated !!!
        :param fn
        :return filename
        '''
        return self.extractFilename(fn)
    
    def getAllFilesInFolder(self, folder):
        from os import listdir
        from os.path import isfile, join
        onlyfiles = [ join(folder,f) for f in listdir(folder) if isfile(join(folder,f)) ]
        return onlyfiles
    
    def getAllFilesInSubfolders(self, folder):
        
        file_list = list()
        for path, subdirs, files in os.walk(folder):
            for name in files:
                file_list.append(os.path.join(path, name))
        return file_list        


    
        
