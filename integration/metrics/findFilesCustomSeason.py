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
from findFilesAbstract import FindFilesAbstract
import murcss_config

class FileError(Exception): pass
class NoFilesFoundError(FileError): pass
class UnexpectedFileFormat(FileError): pass
class NotEnoughYearsInFile(FileError): pass
class WrongDrsStruct(FileError): pass
class LevelNotFound(FileError): pass


class FindFilesCustomSeason(FindFilesAbstract):
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
        #print local_path
        if len(files)==0:
            raise NoFilesFoundError, "Couldn't find files for %s in %s %s %s for starting year %s" % (variable, fileType, model, product, year)
       
        output = list() 
        for fn in files:
            if len(str(year)) == 4:
                year = int(str(year)+'12')
            start_month = datetime.strftime(datetime.strptime(str(year),"%Y%m") + relativedelta(months=1),'%Y-%m-01')
            end_month = datetime.strftime(datetime.strptime(str(year),"%Y%m") + relativedelta(months=maxleadtime+1) - relativedelta(days=1),'%Y-%m-31')
            fileName = str(fn).split('/')[-1]
            output.append(cdo.seldate(','.join([start_month,end_month]), input=fn, 
                              output=self.tmpDir+fileName+str(year+1)+'-'+str(year+maxleadtime)+'.nc', options='-f nc'))       

        #check for curvilinear grid
        if(not hasattr(self,'curvilinearGrid') or self.curvilinearGrid == True):
            output = self.checkGrid(output, model)

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
        
        variable_file = cdo.showname(input=self.observation)[0]
        if variable != variable_file:
            print 'WARNING: Variable in observation file is not %s. \n Variable will be renamed.' % (variable)
            self.observation = cdo.chvar(variable_file+','+variable, input=self.observation, output=self.tmpDir+self.getFilename(self.observation))
        
        if len(str(year)) == 4:
            year = int(str(year)+'12')    
        start_month = datetime.strftime(datetime.strptime(str(year),"%Y%m") + relativedelta(months=1),'%Y-%m-01')
        end_month = datetime.strftime(datetime.strptime(str(year),"%Y%m") + relativedelta(months=maxLeadtime+1) - relativedelta(days=1),'%Y-%m-31')
        
        tmp = cdo.seldate(','.join([start_month,end_month]), input=self.observation, 
                          output=self.tmpDir+'reanalysis_'+experiment+str(year+1)+'-'+str(year+maxLeadtime)+'.nc', options='-f nc')
        return tmp

