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

from evaluation_system.model.file import *
from evaluation_system.model.solr import SolrFindFiles

from tool_abstract import ToolAbstract, unwrap_self_f

class FileError(Exception): pass
class NoFilesFoundError(FileError): pass
class UnexpectedFileFormat(FileError): pass
class NotEnoughYearsInFile(FileError): pass
class WrongDrsStruct(FileError): pass
class LevelNotFound(FileError): pass


class FindFiles(ToolAbstract):
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
        
        super(FindFiles,self).__init__(output_tmp=tmpDir, output_dir=output)
               
    def getFiles(self,year,fileType, model, variable, time_frequency='mon', product='*', ensemblemembers='*', institute='*', exp_prefix='d*', maxleadtime=10):
        '''
        Method to get model files with solr_search.
        
        :param year: decadal starting year
        :param fileType: baseline1, cmip5, historical or...
        :param model: model name i.e. MPI-ESM-LR
        :param variable: CMOR variable
        :param time_frequency: monthly, yearly, daily and so on
        
        :return: list with all ensemblemembers members found
        '''                
        output = list() 
        decStr = exp_prefix+str(year)
        project = fileType.lower()    
        tmpList = list()
        for fn in DRSFile.solr_search(experiment=decStr, latest_version=True, product=product, institute=institute,
                                      variable=variable, time_frequency=time_frequency, model=model, project=project):
            if(str(fn).split('.')[-1] == 'nc'):
                tmpList.append(str(fn))
        try:
            test = tmpList[0]
        except:
            import time
            time.sleep(5) # delays for 5 seconds
            for fn in DRSFile.solr_search(experiment=decStr, latest_version=True, product=product, institute=institute,
                                      variable=variable, time_frequency=time_frequency, model=model, project=project):
                print str(fn)
                if(str(fn).split('.')[-1] == 'nc'):
                    tmpList.append(str(fn))
            try:
                test = tmpList[0]
            except:
                if exp_prefix.find('*') != -1:
                    raise NoFilesFoundError, "Couldn't find files for %s in %s %s %s experiment: %s" % (variable, fileType, model, product, year)
                #OK we can't find files, now try one last time using only the exp_prefix, i.e. "historical"
                for fn in DRSFile.solr_search(experiment=exp_prefix, latest_version=True, product=product, institute=institute,
                                      variable=variable, time_frequency=time_frequency, model=model, project=project):
                    if(str(fn).split('.')[-1] == 'nc'):
                        tmpList.append(str(fn))
                try:
                    test = tmpList[0]
                except:
                    #OK, there are no Files...
                    raise NoFilesFoundError, "Couldn't find files for %s in %s %s %s experiment: %s" % (variable, fileType, model, product, year)  
        
        #select only wanted ensemblemembers
        if type(ensemblemembers) == list and ensemblemembers[0] != '*':
            ensList = list()
            for ens in ensemblemembers:
                onlyfiles =  [f for f in tmpList if f.find(ens) != -1]
                if len(onlyfiles) > 0:
                    ensList.append(onlyfiles[0])

            tmpList = ensList
              
        for fn in tmpList:
            
            years = cdo.showyear(input=str(fn))[0]
            yearList = years.split(' ')
            #print years    
            #print fn       
            if str(year+1) not in yearList or str(year+maxleadtime) not in yearList:
                print year
                raise NotEnoughYearsInFile, "1Not enough years in %s %s %s for starting year %s" % (fileType, model, product, year)
            
            if(len(years.split(' ')) > maxleadtime):
                selStr = ','.join(map(str,range(year+1,year+1+maxleadtime)))
                fileName = str(fn).split('/')[-1]
                output.append(cdo.selyear(selStr, input=str(fn), output=self.tmpDir+fileName+'_'+str(year+1)+'-'+str(year+maxleadtime)))
            else:    
                output.append(str(fn))
                
            if len(cdo.showyear(input=output[-1])[0].split(' ')) < maxleadtime: 
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
    
    def getReanalysis(self,year,fileType, experiment, variable, filePath='', time_frequency='mon', maxLeadtime=10):
        '''
        Wrapper method to find reanalysis file with solr_search.
        
        :param year: startyear
        :param fileType: reanalysis or observation
        :param experiment: i.e. NCEP, HadCrut or MERRA
        :param variable: CMOR Variable
        :param time_frequency: monthly, yearly, daily and so on
        :return: "decadal" file with observations  
        '''
        reanFiles = list()
        if((experiment == 'HadCrut') and (variable == 'tas')):
            return self.getObsFiles(variable, year, maxLeadtime=maxLeadtime)
        
        
        
        #to use your own reanalysis data
        if os.path.isfile(self.observation):
            return self.getObsFiles(variable, year, maxLeadtime=maxLeadtime)
        
        
        
        if(not hasattr(self,'mergedReanFile')):
            #Observation or reanalysis?
            facet = SolrFindFiles.facets(facets='data_type', experiment=experiment, variable=variable, 
                                         time_frequency=time_frequency)
            try:
                if facet['data_type'][0] == 'reanalysis':
                    searchList = DRSFile.solr_search(data_type=['reanalysis','observations'], experiment=experiment, variable=variable, 
                                         time_frequency=time_frequency)
                else:
                    searchList = DRSFile.solr_search(data_type=['reanalysis','observations'], experiment=experiment, variable=variable, 
                                         time_frequency=time_frequency, data_structure='grid')
            except IndexError:
                raise NoFilesFoundError, "Couldn't find files for %s in %s" % (variable, experiment)
            
            for fn in searchList:
                yearTmp = cdo.showyear(input=str(fn))[0]      
                fname = str(fn).split('/')[-1]
                #reanFiles.append(cdo.yearmean(input=str(fn), output=self.tmpDir+fname+'_YEARMEAN'))
                reanFiles.append(str(fn))
                #if more than one year in File we break the loop and expect it to be a observationsfile
                if(len(yearTmp.split(' ')) > 1 ):
                    break
            if(len(reanFiles) == 0):
                raise NoFilesFoundError, "Couldn't find files for %s in %s" % (variable, experiment)    
            mergedFile = cdo.mergetime(input=' '.join(reanFiles), output=self.tmpDir+'mergedREAN_YEARMEAN')
            tmpMean = cdo.timmean(input=mergedFile)
            self.mergedReanFile = cdo.sub(input=' '.join([mergedFile, tmpMean]), output=self.tmpDir+'reananomalies.nc')
            #print self.mergedReanFile
            if self.level is not None:
                self.mergedReanFile = self._selectLevel(self.mergedReanFile)
            
            #print self.mergedReanFile
        
        if(not hasattr(self,'mergedReanFile')):
            raise NoFilesFoundError, "Couldn't find files for %s in %s" % (variable, experiment)
            
        years = cdo.showyear(input=self.mergedReanFile)[0]
        if((years.find(str(year+1)) != -1) and (years.find(str(year+maxLeadtime)) != -1)):
            #create tmp decadal file
            fileStr = ','.join(map(str,range(year+1,year+maxLeadtime+1)))
            tmp= cdo.selyear(fileStr, input=self.mergedReanFile, output=self.tmpDir+'reanalysis_'+experiment+str(year+1)+'-'+str(year+maxLeadtime)+'.nc')
            return tmp
        else:
            raise NotEnoughYearsInFile, "%s-%s are not part of %s reanalysis" % (year+1, year+maxLeadtime, experiment)            
     

    def getObsFiles(self, variable, year, maxLeadtime=10):
        '''
        Get the observation files from an specified folder
        
        :param variable:
        :param year: start year of decadal
        :return tmp file with maxLeadtime years of observation 
        '''    
        if not os.path.isfile(self.observation):
            raise NoFilesFoundError, '%s does not exist.' % (self.observation)
        
        years = cdo.showyear(input=self.observation)[0]
        if(years.find(str(year+1)) != -1) and (years.find(str(year+maxLeadtime)) != -1):
            #create tmp decadal file
            fileStr = ','.join(map(str,range(year+1,year+maxLeadtime+1)))
            tmpFile =  cdo.selyear(fileStr, input=self.observation, 
                                   output=self.tmpDir+self.getFilename(self.observation)+'_'+str(year+1)+'-'+str(year+maxLeadtime))
            if self.level is not None:
                return self._selectLevel(tmpFile)
            else:
                return tmpFile    
        else:
            if years.find(str(year+1)) == -1:
                raise FileError, 'Can\'t find data for year %s in observational data! \n%s' % (year+1, self.observation)
            if years.find(str(year+maxLeadtime)) == -1:
                raise FileError, 'Can\'t find data for year %s in observational data! \n%s' % (year+maxLeadtime, self.observation)
         
    def checkGrid(self,fList,model):
        '''
        Checks if the file has a curvlinear grid. And remaps to lonlat grid after 
        
        '''
        gridInfo = cdo.griddes(input=fList[0])
        gridType = gridInfo[3]
        if gridType.find('curvilinear') == -1:
            self.curvilinearGrid = False
            return fList
        else:
            self.curvilinearGrid = True
            lon = self.__str2int(gridInfo[11])
            lat = self.__str2int(gridInfo[12])
        #single process, becaus multiproccessing caused memory problems
        result = list()
        for fn in fList:
            result.append(self._ceckGrid(fn, model, lon, lat))
        return result
        
    def _ceckGrid(self, f, model, lon, lat):
        
        if model.find('MPI-ESM') != -1:
            lon=lon-1
            sel_str = '2,%s,1,%s' % (lon,lat)
            f = cdo.selindexbox(sel_str, input=f, output=self.tmpDir+self.getFilename(f)+'_sel_box')
                
        grid_str = 'r%sx%s' % (lon,lat)
        return cdo.remapbil(grid_str, input=f, output=self.tmpDir+self.getFilename(f)+'_lonlat')
        
    
    def getFilename(self, fn):
        '''
        Helper to extract a filename out of a path
        :deprecated !!!
        :param fn
        :return filename
        '''
        return self.extractFilename(fn)
    
    def __str2int(self, str):
        '''
        Filter digits and convert str to int
        
        :param str:
        :return int
        '''
        all = maketrans('', '')
        nodigs = all.translate(all, string.digits)
        return int(str.translate(all, nodigs))
    
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
        
    def _selectLevel(self, files):

        try:
            return cdo.sellevel(self.level, input=files, output=self.tmpDir+self.getFilename(files)+'_'+str(self.level)+'.nc')
        except:
            raise LevelNotFound, 'Level %s not found in %s' %(self.level, files)

            
    def selectLevel(self,fileList):
        '''
        Select a specific level from the files
        '''        
        #multi processing
        num_proc = len(fileList)
        pool = multiprocessing.Pool(processes=min([num_proc,24]))
        poolArgs = zip([self]*num_proc, fileList, ['_selectLevel']*num_proc)
        result =  pool.map(unwrap_self_f, poolArgs)
        pool.terminate()
        pool.close()
        
        return result

    
        
