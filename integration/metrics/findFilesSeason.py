'''
Created on 16.03.2014

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
from datetime import datetime
from dateutil.relativedelta import relativedelta

from findFiles import *


class FindFilesSeason(FindFiles):
    '''
    Collects files for "Seasonal" evaluation
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
                print str(fn)
                if(str(fn).split('.')[-1] == 'nc'):
                    tmpList.append(str(fn))
            try:
                test = tmpList[0]
            except:
                if exp_prefix.find('*') != -1:
                    raise NoFilesFoundError, "Couldn't find files for %s in %s %s %s experiment: %s" % (variable, fileType, model, product, year)
                #OK we can't find files, now try one last time using only the exp_prefix, i.e. "historical"
                for fn in SolrFindFiles.search(experiment=exp_prefix, latest_version=True, product=product, institute=institute,
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
            #TODO: Throw exepvtion if date is not in file
            if len(str(year)) == 4:
                year = int(str(year)+'12')
            start_month = datetime.strftime(datetime.strptime(str(year),"%Y%m") + relativedelta(months=1),'%Y-%m-01')
            end_month = datetime.strftime(datetime.strptime(str(year),"%Y%m") + relativedelta(months=maxleadtime+1) - relativedelta(days=1),'%Y-%m-31')
            fileName = str(fn).split('/')[-1]
            output.append(cdo.seldate(','.join([start_month,end_month]), input=fn, 
                              output=self.tmpDir+fileName+str(year+1)+'-'+str(year+maxleadtime)+'.nc'))
                
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
        reanFiles = list()

        #to use your own reanalysis data
        if os.path.isfile(self.observation):
            return self.getObsFiles(variable, year, maxLeadtime=maxLeadtime)

        if(not hasattr(self,'mergedReanFile')):
            #Observation or reanalysis?
            facet = SolrFindFiles.facets(facets='data_type', experiment=experiment, variable=variable, 
                                         time_frequency=time_frequency)
            try:
                if facet['data_type'][0] == 'reanalysis':
                    searchList = SolrFindFiles.search(data_type=['reanalysis','observations'], experiment=experiment, variable=variable, 
                                         time_frequency=time_frequency)
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
                #if more than one year in File we break the loop and expect it to be a observationsfile
                if(len(yearTmp.split(' ')) > 1 ):
                    break
            if(len(reanFiles) == 0):
                raise NoFilesFoundError, "Couldn't find files for %s in %s" % (variable, experiment)    
            mergedFile = cdo.mergetime(input=' '.join(reanFiles), output=self.tmpDir+'mergedREAN_YEARMEAN')
            tmpMean = cdo.timmean(input=mergedFile)
            self.mergedReanFile = cdo.sub(input=' '.join([mergedFile, tmpMean]), output=self.tmpDir+'reananomalies.nc')
            if self.level is not None:
                self.mergedReanFile = self._selectLevel(self.mergedReanFile)
        
        if(not hasattr(self,'mergedReanFile')):
            raise NoFilesFoundError, "Couldn't find files for %s in %s" % (variable, experiment)
        
        if len(str(year)) == 4:
            year = int(str(year)+'12')    
        start_month = datetime.strftime(datetime.strptime(str(year),"%Y%m") + relativedelta(months=1),'%Y-%m-01')
        end_month = datetime.strftime(datetime.strptime(str(year),"%Y%m") + relativedelta(months=maxLeadtime+1) - relativedelta(days=1),'%Y-%m-31')
        
        tmp = cdo.seldate(','.join([start_month,end_month]), input=self.mergedReanFile, 
                          output=self.tmpDir+'reanalysis_'+experiment+str(year+1)+'-'+str(year+maxLeadtime)+'.nc')
        return tmp
       
     



    
        
