"""
Created on 25.10.2013

:author: Sebastian Illing
sebastian.illing@met.fu-berlin.de

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

import os
import time

from nondeamonpool import MyPool
import murcss_config


def unwrap_self_f(arg, **kwarg):
    """
    helper function to use multiprocessing inside a Class
    """
    functionToCall = arg[-1]
    arg = arg[:-1]
    if hasattr(arg[0], functionToCall):
        method = getattr(arg[0], functionToCall)
        arg = arg[1:]
        return method(*arg, **kwarg)
    else:
        raise NameError, '%s is not defined' % functionToCall


class ToolAbstract(object):
    
    
    def __init__(self, *args, **kwargs):
        
        # create folders for all kwargs starting with "output" or "cache"
        for key, value in kwargs.iteritems():
            if key.find('output') == 0 or key.find('cache') == 0:
                if not os.path.isdir(value):
                    os.makedirs(value)
                    
    def extractFilename(self, path):
        
        return path.split('/')[-1]
            
    def getYearRange(self):
        """
        Returns a list with start and endyear
       
        :return list of tuples with start and endyear
        """
        result = list()
        splitList = self.leadtimes.split(',')
        maxValue = list()
        minValue = list()
        for yrange in splitList:
            try:
                result.append((int(yrange), int(yrange)))
                maxValue.append(int(yrange))
                minValue.append(int(yrange))
            except:
                tmp_range = yrange.split('-')
                result.append((int(tmp_range[0]), int(tmp_range[1])))
                maxValue.append(int(tmp_range[1]))
                minValue.append(int(tmp_range[0]))    
        self.maxLeadtime = max(maxValue)           
        self.minLeadtime = min(minValue) 
        return result
    
    def multiProcess(self, poolArgs):    
        
        pool = MyPool(processes=min([len(poolArgs), murcss_config.proc_count]))
        result = pool.map(unwrap_self_f, poolArgs)
        pool.terminate()
        pool.close()
        return result
    
    def listToYeardict(self, liste, exp_list=None):
        if exp_list is None:
            exp_list = self.decadals
        
        resultDict = dict()
        for i in range(0, len(liste)):    
            resultDict[exp_list[i]] = liste[i]
        return resultDict
    
    def yeardictToList(self, yeardict):    
        tmpList = list()
        for year in self.decadals:
            tmpList.append(yeardict[year])
        return tmpList
    
    def deleteCache(self):
        """
        Delete all files in cacheDir
        """
        import shutil
        shutil.rmtree(self.tmpDir)
    
    def checkPath(self, filePath):
        """'
        """
        if(filePath[-1] != '/'):
            return filePath + '/'
        else:
            return filePath
        
    def constructName(self,structure, exp='', startYear='',endYear='',extra=''):
        output = list()
        
        for var in structure:
            if hasattr(self, var+exp):
                temp_val = getattr(self, var+exp)
                val = temp_val.replace('*','')
                output.append(val)
            elif hasattr(self, var):
                temp_val = getattr(self, var)
                val = temp_val.replace('*','')
                output.append(val)
            elif var == 'TIME':
                output.append(time.strftime('%Y%m%d-%H%M%S'))
            elif var == 'YEARRANGE':
                try:
                    output.append(str(self.decadals[0])+'-'+str(self.decadals[-1]))
                except:
                    output.append(str(self.experiments[0])+'-'+str(self.experiments[-1]))
            elif var == 'LEVEL':
                if self.level is not None:
                    output.append(str(self.level))
            elif var == 'SELLONLAT':
                if self.lonlatbox is not None:
                    output.append(self.lonlatbox)                
            elif var == 'LEADTIMES':
                output.append(startYear)
                output.append(endYear)
            else:
                escaped_var = var.replace('*', '')
                output.append(escaped_var)                
        if extra != '':
            output.append(extra)   
        return '_'.join(output)    

    def makeFolder(self, path):
        
        if not os.path.isdir(path):
            os.makedirs(path)
        return self.checkPath(path)
    
    def escapeAttributes(self):
        
        att = [a for a in dir(self) if not a.startswith('__') and not callable(getattr(self, a))]
        for name in att:
            val = getattr(self, name)
            if isinstance(val, str):
                val = val.replace('*', '')
                setattr(self, name, val)

    def getRandomStr(self, size=6):
        import random
        import string
        chars = string.ascii_uppercase + string.digits
        return ''.join(random.choice(chars) for x in range(size))
    
    def getPoolArgs(self, count, *args):
        
        new_args = list()
        for arg in args:
            if type(arg) is list and len(arg) == count:
                new_args.append(arg)
            else:
                new_args.append([arg]*count)
        return zip(*new_args)
