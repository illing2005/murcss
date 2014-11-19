'''
Created on 28.10.2014

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
from tool_abstract import ToolAbstract, unwrap_self_f
from cdo import *
cdo = Cdo()
from string import lowercase, translate, maketrans

class FileError(Exception): pass
class LevelNotFound(FileError): pass

class FindFilesAbstract(ToolAbstract):
    '''
    Abstract class containing the common methods of the FindFiles Types
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
        
        super(FindFilesAbstract,self).__init__(output_tmp=tmpDir, output_dir=output)
    
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

    def mergeSplittedFiles(self,fileList):
        '''
        If ensemblemember is not in one file we merge them together.
        '''
        fl_copy = list(fileList)
        result_list = list()
        while len(fl_copy) != 0:
            #get file trunk of first element
            fn = self.getFilename(fl_copy[0])
            fn = fn.split('_')
            fn = fn[:-1]
            #print fn
            fn = '_'.join(fn)
            #search for files with same file trunk
            merge_list = [f for f in fl_copy if f.find(fn) != -1]
            index_list = [i for i,f in enumerate(fl_copy) if f.find(fn) != -1]
            #print merge_list
            #merge found files
            if len(merge_list) > 1:
                #print 'merged'
                result_list.append(cdo.mergetime(input=' '.join(merge_list),
                                                 output=self.tmpDir+self.getFilename(merge_list[0]+self.getRandomStr()+'_mergedTimes')))
            else:
                result_list.append(merge_list[0])
            #delete already used files from list
            for offset, index in enumerate(index_list):
                index -= offset
                del fl_copy[index]    
        return result_list        