'''
Created on 29.10.2013

@author: Sebastian Illing
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
'''
import unittest
import os, time
from metrics.tool_abstract import ToolAbstract

if 1!=1:
    from integration.metrics.tool_abstract import ToolAbstract

def wait(dt,pid):
    time.sleep(dt*pid)
    print pid
    return pid

class TestMetricAbstract(unittest.TestCase):


    def setUp(self):
        self.experiments = ['1960','1961','1962']
        self.tool_abstract = ToolAbstract()
          
    def tearDown(self):
        pass
     
    def testFolderCreation(self):
        #testfolder to create through __init__
        folders = ['/tmp/test_output/', 'tmp/test_output2/', '/tmp/no_output/']
        tmp = ToolAbstract(output1=folders[0],cache=folders[1],no_output=folders[2])
        self.assertTrue(os.path.isdir(folders[0]))
        self.assertTrue(os.path.isdir(folders[1]))
        self.assertFalse(os.path.isdir(folders[2]))
        
        folders.pop(2)
        for folder in folders: 
            tmp.tmpDir = folder
            tmp.deleteCache()
            self.assertFalse(os.path.isdir(folder))
            
    def testExtractFileName(self):
        path = '/miklip/integration/data4miklip/model/baseline1/output/MPI-M/MPI-ESM-LR/decs4e1993/mon/atmos/tas/r8i1p1/tas_Amon_MPI-ESM-LR_decs4e1993_r8i1p1_199401-200312.nc'
        self.assertEqual(self.tool_abstract.extractFilename(path), 'tas_Amon_MPI-ESM-LR_decs4e1993_r8i1p1_199401-200312.nc')
        
    def testGetYearRange(self):
        year_range = [(1,1),(2,5),(2,9)]
        self.tool_abstract.leadtimes = '1,2-5,2-9'
        self.assertEqual(self.tool_abstract.getYearRange(), year_range)
    
    def testListToYeardict(self):
        liste = ['a', 'b', 'c']
        self.failUnlessRaises(AttributeError, self.tool_abstract.listToYeardict, liste)
        self.failUnlessRaises(IndexError, self.tool_abstract.listToYeardict, liste,['1960','1961'])    
        
        result = self.tool_abstract.listToYeardict(liste, self.experiments)
        result2be = {'1960': 'a', '1961': 'b', '1962': 'c'}
        self.assertEqual(result,result2be)
         
    def testCheckPath(self):
        
        fn = 'miklip/integration/data4miklip/model/baseline1/output/MPI-M/MPI-ESM-LR/decs4e1993/mon/atmos/tas/r10i1p1/'
        self.assertEqual(self.tool_abstract.checkPath(fn),fn)
        fn2 = fn[:-1]  
        self.assertNotEqual(self.tool_abstract.checkPath(fn2),fn2)
        self.assertEqual(self.tool_abstract.checkPath(fn2),fn)
        
    def testConstructName(self):
        structure = 'variable_fileType_product_model_YEARRANGE'.split('_')
        
        self.tool_abstract.experiments = self.experiments
        self.assertEqual(self.tool_abstract.constructName(structure),'variable_fileType_product_model_1960-1962')
        
        self.tool_abstract.variable = 'tas'
        self.tool_abstract.fileType = 'baseline1'
        self.tool_abstract.product = 'output'
        self.tool_abstract.model = 'mpi-esm-lr'
        self.assertEqual(self.tool_abstract.constructName(structure), 'tas_baseline1_output_mpi-esm-lr_1960-1962')
        
    def testMultiProcessing(self):
        
        self.tool_abstract.wait = wait #add wait method to tool_abstract
        poolArgs = zip([self.tool_abstract]*5, [0.1]*5, range(1,6), ['wait']*5)
        self.assertEqual(self.tool_abstract.multiProcess(poolArgs), [1,2,3,4,5])
            
    def testMakeFolder(self):
        
        tmp_path = '/tmp/tool_abstract'
        self.tool_abstract.makeFolder(tmp_path)
        self.assertTrue(os.path.isdir(tmp_path))
        self.tool_abstract.tmpDir = tmp_path
        self.tool_abstract.deleteCache()    
        self.assertFalse(os.path.isdir(tmp_path))
        
    def testGetRandomStr(self):
        
        tmp_random = ''
        for i in range(1,11):
            tmp_random2 = self.tool_abstract.getRandomStr()
            self.assertNotEqual(tmp_random, tmp_random2)
            tmp_random = tmp_random2
        
if __name__ == "__main__":
    unittest.main()