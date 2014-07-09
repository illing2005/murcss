'''
Created on 11.03.2013

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
from metrics.findFilesCustom import *
from metrics.murcss_config import *

if 1 != 1:
    from integration.metrics.findFilesCustom import *
    from integration.metrics.murcss_config import *

class TestFindFilesCustom(unittest.TestCase):


    def setUp(self):
        self.findFiles = FindFilesCustom('/tmp/output_test/')
        
    def testGetFiles(self):
        
        root_dir = DRS_STRUCTURE['root_dir']
        
        files_to_be = [root_dir+'miklip/initialized/mpi-m/mpi-esm-lr/initialized1960/mon/atmos/tas/r1i1p1/tas_Amon_mpi-esm-lr_initialized1960_r1i1p1_196101-196512.nc',
                       root_dir+'miklip/initialized/mpi-m/mpi-esm-lr/initialized1960/mon/atmos/tas/r2i1p1/tas_Amon_mpi-esm-lr_initialized1960_r2i1p1_196101-196512.nc']
        #test experiment 1960
        files = self.findFiles.getFiles(1960, 'miklip', 'mpi-esm-lr', 'tas', 'mon', 'initialized', '*', 'mpi-m', 'initialized', 5, 1)
        self.assertEqual(files_to_be, files, 'Wrong DRS_Structure. Please check murcss_config.py and your file system.')
        
        #test singele ensemble
        files = self.findFiles.getFiles(1960, 'miklip', 'mpi-esm-lr', 'tas', 'mon', 'initialized', 'r1i1p1', 'mpi-m', 'initialized', 5, 1)
        self.assertEqual([files_to_be[0]], files)
        
        #test uninitialized files
        files_to_be = ['/tmp/output_test/tas_Amon_mpi-esm-lr_uninitialized_r2i1p1_196001-200512.nc_1961-1965']
        files = self.findFiles.getFiles(1960, 'miklip', 'mpi-esm-lr', 'tas', 'mon', 'uninitialized', 'r2i1p1', 'mpi-m', 'uninitialized', 5, 1)
        self.assertEqual(files_to_be, files, 'Wrong DRS_Structure. Please check murcss_config.py and your file system.')
        
    def testErrors(self):
        
        self.assertRaises(NoFilesFoundError, self.findFiles.getFiles,1961, 'miklip', 'mpi-esm-lr', 'tas', 'mon', 'initialized', '*', 'mpi-m', 'initialized', 5, 1)
                
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
