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
from metrics.findFiles import *

if 1 != 1:
    from integration.metrics.findFiles import *


class TestFindFiles(unittest.TestCase):


    def setUp(self):
        self.findFiles = FindFiles('/tmp/output_test/')
#        self.ensList = ['/miklip/integration/data4miklip/model/baseline1/output/MPI-M/MPI-ESM-LR/decs4e2000/mon/atmos/tas/r2i1p1/tas_Amon_MPI-ESM-LR_decs4e2000_r2i1p1_200101-201012.nc', '/miklip/integration/data4miklip/model/baseline1/output/MPI-M/MPI-ESM-LR/decs4e2000/mon/atmos/tas/r3i1p1/tas_Amon_MPI-ESM-LR_decs4e2000_r3i1p1_200101-201012.nc', '/miklip/integration/data4miklip/model/baseline1/output/MPI-M/MPI-ESM-LR/decs4e2000/mon/atmos/tas/r1i1p1/tas_Amon_MPI-ESM-LR_decs4e2000_r1i1p1_200101-201012.nc', '/miklip/integration/data4miklip/model/baseline1/output/MPI-M/MPI-ESM-LR/decs4e2000/mon/atmos/tas/r10i1p1/tas_Amon_MPI-ESM-LR_decs4e2000_r10i1p1_200101-201012.nc', '/miklip/integration/data4miklip/model/baseline1/output/MPI-M/MPI-ESM-LR/decs4e2000/mon/atmos/tas/r4i1p1/tas_Amon_MPI-ESM-LR_decs4e2000_r4i1p1_200101-201012.nc', '/miklip/integration/data4miklip/model/baseline1/output/MPI-M/MPI-ESM-LR/decs4e2000/mon/atmos/tas/r5i1p1/tas_Amon_MPI-ESM-LR_decs4e2000_r5i1p1_200101-201012.nc', '/miklip/integration/data4miklip/model/baseline1/output/MPI-M/MPI-ESM-LR/decs4e2000/mon/atmos/tas/r6i1p1/tas_Amon_MPI-ESM-LR_decs4e2000_r6i1p1_200101-201012.nc', '/miklip/integration/data4miklip/model/baseline1/output/MPI-M/MPI-ESM-LR/decs4e2000/mon/atmos/tas/r7i1p1/tas_Amon_MPI-ESM-LR_decs4e2000_r7i1p1_200101-201012.nc', '/miklip/integration/data4miklip/model/baseline1/output/MPI-M/MPI-ESM-LR/decs4e2000/mon/atmos/tas/r8i1p1/tas_Amon_MPI-ESM-LR_decs4e2000_r8i1p1_200101-201012.nc', '/miklip/integration/data4miklip/model/baseline1/output/MPI-M/MPI-ESM-LR/decs4e2000/mon/atmos/tas/r9i1p1/tas_Amon_MPI-ESM-LR_decs4e2000_r9i1p1_200101-201012.nc']
#        self.histList = ['./tmp/tas_Amon_MPI-ESM-LR_historical_r1i1p1_185001-200512.nc_1991-2000', './tmp/tas_Amon_MPI-ESM-LR_historical_r2i1p1_185001-200512.nc_1991-2000', './tmp/tas_Amon_MPI-ESM-LR_historical_r3i1p1_185001-200512.nc_1991-2000']
            
    def testGetFiles(self):
        
        #BASELINE 1 Tests
        listToBe = ['/miklip/integration/data4miklip/model/global/miklip/baseline1/output/MPI-M/MPI-ESM-LR/decs4e2000/mon/atmos/tas/r9i1p1/tas_Amon_MPI-ESM-LR_decs4e2000_r9i1p1_200101-201012.nc', '/miklip/integration/data4miklip/model/global/miklip/baseline1/output/MPI-M/MPI-ESM-LR/decs4e2000/mon/atmos/tas/r8i1p1/tas_Amon_MPI-ESM-LR_decs4e2000_r8i1p1_200101-201012.nc', '/miklip/integration/data4miklip/model/global/miklip/baseline1/output/MPI-M/MPI-ESM-LR/decs4e2000/mon/atmos/tas/r7i1p1/tas_Amon_MPI-ESM-LR_decs4e2000_r7i1p1_200101-201012.nc', '/miklip/integration/data4miklip/model/global/miklip/baseline1/output/MPI-M/MPI-ESM-LR/decs4e2000/mon/atmos/tas/r6i1p1/tas_Amon_MPI-ESM-LR_decs4e2000_r6i1p1_200101-201012.nc', '/miklip/integration/data4miklip/model/global/miklip/baseline1/output/MPI-M/MPI-ESM-LR/decs4e2000/mon/atmos/tas/r5i1p1/tas_Amon_MPI-ESM-LR_decs4e2000_r5i1p1_200101-201012.nc', '/miklip/integration/data4miklip/model/global/miklip/baseline1/output/MPI-M/MPI-ESM-LR/decs4e2000/mon/atmos/tas/r4i1p1/tas_Amon_MPI-ESM-LR_decs4e2000_r4i1p1_200101-201012.nc', '/miklip/integration/data4miklip/model/global/miklip/baseline1/output/MPI-M/MPI-ESM-LR/decs4e2000/mon/atmos/tas/r3i1p1/tas_Amon_MPI-ESM-LR_decs4e2000_r3i1p1_200101-201012.nc', '/miklip/integration/data4miklip/model/global/miklip/baseline1/output/MPI-M/MPI-ESM-LR/decs4e2000/mon/atmos/tas/r2i1p1/tas_Amon_MPI-ESM-LR_decs4e2000_r2i1p1_200101-201012.nc', '/miklip/integration/data4miklip/model/global/miklip/baseline1/output/MPI-M/MPI-ESM-LR/decs4e2000/mon/atmos/tas/r1i1p1/tas_Amon_MPI-ESM-LR_decs4e2000_r1i1p1_200101-201012.nc', '/miklip/integration/data4miklip/model/global/miklip/baseline1/output/MPI-M/MPI-ESM-LR/decs4e2000/mon/atmos/tas/r10i1p1/tas_Amon_MPI-ESM-LR_decs4e2000_r10i1p1_200101-201012.nc']
        result = self.findFiles.getFiles(2000,'BASELINE1', 'MPI-ESM-LR', 'tas', 'mon', 'output', ensemblemembers='*', institute='mpi-m')
        #print result
        self.assertFalse(not result, 'Variable is empty')                #result is empty
        self.assertTrue(isinstance(result, list), 'Variable is not a list object')   #is list
        self.assertTrue(sorted(listToBe) == sorted(result), 'Maybe Paths have changed')
        
        #BASELINE 0 Tests
        result = self.findFiles.getFiles(2000,'BASELINE0', 'MPI-ESM-LR', 'tas', 'mon', 'output1', ['r1i1p1','r2i1p1'], 'mpi-m')
        self.assertFalse(not result, 'Variable is empty')                #result is empty
        self.assertTrue(isinstance(result, list), 'Variable is not a list object')   #is list
        listToBe = ['/miklip/integration/data4miklip/model/global/miklip/baseline0/output1/MPI-M/MPI-ESM-LR/decadal2000/mon/atmos/Amon/r1i1p1/v20120529/tas/tas_Amon_MPI-ESM-LR_decadal2000_r1i1p1_200101-201012.nc', 
                    '/miklip/integration/data4miklip/model/global/miklip/baseline0/output1/MPI-M/MPI-ESM-LR/decadal2000/mon/atmos/Amon/r2i1p1/v20120529/tas/tas_Amon_MPI-ESM-LR_decadal2000_r2i1p1_200101-201012.nc']
        self.assertTrue(len(listToBe) == len(result), 'Ensemble count has changed')
        for item in result:
            self.assertTrue(isinstance(listToBe.index(item), int))
    
        #HISTORICAL Tests
        result = self.findFiles.getFiles(1990,'cmip5', 'MPI-ESM-LR', 'tas', 'mon', 'output1', '*', 'mpi-m', 'historical')
        self.assertFalse(not result, 'Variable is empty')                #result is empty
        histList = ['/tmp/output_test/tas_Amon_MPI-ESM-LR_historical_r3i1p1_185001-200512.nc_1991-2000', '/tmp/output_test/tas_Amon_MPI-ESM-LR_historical_r2i1p1_185001-200512.nc_1991-2000', '/tmp/output_test/tas_Amon_MPI-ESM-LR_historical_r1i1p1_185001-200512.nc_1991-2000', '/tmp/output_test/tas_Amon_MPI-ESM-LR_historical_r9i1p1_185001-200512.nc_1991-2000', '/tmp/output_test/tas_Amon_MPI-ESM-LR_historical_r8i1p1_185001-200512.nc_1991-2000', '/tmp/output_test/tas_Amon_MPI-ESM-LR_historical_r7i1p1_185001-200512.nc_1991-2000', '/tmp/output_test/tas_Amon_MPI-ESM-LR_historical_r6i1p1_185001-200512.nc_1991-2000', '/tmp/output_test/tas_Amon_MPI-ESM-LR_historical_r5i1p1_185001-200512.nc_1991-2000', '/tmp/output_test/tas_Amon_MPI-ESM-LR_historical_r4i1p1_185001-200512.nc_1991-2000', '/tmp/output_test/tas_Amon_MPI-ESM-LR_historical_r10i1p1_185001-200512.nc_1991-2000']#['/tmp/output_test/tas_Amon_MPI-ESM-LR_historical_r3i1p1_185001-200512.nc_1991-2000', '/tmp/output_test/tas_Amon_MPI-ESM-LR_historical_r2i1p1_185001-200512.nc_1991-2000', '/tmp/output_test/tas_Amon_MPI-ESM-LR_historical_r1i1p1_185001-200512.nc_1991-2000']
        self.assertTrue(isinstance(result, list), 'Variable is not a list object')   #is list
        self.assertTrue(len(histList) == len(result), 'Ensemble count has changed')
        for item in result:
            self.assertTrue(isinstance(histList.index(item), int))
        #test error if years are not in file
        self.assertRaises(NotEnoughYearsInFile, self.findFiles.getFiles, 2004,'cmip5', 'MPI-ESM-LR', 'tas', 'mon', 'output1', '*', 'mpi-m', 'historical')
        #test maxleadyear "feature"
        listToBe = ['/tmp/output_test/tas_Amon_MPI-ESM-LR_historical_r3i1p1_185001-200512.nc_2005-2005', '/tmp/output_test/tas_Amon_MPI-ESM-LR_historical_r2i1p1_185001-200512.nc_2005-2005', '/tmp/output_test/tas_Amon_MPI-ESM-LR_historical_r1i1p1_185001-200512.nc_2005-2005']
        result = self.findFiles.getFiles(2004,'cmip5', 'MPI-ESM-LR', 'tas', 'mon', 'output1', ['r1i1p1','r2i1p1','r3i1p1'], 'mpi-m', 'historical', 1)
        print result
        self.assertTrue(isinstance(result, list), 'Variable is not a list object')   #is list
        self.assertTrue(len(listToBe) == len(result), 'Ensemble count has changed')
        self.assertTrue(sorted(listToBe) == sorted(result), 'Maybe Paths have changed')
        
        #PROJECTDATA Tests
        result = self.findFiles.getFiles(2000, 'droughtclip', 'mpi-esm-lr', 'tas', 'mon', )
        self.assertFalse(not result, 'Variable is empty')                #result is empty
        listToBe = ['/miklip/integration/data4miklip/projectdata/DroughtClip/output/MPI-M/MPI-ESM-LR/dec26o2000/mon/atmos/tas/r1i1p1/tas_Amon_MPI-ESM-LR_dec26o2000_r1i1p1_200101-201012.nc', '/miklip/integration/data4miklip/projectdata/DroughtClip/output/MPI-M/MPI-ESM-LR/dec22o2000/mon/atmos/tas/r1i1p1/tas_Amon_MPI-ESM-LR_dec22o2000_r1i1p1_200101-201012.nc', '/miklip/integration/data4miklip/projectdata/DroughtClip/output/MPI-M/MPI-ESM-LR/dec08o2000/mon/atmos/tas/r1i1p1/tas_Amon_MPI-ESM-LR_dec08o2000_r1i1p1_200101-201012.nc']
        self.assertTrue(isinstance(result, list), 'Variable is not a list object')   #is list
        self.assertTrue(len(listToBe) == len(result), 'Ensemble count has changed')
        for item in result:
            self.assertTrue(isinstance(listToBe.index(item), int))
        self.assertTrue(sorted(listToBe)==sorted(result), 'Something in the files changed')
        
        #RAISE CHECK
        self.assertRaises(NoFilesFoundError, self.findFiles.getFiles, 2020, 'BASELINE0', 'MPI-ESM-LR', 'tas')
        
    def testGetReanalysis(self):
        
        #Raise error if no files found
        self.failUnlessRaises(NoFilesFoundError,self.findFiles.getReanalysis, 2000, 'reanalysis', 'merra2', 'tas')
        self.findFiles.deleteCache()
        self.findFiles = FindFiles('/tmp/output_test/')
        #Raise error if year is not in included in reanalysis
        self.failUnlessRaises(NotEnoughYearsInFile,self.findFiles.getReanalysis, 1960, 'reanalysis', 'merra', 'tas')
        self.findFiles.deleteCache()
        
        self.findFiles = FindFiles('/tmp/output_test/')
        fn = '/tmp/output_test/reanalysis_merra2001-2010.nc'
        self.assertEqual(fn,self.findFiles.getReanalysis( 2000, 'reanalysis', 'merra', 'tas'))
        
        self.assertTrue(hasattr(self.findFiles, 'mergedReanFile'))
        fn = '/tmp/output_test/reanalysis_merra2002-2011.nc'
        self.assertEqual(fn,self.findFiles.getReanalysis( 2001, 'reanalysis', 'merra', 'tas'))
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
