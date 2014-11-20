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
import os, time, sys
import numpy as np

from cdo import Cdo
cdo = Cdo()

from metrics.msss import Msss
from metrics.filehandler import FileHandler
from metrics.msssBootstrap import main
from metrics.findFiles import NoFilesFoundError, NotEnoughYearsInFile

if 1 != 1:
    from integration.metrics.msss import Msss
    from integration.metrics.filehandler import FileHandler
    from integration.metrics.msssBootstrap import main
    from integration.metrics.findFiles import NoFilesFoundError, NotEnoughYearsInFile
    
class TestMsss(unittest.TestCase):
    
    tmp_dir = '/tmp/msss_test'
    
    
    def getClassBaseDir(self):
        """:returns: the absolute path to the module defining the class implementing this plug-in."""
        subclass_file = os.path.abspath(sys.modules[self.__module__].__file__)
        return subclass_file#os.path.join(*self._splitPath(subclass_file)[:-len(self.__module__.split('.'))])
    
    def msssDict(self):
        return dict(output=self.tmp_dir+'/output', output_plots=self.tmp_dir+'/plots', decadals=range(1960,1995,5), 
                    variable='tas', cache=self.tmp_dir+'/cache', baseDir=self.base+'/..', maskMissingValues=True, model1='mpi-esm-lr', 
                    model2='mpi-esm-lr', project1='baseline1', project2='baseline0', observation='hadcrut3v', product1='output', product2='output1', 
                    ensemblemembers1='*', ensemblemembers2='*', institute1='mpi-m', institute2='mpi-m', leadtimes='1,2-9', experiment1='decs4e', experiment2='decadal',
                    result_grid=None, level=None, lonlatbox=None, fieldmean=None,basic_output=False)
    
    def tearDown(self):
        if os.path.isdir(self.tmp_dir):
            import shutil
            shutil.rmtree(self.tmp_dir)
    
    def setUp(self):
        
        self.base= '/'.join(self.getClassBaseDir().split('/')[:-1])
        self.test_path = self.base+'/../../test_files/msss'
     
    def testTestfiles(self):
        
        from create_testfiles import createMsssInput
        std1_goal=1. 
        std2_goal=2. 
        corr_goal=0.5
        std_ratio = std2_goal/std1_goal 
        cond_bias = corr_goal - std1_goal/std2_goal
        msss_goal = corr_goal**2 + cond_bias**2
        length = 50
        (hindcasts,observations) = createMsssInput(1960, std1_goal, std2_goal, corr_goal, length)
    
        msssDict = self.msssDict()
        msssDict['decadals'] = hindcasts.keys()
        msssDict['observation'] = 'hadcrut3v'
        msssDict['leadtimes'] = '1'
        msss = Msss(**msssDict)
        msss.input1Remapped = hindcasts
        msss.input2Remapped = hindcasts
        msss.observationRemapped = observations
        
        msss.analyze()
    
        goal_dict = {'msss_goal': '1-1/baseline1_output_mpi-esm-lr_decs4e_input1/accuracy/1_1_tas_baseline1_output_mpi-esm-lr_decs4e_1960-2009_msss.nc',
                     'corr_goal': '1-1/baseline1_output_mpi-esm-lr_decs4e_input1/accuracy/1_1_tas_baseline1_output_mpi-esm-lr_decs4e_1960-2009_correlation.nc',
                     'cond_bias': '1-1/baseline1_output_mpi-esm-lr_decs4e_input1/accuracy/1_1_tas_baseline1_output_mpi-esm-lr_decs4e_1960-2009_conditional_bias.nc',
                     'std_ratio': '1-1/baseline1_output_mpi-esm-lr_decs4e_input1/accuracy/1_1_tas_baseline1_output_mpi-esm-lr_decs4e_1960-2009_std_ratio.nc'}
        
        for key,val in goal_dict.iteritems():
            t1 = FileHandler.openNetCDFFile(msss.outputDir+val,mode='var')
            self.assertAlmostEqual(np.round(t1[12,12],2), locals()[key], 1)
            
    
    def testMsss(self):
        msssDict = self.msssDict()
        msss = Msss(**msssDict)
        msss.prepareInput()
        msss.analyze()
         
        #test file names
        res_files = [f for r,_,files in os.walk(self.tmp_dir+'/output') for f in files]
        test_files = [f for r,_,files in os.walk(self.test_path+'/output') for f in files]
        self.assertEqual(sorted(res_files),sorted(test_files))
        res_files = [f for r,_,files in os.walk(self.tmp_dir+'/plots') for f in files]
        test_files = [f for r,_,files in os.walk(self.test_path+'/plots') for f in files]
        self.assertEqual(sorted(res_files),sorted(test_files))
 
        #test values of all .nc files
        res_files = sorted([os.path.join(r,f) for r,_,files in os.walk(self.tmp_dir+'/output') for f in files])
        test_files = sorted([os.path.join(r,f) for r,_,files in os.walk(self.test_path+'/output') for f in files])
        #test_files = sorted([os.path.join(r,f) for r,_,files in os.walk(tets_p) for f in files])
#        print len(res_files)
#        print len(test_files)
        for i,f in enumerate(res_files):
            tmp_test = cdo.copy(input=test_files[i],options='-f nc') 
            t1 = FileHandler.openNetCDFFile(res_files[i])
            t2 = FileHandler.openNetCDFFile(tmp_test)#test_files[i])
#            print '_______________________________________-'
#            print res_files[i] 
#            print test_files[i]
	    np.testing.assert_array_almost_equal(t1['variable'], t2['variable'], 0)
            self.assertTrue((t1['lon']==t2['lon']).all())
            self.assertTrue((t1['lat']==t2['lat']).all())
             
    def testEnsemblemembers(self):
          
        msssDict = self.msssDict()
        msssDict['ensemblemembers1'] = 'r1i1p1'
        msssDict['ensemblemembers2'] = 'r2i1p1,r3i1p1'
        msss = Msss(**msssDict)
        msss.prepareInput()
        msss.analyze()
          
        for year in msssDict['decadals']:
            self.assertTrue(len(msss.ensList[year]) == 1, 'Too many ensemblemembers found!')
            self.assertTrue(len(msss.histList[year]) == 2, 'Too many ensemblemembers found!')
   
    def testSameModels(self):
          
        msssDict = self.msssDict()
        msssDict['project2'] = msssDict['project1']
        msssDict['product2'] = msssDict['product1']
        msssDict['experiment2'] = msssDict['experiment1']
        msssDict['experiment2'] = msssDict['experiment1']
        msss = Msss(**msssDict)
        msss.prepareInput()
        msss.analyze()
              
    def testMaxLeadYearHistorical(self):
        msssDict = self.msssDict()
        msssDict['decadals'] = range(1990,2000,5)
        msssDict['project2'] = 'cmip5'
        msssDict['leadtimes'] = '1,2'
        msssDict['experiment2'] = 'historical'
        msss = Msss(**msssDict)
        msss.prepareInput()
        msss.analyze()
     
    def testErrors(self):
        msss_dict = self.msssDict()
          
        msss_dict['variable'] = 'pr'
        msss = Msss(**msss_dict)
        self.assertRaises(NoFilesFoundError, msss.prepareInput)
        msss_dict['variable'] = 'tas'    
        msss_dict['observation'] = 'merra'
        msss = Msss(**msss_dict)
        self.assertRaises(NotEnoughYearsInFile, msss.prepareInput)
        msss_dict['observation'] = 'HadCrut'
        msss_dict['project2'] = 'cmip5'
        msss_dict['experiment2'] = 'hist*'
        msss = Msss(**msss_dict)
        self.assertRaises(NoFilesFoundError, msss.prepareInput)
             
    def testMsssBootstrap(self):    
        config_dict = self.msssDict()
        config_dict['significance'] = True
        config_dict['bootstrap_number'] = 2
        #run the bootstrapping
        main(config_dict, self.base+'/..')
        pass
    

if __name__ == "__main__":
    unittest.main()
