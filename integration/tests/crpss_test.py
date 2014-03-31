'''
Created on 21.02.2014

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
import scipy.stats as stats

from metrics.crpss import Crpss, NotEnoughEnsemblemembersFound
from metrics.filehandler import FileHandler
from metrics.msssBootstrap import main
from metrics.findFiles import NoFilesFoundError, NotEnoughYearsInFile

if 1 != 1:
    from integration.crpss import Crpss, NotEnoughEnsemblemembersFound
    from integration.metrics.filehandler import FileHandler
    from integration.metrics.msssBootstrap import main
    from integration.metrics.findFiles import NoFilesFoundError, NotEnoughYearsInFile
    
    
class TestCrpss(unittest.TestCase):
    
    tmp_dir = '/tmp/crpss_test'
    
    def getClassBaseDir(self):
        """:returns: the absolute path to the module defining the class implementing this plug-in."""
        subclass_file = os.path.abspath(sys.modules[self.__module__].__file__)
        return subclass_file
    
    def crpssDict(self):
        return dict(output=self.tmp_dir+'/output', 
                    output_plots=self.tmp_dir+'/plots', 
                    decadals=range(1960,1995,5),
                    variable='tas', 
                    project='BASELINE1', 
                    product1='output', 
                    institute1='mpi-m', 
                    model = 'mpi-esm-lr',
                    experiment='decs4e',
                    leadtimes='1,2-9', 
                    observation = 'HadCrut',
                    ensemblemembers='*', 
                    maskMissingValues = True,
                    result_grid=None, 
                    cache=self.tmp_dir+'/cache', 
                    baseDir = self.base+'/..', 
                    timeFreq = 'mon',  
                    bootstrapSwitch=False,
                    bootstrap_number = 500,
                    level=None, 
                    lonlatbox=None)
    
    def tearDown(self):
        if os.path.isdir(self.tmp_dir):
            import shutil
            shutil.rmtree(self.tmp_dir)
            
    def setUp(self):
        
        self.base= '/'.join(self.getClassBaseDir().split('/')[:-1])
        self.test_path = self.base+'/../../test_files/crpss'
        
    
    def testTestfiles(self):
        
        def crps(var,H,O):
            x = (H-O)/var
            crps = -var * (1./np.sqrt(np.pi) - 2. * stats.norm.pdf(x) - x * (2. * stats.norm.cdf(x) - 1.))
            return crps
        
        def condBias(H,O):  
            H_ensmean = np.mean(H, axis=1)
            r = np.corrcoef(H_ensmean, O)[0,1]
            std_H = np.std(H_ensmean)
            std_O = np.std(O)            
            cond_bias = r * std_O/std_H
            return cond_bias
            
        from create_testfiles import createCrpssInput, getCrpssTimeseries
        mse_goal = 4
        ensspread_goal = 3
        ensemble_size = 10
        length = 20    
        #create timeseries
        (t,hind,obs) = getCrpssTimeseries(mse_goal, ensspread_goal, ensemble_size, length)    
        
        (hindcasts,observations) = createCrpssInput(1960, hind, obs, ensemble_size, length)
        
        crpssDict = self.crpssDict()
        crpssDict['decadals'] = hindcasts.keys()
        crpssDict['observation'] = 'HadCrut'
        crpssDict['leadtimes'] = '1'
        crpss = Crpss(**crpssDict)
        crpss.inputRemapped = hindcasts
        crpss.observationRemapped = observations
        crpss.analyze()
        
        factor = length/(length-2.)
        factor2 = ensemble_size/(ensemble_size-1.)

        ensspread = np.mean((factor2*np.var(hind,axis=1)))**0.5        
        hind_mean = np.mean(hind,axis=1)
        
        #calculte crossvalidated mean
        hind_mean_anom = hind_mean.copy()
        crossval = np.zeros((length))
        for i in range(length):
            tmp = np.delete(hind_mean,i)
            crossval[i] =  np.mean(tmp)
            hind_mean_anom[i] = hind_mean[i] - crossval[i]
        hind_mean = hind_mean_anom.copy()
        
        #remove conditional bias
        bias = condBias(hind, obs)
        for i in range(length):
            hind[i,:] = hind[i,:] - crossval[i]
            hind[i,:] = hind[i,:] + (bias-1)*hind_mean[i]
        hind_mean_corrected = np.mean(hind,axis=1)
        RMSE = (factor*np.mean((np.mean(hind,axis=1)-obs)**2))**0.5
        
        #calculate crpss_goal
        crps1 = np.zeros((length))
        crps2 = np.zeros((length))
        for i in range(length):
            crps1[i] = crps(ensspread, hind_mean_corrected[i], obs[i])
            crps2[i] = crps(RMSE,hind_mean_corrected[i],obs[i]) 
        crpss_goal = 1 - np.mean(crps1)/np.mean(crps2)
        
        #test spreadscore
        tf = crpss.outputDir+'1-1/BASELINE1_output_mpi-esm-lr_decs4e/crpss/1_1_tas_BASELINE1_output_mpi-esm-lr_decs4e_1960-1979_ensspread_vs_referror.nc'
        t1 = FileHandler.openNetCDFFile(tf,mode='var')
        self.assertAlmostEqual(t1[12,12],ensspread/RMSE,2)

        #test crpss ens vs ref
        tf = crpss.outputDir+'1-1/BASELINE1_output_mpi-esm-lr_decs4e/crpss/1_1_tas_BASELINE1_output_mpi-esm-lr_decs4e_1960-1979_ens-vs-ref_crpss.nc'
        t1 = FileHandler.openNetCDFFile(tf,mode='var')
        self.assertAlmostEqual(t1[12,12],crpss_goal,2)
        
    def testCrpss(self):
        crpss_dict = self.crpssDict()
        crpss = Crpss(**crpss_dict)
        crpss.prepareInput()
        crpss.analyze()
        
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
        for i,f in enumerate(res_files):
            
            t1 = FileHandler.openNetCDFFile(res_files[i])
            t2 = FileHandler.openNetCDFFile(test_files[i])
            
            np.testing.assert_array_almost_equal(t1['variable'], t2['variable'], 3)
            #elf.assertTrue((t1['variable']==t2['variable']).all())
            self.assertTrue((t1['lon']==t2['lon']).all())
            self.assertTrue((t1['lat']==t2['lat']).all())
        
    def testBootstrap(self):
        crpss_dict = self.crpssDict()
        crpss_dict['bootstrapSwitch'] = True
        crpss_dict['bootstrap_number'] = 5
        crpss = Crpss(**crpss_dict)
        crpss.prepareInput()
        crpss.analyze()
          
    def testEnsemblemembers(self):
         
        crpss_dict = self.crpssDict()
        crpss_dict['ensemblemembers'] = 'r2i1p1,r3i1p1'
        crpss = Crpss(**crpss_dict)
        crpss.prepareInput()
        crpss.analyze()
         
        for year in crpss_dict['decadals']:
            self.assertTrue(len(crpss.inputDict[year]) == 2, 'Too many ensemblemembers found!')
             
    def testMaxLeadYearHistorical(self):
        crpss_dict = self.crpssDict()
        crpss_dict['decadals'] = range(1990,2000,5)
        crpss_dict['project'] = 'cmip5'
        crpss_dict['product1'] = 'output1'
        crpss_dict['leadtimes'] = '1,2'
        crpss_dict['experiment'] = 'historical'
        crpss = Crpss(**crpss_dict)
        crpss.prepareInput()
        crpss.analyze()
     
    def testErrors(self):
        crpss_dict = self.crpssDict()                 
        crpss_dict['variable'] = 'pr'
        crpss = Crpss(**crpss_dict)
        self.assertRaises(NoFilesFoundError, crpss.prepareInput)
        crpss_dict['variable'] = 'tas'    
        crpss_dict['observation'] = 'merra'
        crpss = Crpss(**crpss_dict)
        self.assertRaises(NotEnoughYearsInFile, crpss.prepareInput)
        crpss_dict['observation'] = 'HadCrut'
        crpss_dict['project'] = 'cmip5'
        crpss_dict['experiment'] = 'hist*'
        crpss_dict['product1'] = 'output1'
        crpss = Crpss(**crpss_dict)
        self.assertRaises(NoFilesFoundError, crpss.prepareInput)
 
        crpss_dict = self.crpssDict()
        crpss_dict['ensemblemembers'] = 'r1i1p1'
        crpss = Crpss(**crpss_dict)
        crpss.prepareInput()
        self.assertRaises(NotEnoughEnsemblemembersFound, crpss.analyze)
        
if __name__ == "__main__":
    unittest.main()       