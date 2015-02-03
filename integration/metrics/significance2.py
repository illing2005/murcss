'''
Created on 17.08.2013

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

import numpy as np
import os

import murcss_config
if murcss_config.file_system == 'miklip':
    from findFiles import FindFiles
else:
    from findFilesCustom import FindFilesCustom as FindFiles
from filehandler import FileHandler


class WrongArgument(Exception): pass

class Significance(object):
    '''
    Class for significance calculation
    '''
    def __init__(self, tmp_dir, output_dir):
        self.find_files = FindFiles(tmpDir = tmp_dir)
        
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        self.output_dir = output_dir
        
    def __getQuantile(self, mVar, q1,q2, precision=1.0):
            """
            Returns the q'th percentile of the distribution given in the argument
            'data'. Uses the 'precision' parameter to control the noise level.
            
            :param mvar: numpy variable
            :param q1,q2: lower and upper level
            :return: q'th percentile 
            """
            data = np.reshape(mVar, np.size(mVar))
            
            k, bins = np.histogram(data, bins=precision*np.sqrt(len(data)))
            norm_cumul = 1.0*k.cumsum() / len(data)
            return [bins[norm_cumul > q2][0], bins[norm_cumul > q1][0]]
    
    def __getLonLat(self, ifile):
        '''
        Get vectors with lon and lat values from a netdf file using cdo.griddes
        Was introduced because we were using a damaged grid
        lon,lat valued are located in the center of a gridbox
        
        :param ifile: netcdf fn
        :result: lon,lat vectors
        '''
        def searchGriddes(grid, needle):
            tmp = [x for x in grid if x.find(needle) != -1]
            return float(tmp[0].split(' ')[-1])
        
        from cdo import Cdo
        cdo = Cdo()  
        grid =cdo.griddes(input=ifile)
        try:
            xinc = searchGriddes(grid,'xinc')
            xsize = searchGriddes(grid,'xsize')
            xfirst = searchGriddes(grid,'xfirst') 
        except:
            xinc = searchGriddes(grid,'yinc')
            xsize = searchGriddes(grid,'ysize')
            xfirst = searchGriddes(grid,'yfirst') 
        yfirst = searchGriddes(grid,'yfirst')
        ysize = searchGriddes(grid,'ysize')
        yinc = searchGriddes(grid,'yinc') 
        lon = np.arange(xfirst+xinc/2, xsize*xinc+xfirst+xinc/2, xinc, dtype=float)
        lat = np.arange(yfirst+yinc/2, ysize*yinc+yfirst+yinc/2, yinc, dtype=float)
        lon = np.arange(xfirst, xsize*xinc+xfirst, xinc, dtype=float)
        lat = np.arange(yfirst, ysize*yinc+yfirst, yinc, dtype=float)
#        lon = np.arange(xfirst, xsize*xinc+xfirst, xinc, dtype=float)
#        lat = np.arange(yfirst, ysize*yinc+yfirst, yinc, dtype=float)
        return (lon,lat)
    
    def checkSignificanceFldmean(self, input, result_file, q1=0.05, q2=0.95, check_value=0):
        '''
        
        '''
        if type(input) == list:
            bootstrap_files = input
        elif os.path.isdir(input):
            bootstrap_files = self.find_files.getAllFilesInFolder(input)
        else:
            raise WrongArgument, 'Input has to be a list of files or a folder'
            
        bootstrap_arrays = list()
        
        #load all files to arrays
        for b_file in bootstrap_files:
            tmp_var = FileHandler.openNetCDFFile(b_file)
            bootstrap_arrays.append(tmp_var['variable'])
            
        result = FileHandler.openNetCDFFile(result_file, mode='var')

        test_sample = np.zeros(len(bootstrap_files))
        k=0
        for item in bootstrap_arrays:
            #if item[i,j] < 1e10 and item[i,j] > -1e10: #don't test missing values
            if item < 1e10 and item > -1e10:
                test_sample[k] = item  #construct testsample
                k+=1
        test_mean = np.mean(test_sample)
            
        quant = self.__getQuantile(test_sample, q1, q2)
        #mark significant values
        FileHandler.saveToNetCDF(quant[0], result_file, '_bootstrap_min_val')
        FileHandler.saveToNetCDF(quant[1], result_file, '_bootstrap_max_val')
        
        
        return (quant[0], quant[1])
        #self.save_significance_mask(sig_y, sig_x, result_file)
        #return (sig_lon, sig_lat)
        
    def checkSignificance(self, input, result_file, q1=0.05, q2=0.95, check_value=0):
        '''
        Checks if a value is statistically significant different from zero
        Uses the bootstrapped files
        
        :param input: list with files or folder of bootstraps
        :param result_file: file to check
        :return: sig_lon,sig_lat lists with significant points
        '''
        if type(input) == list:
            bootstrap_files = input
        elif os.path.isdir(input):
            bootstrap_files = self.find_files.getAllFilesInFolder(input)
        else:
            raise WrongArgument, 'Input has to be a list of files or a folder'
            
        bootstrap_arrays = list()
        
        #load all files to arrays
        for b_file in bootstrap_files:
            tmp_var = FileHandler.openNetCDFFile(b_file)
            bootstrap_arrays.append(tmp_var['variable'])
            
        result = FileHandler.openNetCDFFile(result_file, mode='var')
        (imax,jmax) = np.shape(result)    
        #print np.shape(result)    
        lon,lat = self.__getLonLat(result_file)
        sig_lon,sig_lat,sig_x,sig_y = list(),list(),list(),list()

        #loop array
        for i in range(0,imax):
            for j in range(0,jmax):
                test_sample = np.zeros(len(bootstrap_files))
                k=0
                for item in bootstrap_arrays:
                    if item[i,j] < 1e10 and item[i,j] > -1e10: #don't test missing values
                        test_sample[k] = item[i,j]  #construct testsample
                        k+=1
                test_mean = np.mean(test_sample)
                    
                quant = self.__getQuantile(test_sample, q1, q2)
                #mark significant values
                if (quant[0] > check_value and quant[1] > check_value) or (quant[0] < check_value and quant[1] < check_value):
                    if test_mean > -1e19 and test_mean < 1e19:
                        #no missingvalue in result
                        if result[i,j] < 1e10 and result[i,j] > -1e10:
                            sig_lon.append(lon[j])
                            sig_lat.append(lat[i])
                            sig_x.append(j)
                            sig_y.append(i)
       
        self.save_significance_mask(sig_y, sig_x, result_file)
        return (sig_lon, sig_lat)

    def save_significance_mask(self, sign_lon, sign_lat, result_file):
        '''
        Saves the significance mask in a netcdf file
        
        :param sign_lon,sign_lat: significant points
        :param result_file: File to overwrite
        '''
        result = FileHandler.openNetCDFFile(result_file, mode='var')
        res_shape = result.shape
        
        sign_array = np.zeros(res_shape)
        for i in range(0,len(sign_lon)):
            #pass
            sign_array[sign_lon[i], sign_lat[i]] = 1
        
        return FileHandler.saveToNetCDF(sign_array, result_file, '_significance_mask')     
 