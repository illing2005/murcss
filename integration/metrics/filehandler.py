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

import shutil
import Scientific.IO.NetCDF as S
from cdo import *
cdo = Cdo()
from numpy import float32


class FileHandler(object):
    
    @staticmethod
    def saveToNetCDF(field, refFile, tag):
        '''
        Saves a numpy array to an netCDF file. 
        
        :param field: numpy array
        :param refFile: File to "clone"
        :param tag: new filename
        :return filepath
        '''
        fn = refFile + tag
        shutil.copyfile(refFile, fn)
        f = S.NetCDFFile(fn, 'a')
        variable = f.variables[cdo.showname(input=fn)[0]]
        variable[:] = float32(field)
        f.close()    
        return fn  
    
    @staticmethod
    def openNetCDFFile(fn, mode='all'):
        '''
        Open a file as numpy array
        Either get a dict with "var","lon","lat" or only "var"
        '''
        f = S.NetCDFFile(fn, mode='r')
        try:
            lon = f.variables['lon'].getValue()
            lat = f.variables['lat'].getValue()
        except:
            try:
                lon = f.variables['X'].getValue()
                lat = f.variables['Y'].getValue() 
            except:
                print 'Can\'t find lon/lat variables.'
                
        varName = cdo.showname(input = fn)[0]     #TODO: Better way to get the variable name
        mVar = f.variables[varName].getValue()
        
        mVar = mVar.squeeze()
        if(len(np.shape(mVar)) == 3):
            mVar = mVar[0,:,:]  
        elif(len(np.shape(mVar)) == 2):
            mVar = mVar[:,:]
        
        f.close()
        if mode == 'all':   
            return {'variable': mVar, 'lon':lon, 'lat':lat}
        else:
            return mVar