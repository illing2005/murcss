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
#import Scientific.IO.NetCDF as S
from scipy.io import netcdf
from cdo import *
cdo = Cdo()
from numpy import float32
import numpy as np

class FileHandler(object):
    '''
    Simple static class to read and write netcdf files unsing scipy.io.netcdf
    
    '''
    @staticmethod
    def createDimension(f,value,name,var_type='double'):
        '''
        Create a dimension in netcdf file f
        
        :param f: open netcdf file
        :param value: array or list with dimenstion data
        :param name: name of dimension
        :param var_type: type of dimenstion (double, float, ...)
        '''
	f.createDimension(name,len(value))
        dim_var = f.createVariable(name,var_type,(name,))
        dim_var[:] = value
        return dim_var
    
    @staticmethod
    def _copyNetCDF(fn, new_fn,field):
        '''
        Copies a given netcdf file using scipy.io
        Also stores a new field in the netcdf file
        
        :param fn: filename of netcdf to copy
        :param new_fn: new filename
        :param field: Array with new values for result file
        '''
        shutil.copyfile(fn, new_fn)
        f_org = netcdf.netcdf_file(fn,'r')
        try:
            lat_org = FileHandler.getVariableValue(f_org.variables['lat'])#.data
            lon_org = FileHandler.getVariableValue(f_org.variables['lon'])#.data
        except:
            lat_org = FileHandler.getVariableValue(f_org.variables['Y'])#.data
            lon_org = FileHandler.getVariableValue(f_org.variables['X'])#.data
        #time_org = f_org.variables['time']

        varname = cdo.showname(input=fn)[0]

        f = netcdf.netcdf_file(new_fn,'w')
#        f.createDimension('lon',len(lon_org))
#        lon = f.createVariable('lon','double',('lon',))
#        lon[:] = lon_org

#        f.createDimension('lat',len(lat_org))
#        lat = f.createVariable('lat','double',('lat',))
#        lat[:] = lat_org
        dimensions = list()
        try:
            time_org = f_org.variables['time']
            #f.createDimension('time',len(time_org.data))
            #time = f.createVariable('time','double',('time',))
            #time[:] = time_org.data	
            if len(time_org.data) > 1:
                time = FileHandler.createDimension(f, time_org.data, 'time', 'double')
                time.units = time_org.units
                time.calendar = time_org.calendar
                dimensions.append('time')
            #variable = f.createVariable(varname,'float',('time','lat','lon',))
        except:
            pass
        try:
            plev_org = f_org.variables['plev']
            if len(plev_org.data) > 1:
                plev = FileHandler.createDimension(f, plev_org.data, 'plev', 'double')
                plev.units = plev_org.units
                dimensions.append('plev')
        except:
            pass

        if len(lat_org) > 1:
            FileHandler.createDimension(f, lat_org, 'lat', 'double')   
            dimensions.append('lat')
        if len(lon_org) > 1:
            FileHandler.createDimension(f, lon_org, 'lon', 'double')
            dimensions.append('lon')

        
        if dimensions == [] and field.shape == ():
            #THIS IS A HACK FOR FIELDMEAN 
            time = FileHandler.createDimension(f, time_org.data, 'time', 'double')
            time.units = time_org.units
            time.calendar = time_org.calendar
            dimensions.append('time') 
            FileHandler.createDimension(f, lat_org, 'lat', 'double')   
            dimensions.append('lat')
            FileHandler.createDimension(f, lon_org, 'lon', 'double')
            dimensions.append('lon')
            field = np.expand_dims(field,axis=1)
            field = np.expand_dims(field,axis=1)
            field = np.expand_dims(field,axis=1)
        
        variable = f.createVariable(varname, 'float', dimensions)     
        variable._FillValue = 1.e+20
        variable.missing_value = 1.e+20   
        variable[:] = float32(field)
        f.history = f_org.history
        f.close()
        return new_fn

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
        #print refFile
        
        #shutil.copyfile(refFile, fn)
        FileHandler._copyNetCDF(refFile,fn,field)

        #f = S.NetCDFFile(fn, 'a'e)
        #f_org = netcdf.netcdf_file(refFile,'r')
        #f = netcdf.netcdf_file(fn,'w')
        #variable = f.variables[cdo.showname(input=fn)[0]]
        #variable[:] = float32(field)
        #f.close()    
        return fn  
    
    @staticmethod
    def getVariableValue(var):
        '''
        Get values of a variable using scipy.io
        
        :param var: netcdf.variable
        '''
        try:
            var = np.array(var.data)
        except:
            try:
                var = np.array(var.getValue())
            except:
                var = np.array(var.__array_data__)
        return var
    
    @staticmethod
    def openNetCDFFile(fn, mode='all'):
        '''
        Open a file as numpy array
        Either get a dict with "var","lon","lat" or only "var"
        New Option mode='plev' gives also pressure levels
        '''
        #f = S.NetCDFFile(fn, mode='r') 
        f = netcdf.netcdf_file(fn,'r')
        try:
            lon = np.array(FileHandler.getVariableValue(f.variables['lon']))#.getValue()
            lat = np.array(FileHandler.getVariableValue(f.variables['lat']))#.getValue()
        except:
            try:
                lon = np.array(FileHandler.getVariableValue(f.variables['X']))#.getValue()
                lat = np.array(FileHandler.getVariableValue(f.variables['Y']))#.getValue() 
            except:
                pass#print 'Can\'t find lon/lat variables.'
        
        try:
            plev = np.array(FileHandler.getVariableValue(f.variables['plev']))
        except:
            pass
        
        varName = cdo.showname(input = fn)[0]     #TODO: Better way to get the variable name
        mVar = FileHandler.getVariableValue(f.variables[varName])#.getValue()
#        try:
#            mVar = np.array(mVar.data)
#        except:
#            try:
#                mVar = np.array(mVar.getValue())
#            except:
#                mVar = np.array(mVar.__array_data__)
        #print mVar
        #print mVar 
        #print mVar.shape
        mVar = mVar.squeeze()
        if(len(np.shape(mVar)) == 3):
            mVar = mVar[0,:,:]  
        elif(len(np.shape(mVar)) == 2):
            mVar = mVar[:,:]
        #print mVar.shape
        f.close()
        if mode == 'all':   
            return {'variable': mVar, 'lon':lon, 'lat':lat}
        elif mode == 'plev':
            return {'variable': mVar, 'plev': plev, 'lat':lat,'lon':lon}
        else:
            #print mVar
            return mVar
