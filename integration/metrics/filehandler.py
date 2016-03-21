"""
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
"""

import shutil
from scipy.io import netcdf
from cdo import *
from numpy import float32
import numpy as np

cdo = Cdo()


class FileHandler(object):
    """
    Simple static class to read and write netcdf files using scipy.io.netcdf
    """
    @staticmethod
    def createDimension(f, value, name, var_type='double'):
        """
        Create a dimension in netcdf file f
        
        :param f: open netcdf file
        :param value: array or list with dimension data
        :param name: name of dimension
        :param var_type: type of dimension (double, float, ...)
        """
        f.createDimension(name, len(value))
        dim_var = f.createVariable(name, var_type, (name,))
        dim_var[:] = value
        return dim_var
    
    @staticmethod
    def _copyNetCDF(fn, new_fn, field):
        """
        Copies a given netcdf file using scipy.io
        Also stores a new field in the netcdf file
        
        :param fn: filename of netcdf to copy
        :param new_fn: new filename
        :param field: Array with new values for result file
        """
        shutil.copyfile(fn, new_fn)
        f_org = netcdf.netcdf_file(fn, 'r')
        try:
            lat_org = FileHandler.getVariableValue(f_org.variables['lat'])
            lon_org = FileHandler.getVariableValue(f_org.variables['lon'])
        except:
            lat_org = FileHandler.getVariableValue(f_org.variables['Y'])
            lon_org = FileHandler.getVariableValue(f_org.variables['X'])

        varname = cdo.showname(input=fn)[0]

        f = netcdf.netcdf_file(new_fn, 'w')
        dimensions = list()
        try:
            time_org = f_org.variables['time']
            if len(time_org.data) > 1:
                time = FileHandler.createDimension(f, time_org.data, 'time', 'double')
                time.units = time_org.units
                time.calendar = time_org.calendar
                dimensions.append('time')
        except:
            pass
        try:
            plev_org = f_org.variables['plev']
            if len(plev_org.data) > 1:
                plev = FileHandler.createDimension(f, plev_org.data, 'plev', 'double')
                plev.units = plev_org.units
                dimensions.append('plev')
        except:
            time_org = False

        if len(lat_org) > 1:
            FileHandler.createDimension(f, lat_org, 'lat', 'double')   
            dimensions.append('lat')
        if len(lon_org) > 1:
            FileHandler.createDimension(f, lon_org, 'lon', 'double')
            dimensions.append('lon')
        
        if dimensions == [] and field.shape == ():
            # THIS IS A HACK FOR FIELDMEAN
            if time_org:
                time = FileHandler.createDimension(f, time_org.data, 'time', 'double')
                time.units = time_org.units
                time.calendar = time_org.calendar
                dimensions.append('time') 
            FileHandler.createDimension(f, lat_org, 'lat', 'double')   
            dimensions.append('lat')
            FileHandler.createDimension(f, lon_org, 'lon', 'double')
            dimensions.append('lon')
            field = np.expand_dims(field, axis=1)
            field = np.expand_dims(field, axis=1)
            field = np.expand_dims(field, axis=1)
        
        variable = f.createVariable(varname, 'float', dimensions)     
        variable._FillValue = 1.e+20
        variable.missing_value = 1.e+20   
        variable[:] = float32(field)
        f.history = f_org.history
        f.close()
        return new_fn

    @staticmethod
    def saveToNetCDF(field, refFile, tag):
        """
        Saves a numpy array to an netCDF file. 
        
        :param field: numpy array
        :param refFile: File to "clone"
        :param tag: new filename
        :return filepath
        """
        fn = refFile + tag
        FileHandler._copyNetCDF(refFile,fn,field)
        return fn  
    
    @staticmethod
    def getVariableValue(var):
        """
        Get values of a variable using scipy.io
        
        :param var: netcdf.variable
        """
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
        """
        Open a file as numpy array
        Either get a dict with "var","lon","lat" or only "var"
        New Option mode='plev' gives also pressure levels
        """
        f = netcdf.netcdf_file(fn, 'r')
        try:
            lon = np.array(FileHandler.getVariableValue(f.variables['lon']))
            lat = np.array(FileHandler.getVariableValue(f.variables['lat']))
        except:
            try:
                lon = np.array(FileHandler.getVariableValue(f.variables['X']))
                lat = np.array(FileHandler.getVariableValue(f.variables['Y']))
            except:
                pass
        
        try:
            plev = np.array(FileHandler.getVariableValue(f.variables['plev']))
        except:
            pass
        
        varName = cdo.showname(input=fn)[0]     # TODO: Better way to get the variable name
        mVar = FileHandler.getVariableValue(f.variables[varName])

        mVar = mVar.squeeze()
        if len(np.shape(mVar)) == 3:
            mVar = mVar[0, :, :]
        elif len(np.shape(mVar)) == 2:
            mVar = mVar[:, :]

        f.close()
        if mode == 'all':   
            return {'variable': mVar, 'lon': lon, 'lat': lat}
        elif mode == 'plev':
            return {'variable': mVar, 'plev': plev, 'lat': lat, 'lon': lon}
        else:
            return mVar
