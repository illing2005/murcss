'''
Created on 20.06.2014

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

#Max process number
proc_count = 16

#Specify which filesystem you are using. If you are using MurCSS inside of MiKlip set filesystem to 'miklip'
#Otherwise set it to 'custom'. Then the findFilesCustom instance will be used. 
file_system = 'custom' #miklip|custom

#If you are working on your own filesystem specify the folder structure:
DRS_STRUCTURE = {
         "root_dir":"/home/PUTYOUTPATHHERE/murcss/sample_data/",
         "parts_dir":"project/product/institute/model/experiment/time_frequency/realm/variable/ensemble/file_name".split('/'),
         "defaults" : {}
        }
