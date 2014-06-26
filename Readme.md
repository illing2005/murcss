MurCSS
=====
A Tool for Standardized Evaluation of Decadal Hindcast Systems

This Readme is for the MurCSS tool.

The main documentation can be found here /doc/build/index.html and here https://www-miklip.dkrz.de/about/murcss

Installation
-
Download and install MurCSS via pypi
```
pip install murcss
```

Requirements
-
MurCSS needs some additional python packages. If you install MurCSS using pypi, they will be downloaded and installed automatically. 
* Matplotlib >= 1.1.0
* Basemap
* NumPy >= 1.5.0
* SciPy >= 0.8.0
* CDOpy

Additionally MurCSS requires a working binary of the Climate Data Operators (CDO >= 1.5.4).

Usage
-
After installation you can use the tool typing `murcss` in your terminal.
* Help:
```
murcss --help 
```

Now you can see all parameters accepted by MurCSS.

* Outside of MiKlip:

If you are using MurCSS outside of MiKlip, you have to set the `file_system = 'custom'` and specify your data structure in `murcss_config.py`. 

```
#If you are working on your own filesystem specify the folder structure:
DRS_STRUCTURE = {
         "root_dir":"/scratch/ROOTDIR_OF_YOUR_DATA/",
         "parts_dir":"project/product/institute/model/experiment/time_frequency/realm/variable/ensemble/file_name".split('/'),
         "defaults" : {}
        }

```

* Generate Sample output:

Download the files in /sample_data/ and the files in /sample_output/. Adjust the DATA_STRUCTURE to your system. Now you can generate the sample output using
```
murcss variable=tas project1=baseline1 ...
```

Unittests
-
Download the files in /integration/tests/ and navigate to the directory.
Now you can run the tests using
```
python -m unittest discover . '*_test.py'
```
Currently some of the tests are designed for the MiKlip file system. Therfore they are commented out.

Support, Issues, Bugs
-
Please open an issue on GitHub or write an email to sebastian.illing@met.fu-berlin.de



License
-
Copyright (C) 2014 Sebastian Illing This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/
