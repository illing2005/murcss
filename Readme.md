MurCSS
=====
A Tool for Standardized Evaluation of Decadal Hindcast Systems

The tool calculates the Mean Squared Error Skill Score (MSESS) its
decomposition (Correlation + Conditional Bias) and the Continuous
Ranked Probability Skill Score (CRPSS) as proposed by Goddard et
al. [2013]. The MSESS of both models and the MSESS "between" the two
models (model versions) are calculated for different leadtimes. The
CRPSS is calculated for both models defined by the input parameters.

The main documentation can be found [here][local-docs] and [here][homepage].

Installation
-
Download and install MurCSS via [PyPI][]
```
pip install murcss
```

Or you can clone this repository:
```
git clone https://github.com/illing2005/murcss.git
```

Requirements
-
MurCSS needs some additional python packages. If you install MurCSS using [PyPI][], they will be downloaded and installed automatically. 
* Matplotlib >= 1.1.0
* Basemap
* NumPy >= 1.5.0
* SciPy >= 0.8.0
* CDOpy
* scientificpython

Additionally MurCSS requires a working binary of the Climate Data Operators (CDO >= 1.5.4).

Usage
-
After installation you can use the tool typing `murcss` in your terminal.
#### Help:
```
murcss --help 
```

Now you can see all parameters accepted by MurCSS.

#### Outside of MiKlip:

If you are using MurCSS outside of MiKlip, you have to set the `file_system = 'custom'` and specify your data structure in `murcss_config.py`. 

```
#If you are working on your own filesystem specify the folder structure:
DRS_STRUCTURE = {
         "root_dir":"/scratch/ROOTDIR_OF_YOUR_DATA/",
         "parts_dir":"project/product/institute/model/experiment/time_frequency/realm/variable/ensemble/file_name".split('/'),
         "defaults" : {}
        }

```

#### Generate Sample output:

To generate the sample files in `/sample_output/` download the files in `/sample_data/`. Adjust the `DRS_STRUCTURE` in `murcss_config.py` to your system. For the comparison with real observations you need to download e.g. a HadCRUT dataset from the Met Office Hadley Centre. [Webpage and user information][HadCRUT].
Download and unzip e.g. the median of the HadCRUT4 dataset 
```
cd /observations/path/
wget http://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/gridded_fields/HadCRUT.4.2.0.0.median_netcdf.zip .
unzip HadCRUT.4.2.0.0.median_netcdf.zip 
```
Now you can generate the sample output using
```
murcss variable=tas project1=miklip product1=initialized institute1=mpi-m model1=mpi-esm-lr experiment1=initialized project2=miklip product2=uninitialized institute2=mpi-m model2=mpi-esm-lr experiment2=uninitialized decadals=1960,1965,1970,1975,1980,1985,1990,1995,2000 leadtimes=1,2-5 result_grid=r72x36 observation=/PATH/TO/OBSERVATION_FILE metrics=all
```

####Inside a Python Script:

MurCSS consists of two types of metrics, accuracy and ensemble_spread. The first is located in `metrics.msss` and the second in `metrics.crpss`.
Here is a simple script to calculate the accuracy part:
```
   from metrics.msss import Msss 
   msss = Msss(**parameter_dict)   #Create an instance of Msss Class 
   msss.prepareInput()             #Search and prepare input Files (ie. remapping)
   msss.analyze()                  #Calculates the actual skill scores
   msss.deleteCache()              #Delete all temporary files 
```
First we create an instance of Msss. Msss takes a bunch of keywords arguments, which can be separated in two parts: 

* arguments to find the input data
* specific options  for calculation

For a complete list and description of the keyword arguments see `metrics.msss.Msss`. Then we search and prepare the input data for the final calculation. The prepared files will be stored in `msss.input1Remapped`, `msss.input2Remapped` and `msss.observationRemapped` for later use.
The main calculation is done in `msss.analyze()` (see :`metrics.msss.Msss.analyze`). After the calculation is finished we delete the working direcory. 

And here for the ensemble spread part:
```
   from metrics.crpss import Crpss
   crpss = Crpss(**parameter_dict) #Create an instance of Crpss Class
   crpss.prepareInput()            #Search and prepare input Files (ie. remapping)
   crpss.analyze()                 #Calculates the actual skill scores
   crpss.deleteCache()             #Delete all temporary files
```
`crpss.prepareInput()` stores the input files in `crpss.inputRemapped` and `crpss.observationRemapped`. As you can see the scripts looks almost the same. The major difference is that crpss expects only data information about one model (see `metrics.crpss.Crpss`).

Unittests
-
Download the files in `/integration/tests/` and navigate to the directory.
You should also download the sample_data and adapt the DRS_Structure to run the tests for the file input component. 
Now run the tests using
```
python -m unittest discover . '*_test.py'
```
Currently some of the tests are designed for the MiKlip file system. Therfore they are commented out.

Support, Issues, Bugs
-
Please open an issue on GitHub or write an email to `sebastian.illing@met.fu-berlin.de`.


License
-
Copyright (C) 2014 Sebastian Illing This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/

[local-docs]: ./doc/build/index.html
[homepage]: https://www-miklip.dkrz.de/about/murcss
[PyPI]: https://pypi.python.org/pypi/murcss
[HadCRUT]: http://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/download.html
