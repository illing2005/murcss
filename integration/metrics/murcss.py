#!/usr/bin/env python
# encoding: utf-8
'''
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


import sys
import getopt
import os

from metrics import __version__
from metrics.msss import Msss
from metrics.crpss import Crpss

_short_args = "hnd"
_args = ['help', 'config_file=']

def getMurcssDict():
    '''Standard input of Murcss'''
    return dict(output='/tmp/murcss/output', output_plots='/tmp/murcss/plots', decadals='1960,1965,1970,1975,1980,1985', 
            variable='tas', cache='/tmp/murcss/cache', baseDir='..', maskMissingValues=True, model1='mpi-esm-lr', 
            analysis_type='basic', output_type='map', model2='mpi-esm-lr', project1='miklip', project2='miklip', observation='HadCrut', product1='initialized', product2='uninitialized', 
            ensemblemembers1='*', ensemblemembers2='*', institute1='mpi-m', institute2='mpi-m', leadtimes='1,2-5', experiment1='initialized', experiment2='uninitialized',
            result_grid='r72x36', level=None, lonlatbox=None, fieldmean=None, significance=False, bootstrap_number=100, metrics='all')

def printHelp():
    print """MurCSS (%s): \nThe tool calculates the Mean Squared Error Skill Score (MSESS) its decomposition (Correlation + Conditional Bi
as) and the Continuous Ranked Probability Skill Score (CRPSS)
as proposed by Goddard et al. [2013]. 
The MSSS of both models and the MSSS "between" the two models (model versions) are calculated for different leadtimes.
The CRPSS is calculated for both models defined by the input parameters. 

Plots are stored in "output_plots" and reusable NetCDF files are saved in "output" 
    
Options:
output            (default: "/tmp/murcss/output/"
                  The Output directory
output_plots      (default: "/tmp/murcss/plots/)
                  Output directory of produced plots
metrics           (default: all) [mandatory]
                  Here you can specify which metrics you want to calculate.
                  "accuracy": MSESS, Correlation, Conditional Bias.
                  "ensemble_spread": CRPSS, ESS
analysis_type     (default: 'map')
		  You can choose the output type of the analysis. Valid
                  options are "map", "fieldmean", and "zonalmean
output_type       (default: 'basic')
		  Basic only the basic plots like correlation, msess,
                  conditional bias, less, and crpss are produced. If you
                  choose additional, basic and additional plots will be
                  produced. For more information see read the documentation
decadals          Specify the experiments you want to use. I.e.
                  1960,1965,1970,..,1995.
variable          (default: "tas") [mandatory]
                  The name of the variable you want to analyze (CMOR)
project1          (default: "miklip") [mandatory]
                  miklip, or whatever your projectname might be.
product1          (default: "initialized") [mandatory]
                  Product 1 to analyze
institute1        (default: "mpi-m") [mandatory]
                  I.e. mpi-m
model1            (default: "mpi-esm-lr") [mandatory]
                  Model 1 to analyze
experiment1       (default: "decs4e") [mandatory]
                  Prefix for experiments. I.e. "initialized" or "uninitialized"
ensemblemembers1  (default: "*")
                  Here you can specify which ensemble members you want to use.
                  Please insert a comma separated list. I.e. "r1i1p1,r2i1p1,..".
                  If you leave this blank all members are used!
project2          (default: "miklip") [mandatory]
                  miklip, or whatever your projectname might be.
product2          (default: "uninitialized") [mandatory]
                  Product 1 to analyze
institute2        (default: "mpi-m") [mandatory]
                  I.e. mpi-m
model2            (default: "mpi-esm-lr") [mandatory]
                  Model 2 to analyze
experiment2       (default: "uninitialized") [mandatory]
                  Prefix for experiments. I.e. "initialized" or "uninitialized"
ensemblemembers2  (default: "*")
                  Here you can specify which ensemble members you want to use.
                  Please insert a comma separated list. I.e. "r1i1p1,r2i1p1,..".
                  If you leave this blank all members are used!
leadtimes         (default: 1,2-9) [mandatory]
                  Leadtimes to analyze
observation       (default: None) [mandatory]
                  Specify an observation file.
significance      (default: False)
                  Whether you want to calculate significance levels.
                  WARNING: This could take up to 1 day!
bootstrap_number  (default: 100)
                  Number of bootstrap runs.
level             (default: None)
                  Level to select. If you are using 3D-Files
lonlatbox         (default: <undefined>)
                  Here you can specify a region. I.e. -100,40,20,80
maskMissingValues (default: True)
                  Whether you want to mask missing values in observation file or
                  not
cache             (default: "/tmp/murcss/cache/")
                  Workdir
result_grid       (default: <undefined>)
                  You can specify a gridfile or a grid description like r72x36
months            (default: <undefined>)
                  If you analyze "seasonal" files and the experiment name
                  contains a month""" % (__version__)

def main(argv=None): # IGNORE:C0111
    '''Command line options.'''
    if argv is None:
       argv = sys.argv[1:]
    args, lastargs = getopt.getopt(argv, _short_args, _args)
    
    for flag,arg in args:
        if flag=='-h' or flag=='--help':
            return printHelp()
    
    murcss_dict = getMurcssDict()
    for arg in lastargs:
        if '=' not in arg:
            raise CommandError("Invalid format for query: %s" % arg)
            
        items = arg.split('=')
        key, value = items[0], ''.join(items[1:])
        murcss_dict[key] = value
            
    murcss_dict['decadals'] = map(int,murcss_dict['decadals'].split(','))
    murcss_dict['bootstrap_number'] = int(murcss_dict['bootstrap_number'])
    
    analysis_type = murcss_dict.pop('analysis_type')
    if analysis_type == 'map':
        murcss_dict['fieldmean'] = False
        murcss_dict['zonalmean'] = False
    elif analysis_type == 'fieldmean':
        murcss_dict['fieldmean'] = True
        murcss_dict['zonalmean'] = False
    elif analysis_type == 'zonalmean':
        murcss_dict['fieldmean'] = False
        murcss_dict['zonalmean'] = True     

    output_type = murcss_dict.pop('output_type')
    if output_type == 'basic':
        murcss_dict['basic_output'] = True
    else:
        murcss_dict['basic_output'] = False

    #Check for right metics value
    metrics = murcss_dict.pop('metrics')
    if metrics not in ['all','accuracy','ensemble_spread']:
        raise Exception, '%s is not a valid option. Only valid options for metric are: "all", "ensemble_spread", and "accuracy".' % (metrics)
    
    if metrics in ['all','accuracy']:
        print '#######################'
        print 'Calculating the MSSS'
        print '#######################'
        
        #Calculation of MSESS
        if(not murcss_dict['significance']):
            msss_dict = murcss_dict.copy()
            msss_dict.pop('bootstrap_number')
            msss_dict.pop('significance')
            msss = Msss(**msss_dict)
            msss.prepareInput()
            msss.analyze()   
            msss.deleteCache()
        else:
            from metrics.msssBootstrap import main
            msss = main(murcss_dict.copy(),'..')
    
    if metrics in ['all','ensemble_spread']:
        print '#######################'
        print 'Calculating the CRPSS for Model1'
        print '#######################'
        crpss1 = Crpss(output=murcss_dict['output'], 
                           output_plots=murcss_dict['output_plots'], 
                           basic_output=murcss_dict['basic_output'],
                           decadals=murcss_dict['decadals'],
                           variable=murcss_dict['variable'], 
                           
                           project = murcss_dict['project1'],
                           product1=murcss_dict['product1'], 
                           institute1=murcss_dict['institute1'], 
                           model=murcss_dict['model1'], 
                           experiment=murcss_dict['experiment1'],
                           ensemblemembers=murcss_dict['ensemblemembers1'], 
                           
                           observation = murcss_dict['observation'], 
                           leadtimes=murcss_dict['leadtimes'],
                           result_grid=murcss_dict['result_grid'], 
                           
                           maskMissingValues=murcss_dict['maskMissingValues'], 
                           bootstrapSwitch=murcss_dict['significance'], 
                           bootstrap_number=murcss_dict['bootstrap_number'],
                           level=murcss_dict['level'], 
                           lonlatbox=murcss_dict['lonlatbox'], 
                           fieldmean=murcss_dict['fieldmean'],
                           zonalmean=murcss_dict['zonalmean'],
                           cache=murcss_dict['cache'], 
                           
                           baseDir = murcss_dict['baseDir'])
        if metrics != 'ensemble_spread':
            crpss1.outputDir = msss.outputDir
            crpss1.outputPlots = msss.outputPlots
        crpss1.prepareInput()
        crpss1.analyze()
        crpss1.deleteCache()
        print '#######################'
        print 'Calculating the CRPSS for Model2'
        print '#######################'    
        crpss2 = Crpss(output=murcss_dict['output'], 
                           output_plots=murcss_dict['output_plots'], 
                           basic_output=murcss_dict['basic_output'],
                           decadals=murcss_dict['decadals'],
                           variable=murcss_dict['variable'], 
                           
                           project = murcss_dict['project2'],
                           product1=murcss_dict['product2'], 
                           institute1=murcss_dict['institute2'], 
                           model=murcss_dict['model2'], 
                           experiment=murcss_dict['experiment2'],
                           ensemblemembers=murcss_dict['ensemblemembers2'], 
                           
                           observation = murcss_dict['observation'], 
                           leadtimes=murcss_dict['leadtimes'],
                           result_grid=murcss_dict['result_grid'], 
                           
                           maskMissingValues=murcss_dict['maskMissingValues'], 
                           bootstrapSwitch=murcss_dict['significance'], 
                           bootstrap_number=murcss_dict['bootstrap_number'],
                           level=murcss_dict['level'], 
                           lonlatbox=murcss_dict['lonlatbox'], 
                           fieldmean=murcss_dict['fieldmean'],
                           zonalmean=murcss_dict['zonalmean'],
                           cache=murcss_dict['cache'], 
                           input_part='input2',
                           baseDir = murcss_dict['baseDir'])
        if metrics != 'ensemble_spread':
            crpss2.outputDir = msss.outputDir
            crpss2.outputPlots = msss.outputPlots
        crpss2.prepareInput()
        crpss2.analyze()
        crpss2.deleteCache()    

    import metrics.msss, metrics.metricAbstract, metrics.filehandler, metrics.taylorplot, metrics.crpss, metrics.msssBootstrap
    del metrics.msss.cdo, metrics.metricAbstract.cdo, metrics.filehandler.cdo, metrics.taylorplot.cdo, metrics.crpss.cdo, metrics.msssBootstrap.cdo
    print 'Calculation finished.'
    print 'Plots produced in %s' %(msss.outputPlots,)

if __name__ == "__main__":
    sys.exit(main())
