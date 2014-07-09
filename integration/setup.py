from setuptools import setup, find_packages # Always prefer setuptools over distutils
from codecs import open # To use a consistent encoding
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'Readme.md'), encoding='utf-8') as f:
    long_description = f.read()
    
import metrics

setup(
    name='murcss',

    version = metrics.__version__,

    description='MurCSS: A Tool for Standardized Evaluation of Decadal Hindcast Systems',
    long_description=long_description,

    url='https://github.com/illing2005/murcss',

    author='Sebastian Illing',
    author_email='sebastian.illing@met.fu-berlin.de',

    # Choose your license
    license='GPLv3',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        # 3 - Alpha
        # 4 - Beta
        # 5 - Production/Stable
        'Development Status :: 5 - Production/Stable',

#        'Programming Language :: Python :: 2',
#        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
    ],

    keywords='Decadal prediction, skill score, MiKlip, MSESS, CRPSS, model evaluation, verification',

    packages=find_packages(exclude=['contrib', 'docs', 'tests*']),

    install_requires=['matplotlib>=1.1.0',
                      'numpy>=1.5.0',
                      'cdo',
                      'scipy>=0.8.0',
                      'basemap'],

    # If there are data files included in your packages that need to be
    # installed, specify them here. If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
#    package_data={
#        'obs2': ['tas_obs_HadCrut3.nc'],
#    },

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages.
    # see http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
#    data_files=[('my_data', ['data/data_file'])],

    entry_points={
        'console_scripts': [
            'murcss=metrics.murcss:main',
        ],
    },
)
