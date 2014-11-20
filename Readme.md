MurCSS
=====
A Tool for Standardized Evaluation of Decadal Hindcast Systems

The tool calculates the Mean Squared Error Skill Score (MSESS) its
decomposition (Correlation + Conditional Bias) and the Continuous
Ranked Probability Skill Score (CRPSS) as proposed by Goddard et
al. [2013]. The MSESS of both models and the MSESS "between" the two
models (model versions) are calculated for different leadtimes. The
CRPSS is calculated for both models defined by the input parameters.

The main documentation can be found here:

[Usage and API-Description][local-docs] 

[Theoretical Background][homepage].

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

Support, Issues, Bugs
-
Please open an issue on GitHub or write an email to `sebastian.illing@met.fu-berlin.de`.

Changelog
-
1.6.1:
```
New Features included:
    - field mean analysis --> including plots
    - zonal mean analysis --> including plots
    - basic output option
```
License
-
Copyright (C) 2014 Sebastian Illing This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/

[local-docs]: https://rawgit.com/illing2005/murcss/master/doc/build/html/index.html
[sample data]: ./sample_data
[sample output]: ./sample_output
[unittests]: ./integration/tests
[homepage]: https://www-miklip.dkrz.de/about/murcss
[PyPI]: https://pypi.python.org/pypi/murcss
[HadCRUT]: http://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/download.html
