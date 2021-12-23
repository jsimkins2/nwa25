## Conda Environments

These are the Python anaconda environments I used for the code within this repository. The primary environment is gridtools.yml. One may need to use xgcmtools.yml for chlorophyll file generation on OSx. Sometimes, OSx gets confused when loading xESMF from gridtools.

## Installing a Python Environment from YML

`conda env create -f gridtools.yml`

## Installing HCtFlood from Github

A custom package used in my environments is called [HCtFlood](https://github.com/raphaeldussin/HCtFlood) created by Raphael Dussin. When sourcing the above YML files, anaconda won't be able to install HCtFlood without path. Please install HCtFlood after the environment is created using the following...

`pip install git+https://github.com/raphaeldussin/HCtFlood.git`

