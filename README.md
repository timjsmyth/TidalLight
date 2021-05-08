# TidalLight
A model capable of generating the solar, lunar, and artificial lightscape - in PAR (Broadband) OR Skye meter wavelengths R: 620-740nm, G: 495-560nm, B: 400-500nm (Spectral) - at and below surface of (theoretically) any global coastal position. 

The model requires a location with associated tidal gauge data (recommend https://www.psmsl.org/ as a starting point) to permit accuracte modelling of tidal harmonics and subsequent tide levels. 

A reference datum is selected by the user (default 10%).
Defining a datum percetage (line 296, TidalLight_Model.py) determines the "Below Tide" values of the selected light source over the defined time period to calculate the total light received at this depth (Designed to mimic an area within the intertidal zone).

Datum: Specified as a distance below the lowest astronomical tide for the time frame input

Critical depth of the incoming light source can also be determined for validation with in-situ measurements in either Broadband or Spectral groups

*Note:* World_Atlas_2015.tif & NetCDF4 files, containing global ALAN & Kd information respectively, are not included in this repository due to file size limits. 
Use a local copy in the following locations for any spectral model directory chain:

Relevant to lines 81 & 85 of Falchi_Kd_Position.py
datadir = '../../ALAN_Map+Kd/' *# Relative path to directory that contains Kd files on local machine* 

Python 3.6 is required to run this model.

To run this model open a command line in the TidalLight repository (Linux or WSL):

`bash run_model_Zc.sh`

This will run an example of some of the functionalities of the model, to see all possible commands: 

`python3 TidalLight_Model.py -h`. 

**Resources:**

The following modules are used in the model, for a theoretical understanding/explanation of these calculations 
see the relevant materials (referenced in TidalLight_Model.py --> TidalLight/Papers/) and READMEs of the listed modules:

astropy: [![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

pysolar: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1461066.svg)](https://doi.org/10.5281/zenodo.1461066)

UTide: https://github.com/wesleybowman/UTide
