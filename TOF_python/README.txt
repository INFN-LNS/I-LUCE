Here it is a script to perform Time of Flight analysis from .csv files
and extract the differential proton spectrum.
Notice that this analysis is performed in the hypotesis that only one population of ions is present: protons


"uti.py" collects some functions to perform file reading, computation of particle range in material and energy loss

"TOF_proton_spectrum.py" produce the differential proton spectrum

The "SRIM_Outputs" is a collection of files from SRIM simulation. It contains tabs of stopping power.
Those files are essential for the right computatiin of energy loss in Aluminum filters and in detectors of different materials.

