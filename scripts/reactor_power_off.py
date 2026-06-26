# -*- coding: utf-8 -*-
"""
CSPEC detector prototype beamtest on CT2 (April 2026)
Created on Mon Mar 16 16:48:26 2026
@author: marchalj
"""
#%% Module imports
import numpy as np
import sys
sys.path.append('../')
from pathlib import Path
from detectors import ILL # import detector configuration module
from src import PHS as PHS 
from src import POS as POS 
import warnings
warnings.filterwarnings("ignore")


board_number = 0 # check with lstheader128 on .lst files (not relevant for nxs files)

### ------------------------------------------------------------------####
### Selection of analysis functions
### ------------------------------------------------------------------####
do_single_gain_uniformity_analysis  = 1
do_global_gain_uniformity_analysis  = 1
do_counting_stability_analysis    = 0
do_gain_stability_analysis        = 0

### ------------------------------------------------------------------####
### Report output folder
### ------------------------------------------------------------------####
reports_folder=str((Path.cwd().parent).joinpath('reports'))+'/'

### ------------------------------------------------------------------####
### Path to raw input data (ILL NOMAD + CAEN DT1740 modules)
### ------------------------------------------------------------------####
data_folder="//serdon/illdata/data/ct2/exp_TEST-3519/rawdata/"  # Used on Windows by ILL staff
data_folder="/Volumes/illdata/data/ct2/exp_TEST-3519/rawdata/"  # Used on Mac by ILL staff
#data_folder="..." # Used on VISA/linux by external users

### ------------------------------------------------------------------####
### Data file numbers for different analyses
### ------------------------------------------------------------------####

### --------------
# Data files for tube Pulse Height Spectrum stability 
#PHS_stability_nxs = np.arange(44641,44740) 
PHS_stability_nxs = np.arange(45073,45755) 

### --------------
# Data files for tube Pulse Height Spectrum uniformity
PHS_single_uniformity_nxs = 45073
PHS_global_uniformity_nxs = np.arange(45073,45755) 

### --------------
# Data files for counting stability with direct beam on CT2
#POS_stability_nxs = np.arange(44641,44740)
POS_stability_nxs = np.arange(45073,45755) 
threshold_bin = 50
show_figs = False


### Data analysis
print("--> Input data folder:", data_folder)
print("--> Output reports folder:", reports_folder)


if do_single_gain_uniformity_analysis:
   print("--> Computing single gain uniformity analysis of " + str(PHS_single_uniformity_nxs) )
   PHS.uniformity(
       ILL.CSPEC(), 
       PHS_single_uniformity_nxs, 
       data_folder, 
       reports_folder
       )

if do_global_gain_uniformity_analysis:
   print("--> Computing global gain uniformity analysis of " + str(PHS_global_uniformity_nxs[0]) + " to " + str(PHS_global_uniformity_nxs[-1]) )
   PHS.global_uniformity(
       ILL.CSPEC(), 
       PHS_global_uniformity_nxs, 
       data_folder, 
       reports_folder
       )

if do_gain_stability_analysis:
    print("--> Computing gain stability analysis of " + str(PHS_stability_nxs[0]) + " to " + str(PHS_stability_nxs[-1]) )
    PHS.stability (
        ILL.CSPEC(), 
        PHS_stability_nxs, 
        data_folder, 
        reports_folder
        )

if do_counting_stability_analysis:
    print("--> Computing counting stability analysis of " + str(POS_stability_nxs[0]) + " to " + str(POS_stability_nxs[-1]) )
    POS.stability(
        ILL.CSPEC(),
        POS_stability_nxs,
        data_folder,
        reports_folder,
        show_figs,
        threshold_bin
        )



    
    