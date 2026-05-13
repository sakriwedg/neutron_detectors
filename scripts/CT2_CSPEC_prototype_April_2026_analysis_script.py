# -*- coding: utf-8 -*-
"""
CSPEC detector prototype beamtest on CT2 (April 2026)
Created on Mon Mar 16 16:48:26 2026
@author: marchalj
"""
#%% Module imports
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('../')
from src import ndet_lib
from pathlib import Path
from detectors import ILL # import detector configuration module
from src import PHS as PHS 
from src import POS as POS 
from src import gain as gain
from src import Slit as Slit 
from src import image as image 
import src
import warnings
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages
warnings.filterwarnings("ignore")


board_number = 0 # check with lstheader128 on .lst files (not relevant for nxs files)

### ------------------------------------------------------------------####
### Selection of analysis functions
### ------------------------------------------------------------------####
do_gain_analysis                  = True
do_mean_gain_uniformity_analysis  = True
do_slit_analysis                  = True
do_counting_stability_analysis    = True
do_gain_stability_analysis        = True
do_image_analysis                 = True

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

# Data files for tube Pulse Height Spectrum stability 
PHS_stability_nxs = np.arange(44380,44451) 
#PHS_stability_nxs = np.arange(44380,44382)

# Data files for tube Pulse Height Spectrum uniformity
PHS_uniformity_nxs = 44380 

# Data files for counting stability with direct beam on CT2
POS_stability_nxs = np.arange(44242,44319) 

# NAC powder diffraction data file and background (no sample)
NAC_image_number = 44335
BKG_image_number = 44337

# Data files for gain map
gain_file = [44353]
gain_time = ''
n_pos_bins = 16  

# Data files for spatial resolution analysis with slit at different positions
slit_data_nxs = ['44203','44205','44208']
slit_phys_pos = [138.8 + 19.5, 60 + 19.5, 19.5] # cm
slit_roi_min  = [70,295,460]                    # pixels
slit_roi_max  = [120,310,470]                   # pixels 

# slit_data_nxs = ['44456','44453']
# slit_phys_pos=[138.8 + 19.5, 19.5] # cm
# slit_roi_min=[70,460] # pixels
# slit_roi_max=[120,470]# pixels 



### Data analysis
print("--> Input data folder:", data_folder)
print("--> Output reports folder:", reports_folder)

if do_mean_gain_uniformity_analysis:
   print("--> Computing mean gain uniformity analysis of " + str(PHS_uniformity_nxs) )
   PHS.uniformity(
       ILL.CSPEC(), 
       PHS_uniformity_nxs, 
       data_folder, 
       reports_folder
       )

if do_image_analysis:
    print("--> Computing image analysis of " + str(NAC_image_number) )
    image.analysis(
        ILL.CSPEC(),
        NAC_image_number,
        BKG_image_number,
        data_folder,
        reports_folder,
        'NAC_sample'
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
        show_figs=False
        )

if do_gain_analysis:
    print("--> Computing gain map of " + str(gain_file) )
    gain.map(
        ILL.CSPEC(), 
        gain_file, 
        gain_time, 
        data_folder, 
        reports_folder,
        board_number,
        n_pos_bins
        )

if do_slit_analysis:
    print("--> Computing slit analysis of " + str(slit_data_nxs) )
    Slit.FWHMvsPOS(
        ILL.CSPEC(),
        slit_data_nxs,
        data_folder,
        reports_folder,
        slit_phys_pos,
        slit_roi_min,
        slit_roi_max
        )
    


    
    