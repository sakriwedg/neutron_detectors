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
#data_folder="//serdon/illdata/data/ct2/exp_TEST-3519/rawdata/"  # Used on Windows by ILL staff
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
print("Data folder:",data_folder )
print("Reports folder:",reports_folder )

if do_mean_gain_uniformity_analysis:
   PHS.uniformity(
       ILL.CSPEC(), 
       PHS_uniformity_nxs, 
       data_folder, 
       reports_folder
       )

if do_image_analysis:
    image.analysis(
        ILL.CSPEC(),
        NAC_image_number,
        BKG_image_number,
        data_folder,
        reports_folder,
        'NAC_sample'
        )

if do_gain_stability_analysis:
    PHS.stability (
        ILL.CSPEC(), 
        PHS_stability_nxs, 
        data_folder, 
        reports_folder
        )

if do_counting_stability_analysis:
    POS.stability(
        ILL.CSPEC(),
        POS_stability_nxs,
        data_folder,
        reports_folder,
        show_figs=False
        )

# To update with latest version from William
if do_gain_analysis:
    gain.analysis(
        ILL.CSPEC(), 
        gain_file, 
        gain_time, 
        data_folder, 
        reports_folder,
        board_number
        )

if do_slit_analysis:
    Slit.FWHMvsPOS(
        ILL.CSPEC(),
        slit_data_nxs,
        data_folder,
        reports_folder,
        slit_phys_pos,
        slit_roi_min,
        slit_roi_max
        )
    


    
    