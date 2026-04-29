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

#%% Selection of analysis functions
os = 'MAC' # Select the operating system on which the python scripts are run: os = "Windows" or os = "MAC"
do_gain_analysis = False
do_mean_gain_uniformity_analysis=True
do_slit_analysis = True
do_counting_stability_analysis = True
do_gain_stability_analysis=True
do_image_analysis=True

#%% Detector data files with ILL NOMAD + CAEN DT1740 modules
match os:
    case 'Windows':
        reports_folder=str((Path.cwd().parent).joinpath('reports'))+'/' # Reports folders
        data_folder="//serdon/illdata/data/ct2/exp_TEST-3519/rawdata/"  # Data folder
    case 'MAC':
        reports_folder=str((Path.cwd().parent).joinpath('reports'))+'/' # Reports folders
        data_folder="/Volumes/illdata/data/ct2/exp_TEST-3519/rawdata/"  # Data folder

CSPEC_gain_stability_nxs_file_numbers=np.arange(44380,44451) # Data files for tube Pulse Height Spectrum stability 
#CSPEC_gain_stability_nxs_file_numbers=np.arange(44380,44382) # Data files for tube Pulse Height Spectrum stability 
CSPEC_PHS_data_filename=44380 # Data files for tube Pulse Height Spectrum homogeneity
CSPEC_POS_stability_nxs_file_numbers=np.arange(44242,44319) # Data for files for counting stability with direct beam on CT2
NAC_image_number=44335 # NAC powder diffraction data file
BKG_image_number=44337 # Background (No sample) data file
slit_data_file_list=['44203','44205','44208']
slit_phys_pos_list=[138.8 + 19.5, 60 + 19.5, 19.5] # cm
slit_roi_min_list=[70,295,460] # pixels
slit_roi_max_list=[120,310,470]# pixels 

# slit_data_file_list=['44456','44453']
# slit_phys_pos_list=[138.8 + 19.5, 19.5] # cm
# slit_roi_min_list=[70,460] # pixels
# slit_roi_max_list=[120,470]# pixels 



#%% Data analysis
print("Data folder:",data_folder )
print("Reports folder:",reports_folder )

if do_mean_gain_uniformity_analysis:
   PHS.uniformity(ILL.CSPEC(),CSPEC_PHS_data_filename,data_folder,reports_folder)

if do_image_analysis:
    image.analysis(ILL.CSPEC(),NAC_image_number,BKG_image_number, data_folder, reports_folder,'NAC_sample')

if do_gain_stability_analysis:
    PHS.stability (ILL.CSPEC(), CSPEC_gain_stability_nxs_file_numbers, data_folder, reports_folder)

if do_counting_stability_analysis:
    POS.stability(ILL.CSPEC(),CSPEC_POS_stability_nxs_file_numbers,data_folder,reports_folder,show_figs=False)

# To update with latest version from William
if do_gain_analysis:
    gain.analysis(ILL.CSPEC(), CSPEC_gain_file_numbers, CSPEC_gain_time, CSPEC_gain_data_folder, reports_folder)

if do_slit_analysis:
    Slit.FWHMvsPOS(ILL.CSPEC(),slit_data_file_list,data_folder, reports_folder,slit_phys_pos_list,slit_roi_min_list,slit_roi_max_list)
    


    
    