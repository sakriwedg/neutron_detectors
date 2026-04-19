# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 16:48:26 2026

@author: marchalj
"""
import numpy as np
import sys
sys.path.append('../')
from src import ndet_lib
from pathlib import Path
from detectors import ILL # import detector configuration module
from src import PHS as PHS # import data analysis module
from src import POS as POS # import data analysis module
from src import gain as gain
import warnings
warnings.filterwarnings("ignore")


do_gain_analysis = True


current_directory = Path.cwd()
print("Current working directory:", current_directory)
reports_folder=str((current_directory.parent).joinpath('reports'))+'/' 
print("Reports folder:",reports_folder )

data_folder="//serdon/illdata/data/ct2/internalUse/rawdata/"
data_folder="/Volumes/illdata/data/ct2/internalUse/rawdata/"

#CSPEC_gain_stability_nxs_file_numbers=np.arange(43490,43573)
#CSPEC_POS_stability_nxs_file_numbers=np.arange(43648,43729)

CSPEC_gain_stability_nxs_file_numbers=np.arange(43490,43492)
CSPEC_POS_stability_nxs_file_numbers=np.arange(43648,43650)
CSPEC_POS_uniformity_nxs_file_number='43763'


CSPEC_gain_file_numbers = [43645]
CSPEC_gain_time="/Volumes/shares/share_sdn/SDN_python/ndet_data_analysis/data/time_gain_files.log"
CSPEC_gain_data_folder='/Volumes/shares/share_sdn/SDN_python/ndet_data_analysis/data/'


#PHS.stability (ILL.CSPEC(), CSPEC_gain_stability_nxs_file_numbers, data_folder, reports_folder)
#POS.stability (ILL.CSPEC(), CSPEC_POS_stability_nxs_file_numbers , data_folder, reports_folder,show_figs=False)
#POS.uniformity(ILL.CSPEC(), CSPEC_POS_uniformity_nxs_file_number , data_folder, reports_folder)


if do_gain_analysis:
    gain.analysis(ILL.CSPEC(), CSPEC_gain_file_numbers, CSPEC_gain_time, CSPEC_gain_data_folder, reports_folder)