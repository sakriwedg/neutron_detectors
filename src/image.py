# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 16:48:26 2026

@author: marchalj
"""
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.dates as mdates
from termcolor import colored
from matplotlib.backends.backend_pdf import PdfPages
#from scipy.optimize import curve_fit
#from dateutil import parser

def analysis(det,image_number,BKG_number, data_folder, reports_folder,sample):

    data,file_time,BM_sum=det.importPOS(data_folder+str(image_number).zfill(6)+'.nxs')
    NAC=np.int32(data)
    data,file_time,BM_sum=det.importPOS(data_folder+str(BKG_number).zfill(6)+'.nxs')
    BKG=np.int32(data)
    report_name="NACALF sample"
    
    report_name=reports_folder+sample+'.pdf'
    
    with PdfPages(report_name) as pdf:
        
        plt.rc('xtick', labelsize=6)
        plt.rc('ytick', labelsize=6)
        
        fig, axs = plt.subplots(3,1)
        plt.suptitle(det.detname, x=0.512, y=0.99, fontsize=8, ha='center')
        #plt.gcf().text(0.512, 0.92,data_file+'\n'+file_time, ha='center', fontsize=6)
        plt.gcf().text(0.512, 0.93,sample, ha='center', fontsize=6,color='red')
        im0=axs[0].imshow(NAC.T,origin='lower',interpolation = 'none',aspect=1/det.det_aspect_ratio)
        axs[0].set_title(str(image_number), fontsize=6)
        axs[0].set_xlabel('X',fontsize=6)
        axs[0].set_ylabel('Y',fontsize=6)
    
        im1=axs[1].imshow(BKG.T,origin='lower',interpolation = 'none',aspect=1/det.det_aspect_ratio)
        axs[1].set_title(str(BKG_number), fontsize=6)
        axs[1].set_xlabel('X',fontsize=6)
        axs[1].set_ylabel('Y',fontsize=6)
    
        im1=axs[2].imshow((NAC-BKG).T,origin='lower',interpolation = 'none',aspect=1/det.det_aspect_ratio)
        axs[2].set_title(str(image_number)+'-'+str(BKG_number), fontsize=6)
        axs[2].set_xlabel('X',fontsize=6)
        axs[2].set_ylabel('Y',fontsize=6)

        plt.tight_layout()
        pdf.savefig()
        plt.show()
