# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 16:48:26 2026

@author: marchalj
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from termcolor import colored
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit
from dateutil import parser


# Gaussian function
def gaussian(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

def fit(det, run_number, data_folder, reports_folder, show_figs,tubeNumber,roi_min,roi_max):
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    plt.ion()
    plt.style.use('bmh')
    
    data_file = data_folder + str(run_number).zfill(6) + '.nxs'
    data,file_time,BM_sum = det.importPOS(data_file)

    report_name=reports_folder+str(run_number)+'_slit_image.pdf'
    
    counts=(np.sum(data,0))
    if show_figs==True:
        with PdfPages(report_name) as pdf:
            
            plt.rc('xtick', labelsize=6)
            plt.rc('ytick', labelsize=6)
            
            fig, axs = plt.subplots(3,1)
            plt.suptitle(det.detname, x=0.512, y=0.99, fontsize=8, ha='center')
            plt.gcf().text(0.512, 0.92,data_file+'\n'+file_time, ha='center', fontsize=6)
            plt.gcf().text(0.512, 0.90,'Slit image', ha='center', fontsize=6,color='red')
            im0=axs[0].imshow(data.T,origin='lower',interpolation = 'none',aspect=1/det.det_aspect_ratio)
            axs[0].set_title('Linear scale', fontsize=6)
            axs[0].set_xlabel('X',fontsize=6)
            axs[0].set_ylabel('Y',fontsize=6)

            im1=axs[1].imshow(np.log(data.T),origin='lower',interpolation = 'none',aspect=1/det.det_aspect_ratio)
            axs[1].set_title('Log scale', fontsize=6)
            axs[1].set_xlabel('Y',fontsize=6)
            axs[1].set_ylabel('Y',fontsize=6)
    
            y=data[roi_min:roi_max,tubeNumber]
            x=np.arange(roi_min,roi_max)  
        
            axs[2].plot(x,y,'*')
    
            popt, _ = curve_fit(gaussian, x, y,p0=[max(y),roi_min+(roi_max-roi_min)/2.0,30])
            ym = gaussian(x, *popt)

            axs[2].plot(x,ym)
            axs[2].set_title('Projection',fontsize=6)
            
            plt.tight_layout()
            pdf.savefig()
            plt.show()
    
    return popt # To be adapted for indv. read-out 
                

def FWHMvsPOS(det, data_file_list ,data_folder, reports_folder,slit_phys_pos_list,slit_roi_min_list,slit_roi_max_list):

    popt=np.zeros([len(data_file_list),3])
    for n in range((len(data_file_list))):
        popt[n,:]=fit(det,data_file_list[n] , data_folder,reports_folder,show_figs=True,tubeNumber=5,roi_min=slit_roi_min_list[n],roi_max=slit_roi_max_list[n])

    # linear function
    def linear(x, a, b):
        return a * x + b
    
    y=np.array(popt[:,1])
    x=np.array(slit_phys_pos_list)
    
    [a,b], _ = curve_fit(linear, x, y)
    ym = linear(x, *[a,b])
    
    print("linear fit slope =" + str(1/a) + " pixels per cm")
    report_name=reports_folder+'FWHM.pdf'
    
    with PdfPages(report_name) as pdf:
    
        plt.figure()
        plt.rc('xtick', labelsize=6)
        plt.rc('ytick', labelsize=6)
        plt.plot(x,y,'*')
        

        plt.suptitle(det.detname, x=0.512, y=0.99, fontsize=8, ha='center')
        plt.gcf().text(0.512, 0.92,str(data_file_list)+'\n', ha='center', fontsize=6)
        plt.gcf().text(0.512, 0.90,'Slit image', ha='center', fontsize=6,color='red')
    
        plt.plot(x,ym)
        plt.xlabel('Slit position in cm', fontsize=6)
        plt.ylabel('Image slit positon in pixels', fontsize=6)
        plt.legend(['Data','Linear fit'], fontsize=6)
    
        for n in np.arange(len(data_file_list)):
            print('FWHM =' + str(round(abs(2.35*popt[n,2]/a),2)) + ' cm at position ' + str(slit_phys_pos_list[n]) + ' cm')
            plt.gcf().text(0.2, 0.2+n*0.05,'FWHM =' + str(round(abs(2.35*popt[n,2]/a),2)) + ' cm at position ' + str(slit_phys_pos_list[n]) + ' cm', fontsize=6)
        
        plt.tight_layout()
        pdf.savefig()
        plt.show()




    # pos=[138.8 + 19.5, 60 + 19.5, 19.5]
    # popt=np.zeros([3,3])
    # CSPEC_Slit_pos1_nxs_file_number='44203'
    # popt[0,:]=fit(det,CSPEC_Slit_pos1_nxs_file_number , data_folder,reports_folder,show_figs=True,tubeNumber=5,roi_min=70,roi_max=120)
    # CSPEC_Slit_pos2_nxs_file_number='44205'
    # popt[1,:]=fit(det,CSPEC_Slit_pos2_nxs_file_number , data_folder,reports_folder,show_figs=True,tubeNumber=5,roi_min=295,roi_max=310)
    # CSPEC_Slit_pos3_nxs_file_number='44208'
    # popt[2,:]=fit(det,CSPEC_Slit_pos3_nxs_file_number , data_folder,reports_folder,show_figs=True,tubeNumber=5,roi_min=460,roi_max=470)

    # popt[0,:]=fit(det,data_file_list[0] , data_folder,reports_folder,show_figs=True,tubeNumber=5,roi_min=70,roi_max=120)
    # popt[1,:]=fit(det,data_file_list[1] , data_folder,reports_folder,show_figs=True,tubeNumber=5,roi_min=295,roi_max=310)
    # popt[2,:]=fit(det,data_file_list[2] , data_folder,reports_folder,show_figs=True,tubeNumber=5,roi_min=460,roi_max=470)
    
    
    # pos=[138.8 + 19.5, 60 + 19.5, 19.5]
    # popt=np.zeros([3,3])
    # CSPEC_Slit_pos1_nxs_file_number='44456'
    # popt[0,:]=fit(det,CSPEC_Slit_pos1_nxs_file_number , data_folder,reports_folder,show_figs=True,tubeNumber=5,roi_min=70,roi_max=120)
    # CSPEC_Slit_pos2_nxs_file_number='44205'
    # popt[1,:]=fit(det,CSPEC_Slit_pos2_nxs_file_number , data_folder,reports_folder,show_figs=True,tubeNumber=5,roi_min=290,roi_max=320)
    # CSPEC_Slit_pos3_nxs_file_number='44453'
    # popt[2,:]=fit(det,CSPEC_Slit_pos3_nxs_file_number , data_folder,reports_folder,show_figs=True,tubeNumber=5,roi_min=450,roi_max=480)
    
    
