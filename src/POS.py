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


def uniformity(det, run_number, data_folder, reports_folder, show_figs=True):

    # This function processes the data of a single file and produces a counting uniformity analysis

    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    plt.ion()
    plt.style.use('bmh')
    
    data_file = data_folder + str(run_number).zfill(6) + '.nxs'
    data,file_time,BM_sum = det.importPOS(data_file)

    report_name=reports_folder+str(run_number)+'_counting_uniformity.pdf'

    def process(data):
        
        
        counts=(np.sum(data,0))
        if show_figs==True:
        
        
            with PdfPages(report_name) as pdf:
                
                plt.rc('xtick', labelsize=6)
                plt.rc('ytick', labelsize=6)
                
                fig, axs = plt.subplots(1,2)
                plt.suptitle(det.detname, x=0.512, y=0.99, fontsize=8, ha='center')
                plt.gcf().text(0.512, 0.92,data_file+'\n'+file_time, ha='center', fontsize=6)
                plt.gcf().text(0.512, 0.89,'2D image', ha='center', fontsize=6,color='red')
                im0=axs[0].imshow(data,origin='lower',interpolation = 'none',aspect=det.det_aspect_ratio)
                axs[0].set_title('Linear scale', fontsize=6)
                axs[0].set_xlabel('X',fontsize=6)
                axs[0].set_ylabel('Y',fontsize=6)

                
                im1=axs[1].imshow(np.log(data),origin='lower',interpolation = 'none',aspect=det.det_aspect_ratio)
                axs[1].set_title('Log scale', fontsize=6)
                axs[1].set_xlabel('Y',fontsize=6)
                axs[1].set_ylabel('Y',fontsize=6)
            
                plt.tight_layout()
                pdf.savefig()
                plt.show()
                
                # Projections
                
                fig, axs = plt.subplots(1,2)
                plt.suptitle(det.detname, x=0.512, y=0.99, fontsize=8, ha='center')
                plt.gcf().text(0.512, 0.92,data_file+'\n'+file_time, ha='center', fontsize=6)
                plt.gcf().text(0.512, 0.89,'Projections', ha='center', fontsize=6,color='red')
                axs[0].plot(np.sum(data,0))
                axs[0].set_title('Projection along X', fontsize=6)
                axs[0].set_xlabel('X',fontsize=6)
                axs[0].set_ylabel('Y',fontsize=6)

                axs[1].plot(np.sum(data,1))
                axs[1].set_title('Projection along Y', fontsize=6)
                axs[1].set_xlabel('Y',fontsize=6)
                axs[1].set_ylabel('Y',fontsize=6)
            
                plt.tight_layout()
                pdf.savefig()
                plt.show()
                
                # Counting profile per tube
                
                fig, axs = plt.subplots(1,2)
                plt.suptitle(det.detname, x=0.512, y=0.99, fontsize=8, ha='center')
                plt.gcf().text(0.512, 0.92,data_file+'\n'+file_time, ha='center', fontsize=6)
                plt.gcf().text(0.512, 0.89,'Projections', ha='center', fontsize=6,color='red')
                for tube in range(np.size(data,1)):
                    axs[0].plot(data[:,tube])
                    yn=data[:,tube]
                    x=np.arange(np.size(data,0))
                    
                    popt, _ = curve_fit(gaussian, x, yn)
                    ym = gaussian(x, *popt)
            
                                        
                axs[0].set_title('Projection along X', fontsize=6)
                axs[0].set_xlabel('X',fontsize=6)
                axs[0].set_ylabel('Counts',fontsize=6)

                axs[1].plot(np.sum(data,1))
                axs[1].set_title('Projection along Y', fontsize=6)
                axs[1].set_xlabel('Y',fontsize=6)
                axs[1].set_ylabel('Counts',fontsize=6)
            
                plt.tight_layout()
                pdf.savefig()
                plt.show()
                
   
        return counts,file_time

    counts,file_time=process(data)
               
    return counts,file_time,BM_sum # To be adapted for indv. read-out 
                
                

def stability(det,run_numbers,data_folder,reports_folder,show_figs=False):

    counts_list=[]
    file_time_list=[]
    BM_sum_list=[]
    for run in run_numbers:
        counts,file_time,BM_sum=uniformity(det, run, data_folder, reports_folder, show_figs=False)
        counts_list.append(counts)
        file_time_list.append(parser.parse(file_time))
        BM_sum_list.append(BM_sum)
    counts_array=np.array(counts_list)
    file_time_array=np.array(file_time_list)
    BM_sum_array=np.array(BM_sum_list)

    report_name=reports_folder+str(run_numbers[0])+'_'+str(run_numbers[-1])+'_counting_stability.pdf'

    with PdfPages(report_name) as pdf:
        
        plt.rc('xtick', labelsize=6)
        plt.rc('ytick', labelsize=6)

        fig, axs = plt.subplots(2,2)
        plt.suptitle(det.detname + ' Counting stability', x=0.512, y=0.99, fontsize=8, ha='center')
        [axs[0,0].plot(file_time_array,counts_array[:,n]) for n in range(np.size(counts_array,1))]
        axs[0,0].set_title('Counts per channel', fontsize=6)
        axs[0,0].set_xlabel('File number',fontsize=6)
        axs[0,0].set_ylabel('Counts',fontsize=6)
        x = axs[0,0].xaxis
        for item in x.get_ticklabels():
            item.set_rotation(45)
        
        [axs[0,1].plot(file_time_array,((100* (counts_array-np.mean(counts_array,0))) /np.mean(counts_array,0))[:,n]) for n in range(np.size(counts_array,1))]
        axs[0,1].set_title('Counts variation per channel', fontsize=6)
        axs[0,1].set_xlabel('File number',fontsize=6)
        axs[0,1].set_ylabel('%',fontsize=6)
        x = axs[0,1].xaxis
        for item in x.get_ticklabels():
            item.set_rotation(45)
            
        axs[1,0].plot(file_time_array,BM_sum_array)
        axs[1,0].set_title('Beam monitor counts per channel', fontsize=6)
        axs[1,0].set_xlabel('File number',fontsize=6)
        axs[1,0].set_ylabel('Counts',fontsize=6)
        x = axs[1,0].xaxis
        for item in x.get_ticklabels():
            item.set_rotation(45)
        
        axs[1,1].plot(file_time_array,100*(BM_sum_array-np.mean(BM_sum_array)) / np.mean(BM_sum_array)) 
        axs[1,1].set_title('Beam monitor counts variation per channel', fontsize=6)
        axs[1,1].set_xlabel('File number',fontsize=6)
        axs[1,1].set_ylabel('%',fontsize=6)
        x = axs[1,1].xaxis
        for item in x.get_ticklabels():
            item.set_rotation(45)
    
        plt.tight_layout()
        pdf.savefig()
        plt.show()
