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
from dateutil import parser

show_figs=False

def uniformity(det, run_number, PHS_data_folder, reports_folder):

    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    plt.ion()
    plt.style.use('bmh')
    
    data_file = PHS_data_folder + str(run_number).zfill(6) + '.nxs'
    data,file_time=det.importPHS(data_file)
    
    report_name=reports_folder+str(run_number)+'_gain_uniformity_.pdf'

    def process(data,comment=''):
        
        bins=np.arange(np.size(data,0))
        bins_array=np.tile(bins.T,(np.size(data,1),1)).T

        mean_pulse_heights=np.round(sum((data*bins_array),0)/sum((data),0),1)
        counts=(np.sum(data,0))
                         
        # Calculate average gain and counts for the whole detector
        det_gain_mean=np.mean(mean_pulse_heights)
        det_counts_mean=np.mean(counts)
        
        # Calculate relative variance of gain and counts for the whole detector
        det_gain_rel_var=np.round(100*np.std(mean_pulse_heights)/det_gain_mean,2)
        det_counts_rel_var=np.round(100*np.std(np.sum(data,0))/det_counts_mean,2)
        
        if show_figs==True:
        
        
            with PdfPages(reports_folder+report_name+' PHS per '+comment+'.pdf') as pdf:
                
                plt.rc('xtick', labelsize=6)
                plt.rc('ytick', labelsize=6)
                
                fig, axs = plt.subplots(1,2)
                plt.suptitle(det.detname, x=0.512, y=0.99, fontsize=8, ha='center')
                plt.gcf().text(0.512, 0.92,data_file+'\n'+file_time, ha='center', fontsize=6)
                plt.gcf().text(0.512, 0.89,'Pulse Height Spectra (PHS) intensity map', ha='center', fontsize=6,color='red')
                im0=axs[0].imshow(data,origin='lower',interpolation = 'none',aspect=det.det_aspect_ratio)
                #plt.suptitle(det.detname+'\n'+PHS_data_file+'\n'+file_time+'\n'+'\n'+'Pulse Height Spectra (PHS) intensity map', fontsize=8, color='#444', ha='center')
                axs[0].set_title('Linear scale', fontsize=6)
                axs[0].set_xlabel(comment+' number',fontsize=6)
                axs[0].set_ylabel('PHS bin',fontsize=6)
                #plt.colorbar(im0,cax=axs[0])
                
                im1=axs[1].imshow(np.log(data),origin='lower',interpolation = 'none',aspect=det.det_aspect_ratio)
                axs[1].set_title('Log scale', fontsize=6)
                axs[1].set_xlabel(comment+' number',fontsize=6)
                axs[1].set_ylabel('PHS bin',fontsize=6)
            
                plt.tight_layout()
                pdf.savefig()
                plt.show()
                
                fig, axs = plt.subplots(1,2)
                plt.suptitle(det.detname, x=0.512, y=0.99, fontsize=8, ha='center')
                plt.gcf().text(0.512, 0.92,data_file+'\n'+file_time, ha='center', fontsize=6)
                plt.gcf().text(0.512, 0.89,'Pulse Height Spectrum (PHS) per channel', ha='center', fontsize=6,color='red')
                for tube in range(det.det_ncols):
                    axs[0].plot((data[:,tube]),color='m',linewidth=0.5)
                axs[0].set_title('Linear scale', fontsize=6)
                axs[0].set_ylabel('a.u.',fontsize=6)
                axs[0].set_xlabel('PHS bins',fontsize=6)
                axs[0].set_box_aspect(1)
                for tube in range(det.det_ncols):
                    axs[1].plot(np.log((data[:,tube])),color='c',linewidth=0.5)
                axs[1].set_title('Log scale', fontsize=6)
                axs[1].set_ylabel('a.u.',fontsize=6)
                axs[1].set_xlabel('PHS bins',fontsize=6)
                axs[1].set_box_aspect(1)
                plt.tight_layout()
                pdf.savefig()
                plt.show()
                
                fig, axs = plt.subplots(1,2)
                plt.suptitle(det.detname, x=0.512, y=0.99, fontsize=8, ha='center')
                plt.gcf().text(0.512, 0.92,data_file+'\n'+file_time, ha='center', fontsize=6)
                plt.gcf().text(0.512, 0.89,'Mean pulse height per channel', ha='center', fontsize=6,color='red')
                axs[0].plot(mean_pulse_heights,color='tab:blue')
                axs[0].set_title('Mean pulse height', fontsize=6)
                axs[0].set_ylabel('PHS bins',fontsize=6)
                axs[0].set_xlabel(comment +' number',fontsize=6)
                axs[0].set_box_aspect(1)
                axs[1].plot(100*(mean_pulse_heights-np.mean(mean_pulse_heights))/np.mean(mean_pulse_heights),color='tab:blue')
                axs[1].set_title('Gain variations (rel. std. = '+str(det_gain_rel_var) + '%)', fontsize=6)
                axs[1].set_ylabel('%',fontsize=6)
                axs[1].set_xlabel(comment+' number',fontsize=6)
                axs[1].set_box_aspect(1)
                axs[1].set_ylim(-2.5,2.5)
                plt.tight_layout()
                pdf.savefig()
                plt.show()

                
                fig, axs = plt.subplots(1,2)
                plt.suptitle(det.detname, x=0.512, y=0.99, fontsize=8, ha='center')
                plt.gcf().text(0.512, 0.92,data_file+'\n'+file_time, ha='center', fontsize=6)
                plt.gcf().text(0.512, 0.89,'PHS integrated for each channel', ha='center', fontsize=6,color='red')
                axs[0].plot(np.sum(data,0),color='tab:red')
                axs[0].set_title('Counts per tube', fontsize=6)
                axs[0].set_ylabel('Counts',fontsize=6)
                axs[0].set_xlabel(comment+' number',fontsize=6)
                axs[0].set_box_aspect(1)
                axs[1].plot(100*(np.sum(data,0)-np.mean(np.sum(data,0)))/np.mean(np.sum(data,0)),color='tab:red')
                axs[1].set_title('Counting variations (rel. std. = '+str(det_counts_rel_var) + '%)', fontsize=6)
                axs[1].set_ylabel('%',fontsize=6)
                axs[1].set_xlabel(comment+' number',fontsize=6)
                axs[1].set_box_aspect(1)
                axs[1].set_ylim(-50,50)
                plt.tight_layout()
                pdf.savefig()
                plt.show()
                
                
                if det.nmods!=1:
                    # Calculate average gain and average counts for each module
                    gain_mean=[np.mean(mean_pulse_heights[m*det.ncols:m*det.ncols+det.ncols]) for m in range(det.nmods)]
                    counts_mean=[np.mean(np.sum(PHSdata_ROI,0)[m*det.ncols:m*det.ncols+det.ncols]) for m in range(det.nmods)]
                    
                    # Calculate relative variance of gain and counting for each module
                    gain_rel_var=[100*np.std(mean_pulse_heights[m*det.ncols:m*det.ncols+det.ncols])/np.mean(mean_pulse_heights[m*det.ncols:m*det.ncols+det.ncols]) for m in range(det.nmods)]
                    counts_rel_var=[100*np.std(np.sum(PHSdata_ROI,0)[m*det.ncols:m*det.ncols+det.ncols])/np.mean(np.sum(PHSdata_ROI,0)[m*det.ncols:m*det.ncols+det.ncols]) for m in range(det.nmods)]
                  
            
                    fig, axs = plt.subplots(1,2)
                    plt.suptitle(det.detname, x=0.512, y=0.99, fontsize=8, ha='center')
                    plt.gcf().text(0.512, 0.92,PHS_data_file+'\n'+file_time, ha='center', fontsize=6)
                    plt.gcf().text(0.512, 0.89,'Gain and counting variations between modules', ha='center', fontsize=6,color='red')
                    
                    bars=axs[0].bar(np.arange(det.nmods)+0.5,100*((gain_mean-np.mean(gain_mean))/np.mean(gain_mean)),color='tab:blue', align="edge", animated=0.4)
                    for rect in bars:
                        height = rect.get_height()
                        axs[0].text(rect.get_x() + rect.get_width()/2.0, height, str(round(height,1)), ha='center', va='bottom',fontsize=6)
                        axs[0].set_title('Gain', fontsize=6)
                    #axs[0].set_ylim(0,1)
                    axs[0].set_ylabel('%',fontsize=6)
                    axs[0].set_xlabel('Module number',fontsize=6)
                    axs[0].set_box_aspect(1)
                    
                    bars=axs[1].bar(np.arange(det.nmods)+0.5,100*((counts_mean-np.mean(counts_mean))/np.mean(counts_mean)),color='tab:red', align="edge", animated=0.4)
                    for rect in bars:
                        height = rect.get_height()
                        axs[1].text(rect.get_x() + rect.get_width()/2.0, height, str(round(height,1)), ha='center', va='bottom',fontsize=6)
                        axs[1].set_title('Counts', fontsize=6)
                    axs[1].set_ylabel('%',fontsize=6)
                    axs[1].set_xlabel('Module number',fontsize=6)
                    axs[1].set_box_aspect(1)
                    plt.tight_layout()
                    pdf.savefig()
                    plt.show()
            
                    for m in range(det.nmods):
                        fig, axs = plt.subplots(1,2)
                        plt.suptitle(det.detname, x=0.512, y=0.99, fontsize=8, ha='center')
                        plt.gcf().text(0.512, 0.92,PHS_data_file+'\n'+file_time, ha='center', fontsize=6)
                        plt.gcf().text(0.512, 0.89,'Gain and counting variations within module #'+str(m+1), ha='center', fontsize=6,color='red')
                        axs[0].plot(100*((mean_pulse_heights[m*det.ncols:m*det.ncols+det.ncols]-np.mean(mean_pulse_heights[m*det.ncols:m*det.ncols+det.ncols]))/np.mean(mean_pulse_heights[m*det.ncols:m*det.ncols+det.ncols])),color='tab:blue')
                        axs[0].set_title('Gain variations (rel. std. = '+str(round(gain_rel_var[m],1)) + '%)', fontsize=6)
                        axs[0].set_ylabel('%',fontsize=6)
                        axs[0].set_xlabel('Tube number',fontsize=6)
                        axs[0].set_box_aspect(1)
                        axs[0].set_ylim(-10,10)
                        axs[1].plot(100*((np.sum(PHSdata_ROI,0)[m*det.ncols:m*det.ncols+det.ncols]-np.mean(np.sum(PHSdata_ROI,0)[m*det.ncols:m*det.ncols+det.ncols]))/np.mean(np.sum(PHSdata_ROI,0)[m*det.ncols:m*det.ncols+det.ncols])),color='tab:red')
                        axs[1].set_title('Counting variations (rel. std. = '+str(round(counts_rel_var[m],1)) + '%)', fontsize=6)
                        axs[1].set_ylabel('%',fontsize=6)
                        axs[1].set_xlabel('Tube number',fontsize=6)
                        axs[1].set_box_aspect(1)
                        axs[1].set_ylim(-10,10)
                        plt.tight_layout()
                        pdf.savefig()
                        plt.show()
            
                   
                    fig, axs = plt.subplots(1,2)
                    plt.suptitle(det.detname, x=0.512, y=0.99, fontsize=8, ha='center')
                    plt.gcf().text(0.512, 0.92,PHS_data_file+'\n'+file_time, ha='center', fontsize=6)
                    plt.gcf().text(0.512, 0.89,'Relative standard deviation of gains and counts within each module', ha='center', fontsize=6,color='red')
                    
                    bars=axs[0].bar(np.arange(det.nmods)+0.5,gain_rel_var,color='tab:blue', align="edge", animated=0.4)
                    for rect in bars:
                        height = rect.get_height()
                        axs[0].text(rect.get_x() + rect.get_width()/2.0, height, str(round(height,1)), ha='center', va='bottom',fontsize=6)
                        axs[0].set_title('Gain variations (rel. std. = '+str(det_gain_rel_var) + '%)', fontsize=6)
                    #axs[0].set_ylim(0,1)
                    axs[0].set_ylabel('rel. std. (%)',fontsize=6)
                    axs[0].set_xlabel('Module number',fontsize=6)
                    axs[0].set_box_aspect(1)
                    
                    bars=axs[1].bar(np.arange(det.nmods)+0.5,counts_rel_var,color='tab:red', align="edge", animated=0.4)
                    for rect in bars:
                        height = rect.get_height()
                        axs[1].text(rect.get_x() + rect.get_width()/2.0, height, str(round(height,1)), ha='center', va='bottom',fontsize=6)
                        axs[1].set_title('Counts variations (rel. std. = '+str(det_counts_rel_var) + '%)', fontsize=6)
                    axs[1].set_ylabel('rel. std. (%)',fontsize=6)
                    axs[1].set_xlabel('Module number',fontsize=6)
                    axs[1].set_box_aspect(1)
                    plt.tight_layout()
                    
                    pdf.savefig()
                    plt.show()
                    
                    del PHSdata_ROI
            

        return mean_pulse_heights,counts,file_time


    match det.readout:
        case 'charge division':
            mean_pulse_heights,counts,file_time=process(data,'Tube')

                
        case 'indiv. readout':
            if det.nrows==0:
                mean_pulse_heights,counts,file_time=process(data,'X channel')     
            else:
                mean_pulse_heights_X,counts_X,file_time=process(data[:,0:det.ncols],'X channel')
                mean_pulse_heights_Y,counts_Y,file_time=process(data[:,det.ncols:det.ncols+det.nrows],'Y channel')
                mean_pulse_heights=np.concatenate(mean_pulse_heights_X,mean_pulse_heights_Y,axis=0)
                counts=np.concatenate(counts_X,counts_Y,axis=0)

                
    return mean_pulse_heights,counts,file_time # To be adapted for indv. read-out 
                
                

def stability(det,run_numbers,data_folder,reports_folder):
    
    mean_pulse_heights_list=[]
    counts_list=[]
    file_time_list=[]
    for run in run_numbers:
        mean_pulse_heights,counts,file_time=uniformity(det, run, data_folder, reports_folder)
        mean_pulse_heights_list.append(mean_pulse_heights)
        counts_list.append(counts)
        file_time_list.append(parser.parse(file_time))
    mean_pulse_heights_array=np.array(mean_pulse_heights_list)
    counts_array=np.array(counts_list)
    file_time_array=np.array(file_time_list)

    report_name=reports_folder+str(run_numbers[0])+'_'+str(run_numbers[-1])+'_gain_stability_.pdf'

    with PdfPages(report_name) as pdf:
        
        plt.rc('xtick', labelsize=6)
        plt.rc('ytick', labelsize=6)
        fig, axs = plt.subplots(1,2)
        plt.suptitle(det.detname + ' Gain stability', x=0.512, y=0.99, fontsize=8, ha='center')
        [axs[0].plot(file_time_array,mean_pulse_heights_array[:,n]) for n in range(np.size(mean_pulse_heights_array,1))]
        axs[0].set_title('Mean pulse height per channel', fontsize=6)
        axs[0].set_xlabel('File number',fontsize=6)
        axs[0].set_ylabel('PHS bin',fontsize=6)
        x = axs[0].xaxis
        for item in x.get_ticklabels():
            item.set_rotation(45)
        
        [axs[1].plot(file_time_array,((100* (mean_pulse_heights_array-np.mean(mean_pulse_heights_array,0))) /np.mean(mean_pulse_heights_array,0))[:,n]) for n in range(np.size(mean_pulse_heights_array,1))]
        axs[1].set_title('Gain variation per channel', fontsize=6)
        axs[1].set_xlabel('File number',fontsize=6)
        axs[1].set_ylabel('%',fontsize=6)
        x = axs[1].xaxis
        for item in x.get_ticklabels():
            item.set_rotation(45)
        plt.tight_layout()
        pdf.savefig()
        plt.show()
    
        fig, axs = plt.subplots(1,2)
        plt.suptitle(det.detname + ' Counting stability', x=0.512, y=0.99, fontsize=8, ha='center')
        [axs[0].plot(file_time_array,counts_array[:,n]) for n in range(np.size(counts_array,1))]
        axs[0].set_title('Counts per channel', fontsize=6)
        axs[0].set_xlabel('File number',fontsize=6)
        axs[0].set_ylabel('Counts',fontsize=6)
        x = axs[0].xaxis
        for item in x.get_ticklabels():
            item.set_rotation(45)
        
        [axs[1].plot(file_time_array,((100* (counts_array-np.mean(counts_array,0))) /np.mean(counts_array,0))[:,n]) for n in range(np.size(mean_pulse_heights_array,1))]
        axs[1].set_title('Counts variation per channel', fontsize=6)
        axs[1].set_xlabel('File number',fontsize=6)
        axs[1].set_ylabel('%',fontsize=6)
        x = axs[1].xaxis
        for item in x.get_ticklabels():
            item.set_rotation(45)
    
        plt.tight_layout()
        pdf.savefig()
        plt.show()
