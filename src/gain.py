# -*- coding: utf-8 -*-
"""
Created on Wed 26, 2026

@author: saenz-arevalo@ill.fr
-->
--> Gain analysis module for neutron detectors
"""

import re
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from termcolor import colored
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit
from matplotlib.colors import BoundaryNorm



def analysis(det, file_numbers, time, data_folder, reports_folder, board_number):
    

    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    plt.ion()
    plt.style.use('bmh')
    

    # counts rate (over PHS threshold) per tube vs time
    threshold_bin_ab = 25  #
    threshold_adc_ab  = threshold_bin_ab * det.ene_bin_size   
    print(colored(f"Using threshold at bin {threshold_bin_ab} → {threshold_adc_ab:.1f} adc","yellow"))


    data_ab, data_aa, data_bb, data_pos, data_gain, time_list = det.importListMode(file_numbers,time,data_folder, board_number)
    # data_ab, data_aa, data_bb  --> (file_time, Pulse Height    , tube)
    # data_pos                   --> (file_time, Position A/(A+B), tube)
    # data_gain                  --> (file_time, Pulse Height, Position, tube)

    report_name=reports_folder+str(file_numbers[0])+'-'+str(file_numbers[-1])+'_gain_maps.pdf'
    
    
    def process():

        print(colored('### Processing gain data...','green'))

        # total counts in all times vs tube (over PHS threshold)
        counts_vs_tube_ab = (data_ab[:, threshold_bin_ab:, :].sum(axis=1)).sum(axis=0)  # (tube,)
        
        # sum over time of pos vs tube 
        counts_vs_tube_pos = data_pos.sum(axis=0)  # (PH, tube)

        # gain vs tube
        gain_vs_tube_ab = np.zeros(det.det_ncols)  # (tube,)
        for tube in range(det.det_ncols):
            hist = data_ab[:, :, tube].sum(axis=0)  # (PH bins,)

            # apply threshold
            hist_sel = hist[threshold_bin_ab:]
            ph_sel = np.arange(threshold_bin_ab, det.n_bins_ene) * det.ene_bin_size

            if hist_sel.sum() > 0:
                gain_vs_tube_ab[tube] = np.sum(ph_sel * hist_sel) / np.sum(hist_sel)
            else:
                gain_vs_tube_ab[tube] = np.nan  

        dev_gain_vs_tube_ab = 100 * (gain_vs_tube_ab - np.nanmean(gain_vs_tube_ab)) / np.nanmean(gain_vs_tube_ab)  

        # average over all tubes:
        avg_gain_vs_tube_ab = np.nanmean(gain_vs_tube_ab)


        # gain vs (tube, position)

        phs_bins = np.arange(data_gain.shape[1])  # PHS bins


        gain_vs_tube_vs_pos   = np.zeros((det.n_bins_pos, det.det_ncols))  # (pos_bin, tube)
        counts_vs_tube_vs_pos = np.zeros((det.n_bins_pos, det.det_ncols))  # (pos_bin, tube)


        for tube in range(det.det_ncols):
            for pos_bin in range(det.n_bins_pos):

                hist = data_gain[:, :, pos_bin, tube].sum(axis=0)  # (PH bins,)

                # apply threshold
                hist_sel = hist[threshold_bin_ab:]
                ph_sel = phs_bins[threshold_bin_ab:]

                if hist_sel.sum() > 0:
                    gain_vs_tube_vs_pos[pos_bin, tube] = (
                        np.sum(ph_sel * det.ene_bin_size * hist_sel)
                        / np.sum(hist_sel)
                    )
                    counts_vs_tube_vs_pos[pos_bin, tube] = hist_sel.sum()
                else:
                    gain_vs_tube_vs_pos[pos_bin, tube] = np.nan
                    counts_vs_tube_vs_pos[pos_bin, tube] = 0

        # deviation from avg_gain_vs_tube_ab
        deviation_from_avg = 100 * (gain_vs_tube_vs_pos - avg_gain_vs_tube_ab) / avg_gain_vs_tube_ab


        with PdfPages(report_name) as pdf:

            plt.rc('xtick', labelsize=6)
            plt.rc('ytick', labelsize=6)

            ## --------------------------------------------------------------------------------##
            ##  1 Pos vs tube for all times
            ## --------------------------------------------------------------------------------##
            print(colored('Plotting pos vs tube for all times...','green'))

            fig, axs = plt.subplots(1,2, sharey=True)
            plt.suptitle(det.detname, x=0.512, y=0.99, fontsize=8, ha='center')
            plt.gcf().text(0.512, 0.92, data_folder+str(file_numbers[0])+'-'+str(file_numbers[-1]), ha='center', fontsize=6)
            plt.gcf().text(0.512, 0.89,'Pos vs tube for all times', ha='center', fontsize=6,color='red')


            im0 = axs[0].imshow(
                counts_vs_tube_pos, 
                origin='lower', 
                interpolation='nearest', 
                aspect='auto',
                extent=[
                    1, counts_vs_tube_pos.shape[1] + 1,                      # x: tube number
                    0, 350 * counts_vs_tube_pos.shape[0] * det.pos_bin_size / det.pos_max   # y: position bins → physical
                 ]
            )
            axs[0].set_title('Linear scale', fontsize=6)
            axs[0].set_xlabel('Tube number',fontsize=6)
            axs[0].set_ylabel('Pos bin',fontsize=6)     
            
            im1 = axs[1].imshow(
                np.log(counts_vs_tube_pos + 1), 
                origin='lower', 
                interpolation='nearest', 
                aspect='auto',
                extent=[
                    1, counts_vs_tube_pos.shape[1] + 1,                      # x: tube number
                    0, 350 * counts_vs_tube_pos.shape[0] * det.pos_bin_size / det.pos_max   # y: position bins → physical
                 ]
            )
    
            axs[1].set_title('Log scale', fontsize=6)
            axs[1].set_xlabel('Tube number',fontsize=6)
        
            plt.tight_layout(rect=[0, 0, 1, 0.95])
            pdf.savefig()
            plt.close(fig)

            ## --------------------------------------------------------------------------------##
            ##  2 Total counts vs tube for A+B
            ## --------------------------------------------------------------------------------##
            print(colored('Plotting total counts vs tube for A+B...','green'))

            fig, axs = plt.subplots(1,1)
            plt.suptitle(det.detname, x=0.512, y=0.99, fontsize=8, ha='center')
            plt.gcf().text(0.512, 0.92, data_folder+str(file_numbers[0])+'-'+str(file_numbers[-1]), ha='center', fontsize=6)
            plt.gcf().text(0.512, 0.89,'Total counts vs tube number', ha='center', fontsize=6,color='red')
        
            xticks = np.arange(0, det.ntubes + 2, 4)
    
            axs.set_xticks(xticks)
            axs.set_xticklabels(xticks + 1)


            axs.plot(np.arange(det.det_ncols), counts_vs_tube_ab, color='blue', alpha=0.7, linewidth=0.5, marker='o', markersize=3)
            axs.plot(counts_vs_tube_ab, color='blue', alpha=0.7, linewidth=0.5, marker='o', markersize=3)
            axs.set_xlabel('Tube number',fontsize=6)
            axs.set_ylabel('Total counts',fontsize=6)


            plt.tight_layout(rect=[0, 0, 1, 0.95])
            pdf.savefig()
            plt.close(fig)



            ## --------------------------------------------------------------------------------##
            ##  3 2D histogram of A+B PHS vs tube for all times
            ## --------------------------------------------------------------------------------##
            print(colored('Plotting 2D histogram of A+B PHS vs tube for all times...','green'))

            fig, axs = plt.subplots(1,2, sharey=True)
            plt.suptitle(det.detname, x=0.512, y=0.99, fontsize=8, ha='center')
            plt.gcf().text(0.512, 0.92, data_folder+str(file_numbers[0])+'-'+str(file_numbers[-1]), ha='center', fontsize=6)
            plt.gcf().text(0.512, 0.89,'A+B Pulse Height Spectra (PHS) for background in all files', ha='center', fontsize=6,color='red')
        
            im0 = axs[0].imshow(
                data_ab[:, :, :].sum(axis=0), 
                origin='lower', 
                interpolation='nearest', 
                aspect='auto',
                extent=[
                    1, data_ab.shape[2] + 1,                      # x: tube number
                    0, data_ab.shape[1] * det.ene_bin_size    # y: energy bins → physical
                 ]
                )
            axs[0].set_title('Linear scale', fontsize=6)
            axs[0].set_xlabel('Tube number',fontsize=6)
            axs[0].set_ylabel('PHS bin',fontsize=6)
            # adding red line to show threshold
            axs[0].axhline(threshold_adc_ab, color='red', linestyle='-', linewidth=0.9)   
            
            im1 = axs[1].imshow(
                np.log(data_ab[:, :, :].sum(axis=0) + 1), 
                origin='lower', 
                interpolation='nearest', 
                aspect='auto',
                extent=[
                    1, data_ab.shape[2] + 1,                      # x: tube number
                    0, data_ab.shape[1] * det.ene_bin_size    # y: energy bins → physical
                 ]
                )
            axs[1].set_title('Log scale', fontsize=6)
            axs[1].set_xlabel('Tube number',fontsize=6)
            axs[1].axhline(threshold_adc_ab, color='red', linestyle='-', linewidth=0.9)   
        
            plt.tight_layout(rect=[0, 0, 1, 0.95])
            pdf.savefig()
            plt.close(fig)


            
            ## --------------------------------------------------------------------------------##
            ##  4 Plot of gain vs tube and deviation from average gain vs tube
            ## --------------------------------------------------------------------------------##
            print(colored('4. Plotting gain vs tube and deviation from average gain vs tube...','green'))
            
            fig, axs = plt.subplots(1,2, sharey=False)
            plt.suptitle(det.detname, x=0.512, y=0.99, fontsize=8, ha='center')
            plt.gcf().text(0.512, 0.92, data_folder+str(file_numbers[0])+'-'+str(file_numbers[-1]), ha='center', fontsize=6)
            plt.gcf().text(0.512, 0.89,'', ha='center', fontsize=6,color='red')

            axs[0].plot(gain_vs_tube_ab, color='blue', alpha=0.7, linewidth=0.5, marker='o', markersize=3)
            axs[0].set_xlabel('Tube number',fontsize=6)
            axs[0].set_ylabel('Gain (mean PHS)',fontsize=6)

            xticks = np.arange(0, det.ntubes + 2, 4)
    
            axs[0].set_xticks(xticks)
            axs[0].set_xticklabels(xticks + 1)

            axs[1].set_xticks(xticks)
            axs[1].set_xticklabels(xticks + 1)

            axs[1].plot(
                dev_gain_vs_tube_ab, 
                color='blue', 
                alpha=0.7, 
                linewidth=0.5, 
                marker='o', 
                markersize=3,
            )
            
            axs[1].set_xlabel('Tube number',fontsize=6)
            axs[1].set_ylabel('Deviation from average gain (%)',fontsize=6)
            plt.tight_layout(rect=[0, 0, 1, 0.95])

            pdf.savefig()
            plt.close(fig)

                        
            ## --------------------------------------------------------------------------------##
            ##  5 Gain map (tube vs pos) for all times and deviation from average gain map. And counts per pixel map
            ## --------------------------------------------------------------------------------##
            print(colored('5. Plotting gain map (tube vs pos) for all times and deviation from average gain map...','green'))
            
            n_colors = 14

            # define limits
            vmin = np.nanmin(12000)
            vmax = np.nanmax(gain_vs_tube_vs_pos)
            bounds = np.linspace(vmin, vmax, n_colors + 1)
            cmap = plt.get_cmap('plasma', n_colors)
            norm = BoundaryNorm(bounds, cmap.N)

            fig, axs = plt.subplots(1,3, sharey=False)
            plt.suptitle(det.detname, x=0.512, y=0.99, fontsize=8, ha='center')
            plt.gcf().text(0.512, 0.92, data_folder+str(file_numbers[0])+'-'+str(file_numbers[-1]), ha='center', fontsize=6)
            plt.gcf().text(0.512, 0.89,'', ha='center', fontsize=6,color='red')

            low_limit = 0

            im0 = axs[0].imshow(
                gain_vs_tube_vs_pos,
                origin='lower',
                interpolation='nearest',
                aspect='auto',
                cmap=cmap,
                norm=norm,
                extent=[
                    1, gain_vs_tube_vs_pos.shape[1] + 1,
                    0, 350 * gain_vs_tube_vs_pos.shape[0] * det.pos_bin_size / det.pos_max
                ]
            )

            cbar = plt.colorbar(im0, ax=axs[0], boundaries=bounds)
            cbar.set_ticks(bounds)
            cbar.ax.tick_params(labelsize=6)

            # set x tics frequency to 4 for better readability
            axs[0].set_xticks(np.arange(1, gain_vs_tube_vs_pos.shape[1] , 8))
            #axs[0].set_ylim(low_limit, 2**16 - low_limit)
            axs[0].set_title('Gain map (mean PHS)', fontsize=6)
            axs[0].set_xlabel('Tube number',fontsize=6)
            axs[0].set_ylabel('Position',fontsize=6)


            # define limits
            vmin = np.nanmin(-7)
            vmax = np.nanmax(7)
            bounds = np.linspace(vmin, vmax, n_colors + 1)
            cmap = plt.get_cmap('gnuplot', n_colors)
            norm = BoundaryNorm(bounds, cmap.N)

            im1 = axs[1].imshow(
                deviation_from_avg, 
                origin='lower',
                interpolation='nearest',
                aspect='auto',
                cmap=cmap,
                norm=norm,
                extent=[
                    1, gain_vs_tube_vs_pos.shape[1] + 1,
                    0, 350 * gain_vs_tube_vs_pos.shape[0] * det.pos_bin_size / det.pos_max
                ]
            )
            cbar = plt.colorbar(im1, ax=axs[1], boundaries=bounds)
            cbar.set_ticks(bounds)
            cbar.ax.tick_params(labelsize=6)

            #axs[1].set_ylim(low_limit, 2**16 - low_limit)
            axs[1].set_title('Deviation from average gain (%)', fontsize=6)
            axs[1].set_xlabel('Tube number',fontsize=6)
            axs[1].set_ylabel('Position',fontsize=6)

            n_colors = 8

            # define limits
            vmin = np.nanmin(0)
            vmax = np.nanmax(40000)
            bounds = np.linspace(vmin, vmax, n_colors + 1)
            cmap = plt.get_cmap('plasma', n_colors)
            norm = BoundaryNorm(bounds, cmap.N)

            im2 = axs[2].imshow(
                counts_vs_tube_vs_pos, 
                origin='lower', 
                interpolation='nearest', 
                aspect='auto',
                cmap=cmap,
                norm=norm,
                extent=[
                    1, counts_vs_tube_vs_pos.shape[1] + 1,                      # x: tube number
                    0, 350 * counts_vs_tube_vs_pos.shape[0] * det.pos_bin_size / det.pos_max    # y: position bins → physical
                 ]
            )

            cbar = plt.colorbar(im2, ax=axs[2], boundaries=bounds)
            cbar.set_ticks(bounds)
            cbar.ax.tick_params(labelsize=6)

            #axs[2].set_ylim(low_limit, 2**16 - low_limit)

            axs[2].set_title('Counts map (over threshold)', fontsize=6)
            axs[2].set_xlabel('Tube number',fontsize=6)
            axs[2].set_ylabel('Position',fontsize=6)

            plt.tight_layout(rect=[0, 0, 1, 0.95])

            pdf.savefig()
            plt.close(fig)

           
    process()
    