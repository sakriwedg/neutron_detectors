# -*- coding: utf-8 -*-
"""
Created on May 13 2026

@author: saenz-arevalo@ill.fr
-->
--> Counting linearity analysis for slit measurement on CT2 with CSPEC prototype (April 2026)
"""

import re
import numpy as np
import sys
from termcolor import colored
import datetime       
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit

def analysis(det, file_numbers, dates, data_folder, reports_folder, board_number, stopping_thickness_mm, acquisition_time_s, roi_min_cm=0, roi_max_cm=0, analysis_tubes=None, position_among_analysis_tubes=0, initial_mu1_guess=0 ):
    

    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    plt.ion()
    plt.style.use('bmh')
    

    slit_data_ab, slit_data_aa, slit_data_bb, slit_data_pos, slit_data_ene, time_list = det.importListMode(file_numbers, dates, data_folder, board_number)
    # data_ab, data_aa, data_bb  --> (file_time, Pulse Height    , tube)
    # data_pos                   --> (file_time, Position A/(A+B), tube)
    # data_gain                  --> (file_time, Pulse Height, Position, tube)

    if slit_data_ab is None:
        print(colored('Error importing data. Skipping analysis.','red'))
        return


    max_number_of_thicknesses_to_plot = min(8, len(stopping_thickness_mm) - 1) # -1 because last one is background
    print(colored(f'Max number of thicknesses to analyze (excluding background): {max_number_of_thicknesses_to_plot}','cyan'))

    # check number of thicknessses are the same as number of files
    if len(stopping_thickness_mm) != len(file_numbers):
        print(colored('Error: number of stopping thicknesses does not match number of files','red'))
        print(colored(f' - stopping_thickness_mm: {stopping_thickness_mm}','yellow'))
        print(colored(f' - file_numbers: {file_numbers}','yellow'))
        return

    report_name=reports_folder+str(file_numbers[0])+'-'+str(file_numbers[-1])+'_counting_linearity.pdf'
        
    def process():


        background_slit_data_pos = slit_data_pos[-1]  # assuming last file is background
        background_slit_data_ab = slit_data_ab[-1]  # assuming last file is background


        with PdfPages(report_name) as pdf:

            plt.rc('xtick', labelsize=6)
            plt.rc('ytick', labelsize=6)

            ## --------------------------------------------------------------------------------##
            ##  1. 2D histogram of pos = A/(A+B) vs tube 
            ## --------------------------------------------------------------------------------##

            fig, axs = plt.subplots(4, 2, sharex=False)
            plt.suptitle(det.detname, x=0.512, y=0.99, fontsize=8, ha='center')
            plt.gcf().text(0.512, 0.92, data_folder+str(file_numbers[0])+'-'+str(file_numbers[-1]), ha='center', fontsize=6)
            plt.gcf().text(0.512, 0.89,'Position (log)', ha='center', fontsize=6,color='red')
        
            # plot first 8 only
            for thickness_idx in range(max_number_of_thicknesses_to_plot): # exclude background
                im0 = axs[thickness_idx // 2, thickness_idx % 2].imshow(
                    np.log( np.maximum((slit_data_pos[thickness_idx]+1) - background_slit_data_pos, 1) ).T, 
                    origin='lower',
                    cmap='inferno', 
                    interpolation='nearest', 
                    aspect='auto',
                    extent=[
                        0, 350 * slit_data_pos[thickness_idx].shape[0] * det.pos_bin_size / det.pos_max,# y: position bin
                        1, slit_data_pos[thickness_idx].shape[1] + 1                                    # x: tube number
                    ]
                )

                axs[thickness_idx // 2, thickness_idx % 2].set_title(f'Thickness: {stopping_thickness_mm[thickness_idx]} mm', fontsize=6)
                axs[thickness_idx // 2, thickness_idx % 2].set_ylabel('Tube number',fontsize=6)
                axs[thickness_idx // 2, thickness_idx % 2].set_xlabel('Pos (cm)',fontsize=6)
            
   
        
            plt.tight_layout(rect=[0, 0, 1, 0.95])
            pdf.savefig()
            plt.close(fig)

          
            ## --------------------------------------------------------------------------------##
            ##  2. 2D histogram of energy vs tube 
            ## --------------------------------------------------------------------------------##

            fig, axs = plt.subplots(4, 2, sharex=True)
            plt.suptitle(det.detname, x=0.512, y=0.99, fontsize=8, ha='center')
            plt.gcf().text(0.512, 0.92, data_folder+str(file_numbers[0])+'-'+str(file_numbers[-1]), ha='center', fontsize=6)
            plt.gcf().text(0.512, 0.89,'Energy', ha='center', fontsize=6,color='red')
        
            for thickness_idx in range(max_number_of_thicknesses_to_plot): # exclude background
                im0 = axs[thickness_idx // 2, thickness_idx % 2].imshow(
                    (slit_data_ab[thickness_idx].T - background_slit_data_ab.T), 
                    origin='lower', 
                    interpolation='nearest', 
                    aspect='auto',
                    extent=[
                        0, slit_data_ab[thickness_idx].shape[0] * det.ene_bin_size,  # y: energy bin
                        1, slit_data_ab[thickness_idx].shape[1] + 1                  # x: tube number
                    ]
                )

                axs[thickness_idx // 2, thickness_idx % 2].set_title(f'Thickness: {stopping_thickness_mm[thickness_idx]} mm', fontsize=6)
                axs[thickness_idx // 2, thickness_idx % 2].set_ylabel('Tube number',fontsize=6)
                axs[thickness_idx // 2, thickness_idx % 2].set_xlabel('Energy (a.u.)',fontsize=6)
            
            plt.tight_layout(rect=[0, 0, 1, 0.95])
            pdf.savefig()
            plt.close(fig)

            
            ## --------------------------------------------------------------------------------##
            ##  3. energy histogram (one per tube)
            ## --------------------------------------------------------------------------------##

            fig, axs = plt.subplots(2, 1, sharex=False)
            plt.suptitle(det.detname, x=0.512, y=0.99, fontsize=8, ha='center')
            plt.gcf().text(0.512, 0.92, data_folder+str(file_numbers[0])+'-'+str(file_numbers[-1]), ha='center', fontsize=6)
            plt.gcf().text(0.512, 0.89,'Energy', ha='center', fontsize=6,color='red')
        

            for thickness_idx in range(slit_data_ene.shape[0]-1): # exclude background
                energy_projection = slit_data_ab[thickness_idx][:, analysis_tubes[position_among_analysis_tubes]] - background_slit_data_ab[:, analysis_tubes[position_among_analysis_tubes]]

                # build energy axis
                energy_bins = np.arange(0, slit_data_ab[thickness_idx].shape[0]) * det.ene_bin_size
                

                axs[0].plot(
                    energy_bins,
                    energy_projection,
                    drawstyle='steps-mid',
                    label=f'{stopping_thickness_mm[thickness_idx]} mm',
                    linewidth=0.75
                )
                axs[0].set_xlabel('Energy (a.u.)', fontsize=6)
                axs[0].set_ylabel('Counts', fontsize=6)
                axs[0].legend(fontsize=6, ncol=2, handlelength=2)

                axs[1].plot(
                    energy_bins,
                    np.log(energy_projection + 1),
                    drawstyle='steps-mid',
                    label=f'{stopping_thickness_mm[thickness_idx]} mm',
                    linewidth=0.75
                )
                axs[1].set_xlabel('Energy (a.u.)', fontsize=6)
                axs[1].set_ylabel('Log(Counts)', fontsize=6)
                axs[1].legend(fontsize=6, ncol=2, handlelength=2)

            plt.tight_layout(rect=[0, 0, 1, 0.95])
            pdf.savefig()
            plt.close(fig)

            ## --------------------------------------------------------------------------------##
            ##  4. count rate over all energies (y-axis) vs tube number (x-axis) for each thickness.
            ## --------------------------------------------------------------------------------##

            fig, axs = plt.subplots(1, 1, sharex=True)
            plt.suptitle(det.detname, x=0.512, y=0.99, fontsize=8, ha='center')
            plt.gcf().text(0.512, 0.92, data_folder+str(file_numbers[0])+'-'+str(file_numbers[-1]), ha='center', fontsize=6)
            plt.gcf().text(0.512, 0.89,'Count rate', ha='center', fontsize=6,color='red')

            for thickness_idx in range(slit_data_ene.shape[0]-1): # exclude background
                count_rate_per_tube = (np.sum(slit_data_ab[thickness_idx], axis=0) - np.sum(background_slit_data_ab, axis=0)) / acquisition_time_s

                axs.plot(
                    np.arange(count_rate_per_tube.size) + 1,
                    count_rate_per_tube,
                    drawstyle='steps-mid',
                    label=f'{stopping_thickness_mm[thickness_idx]} mm'
                )
                axs.set_xlabel('Tube number', fontsize=6)
                axs.set_ylabel('Count rate (Hz)', fontsize=6)
                axs.legend(fontsize=6, ncol=4, handlelength=2)

            axs.set_yscale('log') 
            plt.tight_layout(rect=[0, 0, 1, 0.95])
            pdf.savefig()
            plt.close(fig)

            ## --------------------------------------------------------------------------------##
            ##  5. Count rate (y-axis) in analysis tube vs stopping thickness (x-axis).
            ## --------------------------------------------------------------------------------##
            
            fig, axs = plt.subplots(1, 1, sharex=True)
            plt.suptitle(det.detname, x=0.512, y=0.99, fontsize=8, ha='center')
            plt.gcf().text(0.512, 0.92, data_folder+str(file_numbers[0])+'-'+str(file_numbers[-1]), ha='center', fontsize=8)
            plt.gcf().text(0.512, 0.89,f'Count rate on tube {analysis_tubes[position_among_analysis_tubes]+1}', ha='center', fontsize=8,color='red')

            this_fontsize = 10

            count_rate_in_analysis_tubes = [[] for _ in analysis_tubes]
            count_rate_in_analysis_tubes_error = [[] for _ in analysis_tubes]
            for thickness_idx in range(slit_data_ene.shape[0]-1): # exclude background
                for tube_idx, tube in enumerate(analysis_tubes):
                    count_rate       =        (np.sum(slit_data_ab[thickness_idx][:, tube]) - np.sum(background_slit_data_ab[:, tube])) / acquisition_time_s
                    count_rate_error = np.sqrt(np.sum(slit_data_ab[thickness_idx][:, tube]) + np.sum(background_slit_data_ab[:, tube])) / acquisition_time_s
                    count_rate_in_analysis_tubes[tube_idx].append(count_rate)
                    count_rate_in_analysis_tubes_error[tube_idx].append(count_rate_error)


            for tube_idx, tube in enumerate(analysis_tubes):
                # plot points with error bars
                axs.errorbar(
                    stopping_thickness_mm[:len(count_rate_in_analysis_tubes[tube_idx])],
                    count_rate_in_analysis_tubes[tube_idx],
                    yerr=count_rate_in_analysis_tubes_error[tube_idx],
                    marker='o',
                    linestyle='',
                    label=f'Data'
                    #label=f'Tube {tube+1}'
                )
            
            # Compute exponential fit to count rate vs stopping thickness for each tube
            for tube_idx, tube in enumerate(analysis_tubes):
                if len(count_rate_in_analysis_tubes[tube_idx]) > 1:
                    # Fit an exponential function to the data (including error bars in the fit)
                    def exponential(x, A, k):
                        return A * np.exp(k * x)
                    popt, pcov = curve_fit(
                        exponential,
                        stopping_thickness_mm[:len(count_rate_in_analysis_tubes[tube_idx])],
                        count_rate_in_analysis_tubes[tube_idx],
                        sigma=count_rate_in_analysis_tubes_error[tube_idx],
                        absolute_sigma=True,
                        maxfev=10000
                    )
                    

                    axs.plot(
                        stopping_thickness_mm[:len(count_rate_in_analysis_tubes[tube_idx])],
                        [popt[0] * np.exp(popt[1] * x) for x in stopping_thickness_mm[:len(count_rate_in_analysis_tubes[tube_idx])]],
                        linestyle='--',
                        color='red',
                        label=f'Fit: {popt[0]:.1f} * exp({popt[1]:.1f} * x)'
                    )
            

            axs.set_yscale('log')
            axs.set_xlabel('Stopping thickness (mm)', fontsize=this_fontsize)
            axs.set_ylabel('Count rate (Hz)', fontsize=this_fontsize)
            axs.legend(fontsize=this_fontsize, ncol=1, handlelength=3)
            # Set tics fontsize
            axs.tick_params(axis='both', which='major', labelsize=this_fontsize)
            # set ytics labels to show only 1, 10, 100, 1000, etc
            axs.set_ylim(1e3, 1e5)
  

            plt.tight_layout(rect=[0, 0, 1, 0.95])
            pdf.savefig()
            plt.close(fig)

            
            ## --------------------------------------------------------------------------------##
            ##  6. Position on selected tubes 
            ## --------------------------------------------------------------------------------##
            
            fig, axs = plt.subplots(4, 2, sharex=False)
            plt.suptitle(det.detname, x=0.512, y=0.99, fontsize=8, ha='center')
            plt.gcf().text(0.512, 0.92, data_folder+str(file_numbers[0])+'-'+str(file_numbers[-1]), ha='center', fontsize=6)
            plt.gcf().text(0.512, 0.89,f'Position projection on tube {analysis_tubes[position_among_analysis_tubes]+1}', ha='center', fontsize=6,color='red')


            # plot first 8 only
            for thickness_idx in range(max_number_of_thicknesses_to_plot): # exclude background

                pos_projection = slit_data_pos[thickness_idx][:, analysis_tubes[position_among_analysis_tubes]] - background_slit_data_pos[:, analysis_tubes[position_among_analysis_tubes]]

                # build physical position axis
                pos_axis = 350 * np.arange(pos_projection.size) * det.pos_bin_size / det.pos_max

                axs[thickness_idx // 2, thickness_idx % 2].plot(
                    pos_axis,
                    pos_projection,
                    drawstyle='steps-mid'
                )


                axs[thickness_idx // 2, thickness_idx % 2].set_title(f'Thickness: {stopping_thickness_mm[thickness_idx]} mm', fontsize=6)
                axs[thickness_idx // 2, thickness_idx % 2].set_ylabel('Counts',fontsize=6)
                axs[thickness_idx // 2, thickness_idx % 2].set_xlabel('Pos (cm)',fontsize=6)
            
   

            plt.tight_layout(rect=[0, 0, 1, 0.95])
            pdf.savefig()
            plt.close(fig)

            
            ## --------------------------------------------------------------------------------##
            ##  7. Position in cm on selected tubes 
            ## --------------------------------------------------------------------------------##

            amplitude_fitted_gaussians = []
            total_counts_outside_roi = []
            resolution_list = []
            counts_in_roi_list = []

            # ----------------------------
            # Single Gaussian model
            # ----------------------------
            def single_gaussian(x, A, mu, sigma):
                return (
                    A * np.exp(-(x - mu)**2 / (2 * sigma**2))
                )
            
            # ----------------------------
            # Double Gaussian model
            # ----------------------------
            def double_gaussian(x, A1, mu1, sigma1, A2, mu2, sigma2):
                return (
                    A1 * np.exp(-(x - mu1)**2 / (2 * sigma1**2)) +
                    A2 * np.exp(-(x - mu2)**2 / (2 * sigma2**2))
                )

            # ----------------------------
            # Figure
            # ----------------------------
            fig, axs = plt.subplots(4, 2, sharex=False)
            plt.suptitle(det.detname, x=0.512, y=0.99, fontsize=8, ha='center')
            plt.gcf().text(0.512, 0.92, data_folder + str(file_numbers[0]) + '-' + str(file_numbers[-1]),
                        ha='center', fontsize=6)
            plt.gcf().text(0.512, 0.89,
                        'Position projection on selected tubes (cm)',
                        ha='center', fontsize=6, color='red')

            # ----------------------------
            # Loop over tubes
            # ----------------------------

            for thickness_idx in range(slit_data_pos.shape[0] - 1):  # exclude background

                # ------------------------
                # Data
                # ------------------------
                pos_projection = slit_data_pos[thickness_idx][:, analysis_tubes[position_among_analysis_tubes]] - background_slit_data_pos[:, analysis_tubes[position_among_analysis_tubes]]
                pos_axis = 350 * np.arange(pos_projection.size) * det.pos_bin_size / det.pos_max

                line = axs[thickness_idx // 2, thickness_idx % 2].plot(
                    pos_axis,
                    pos_projection,
                    drawstyle='steps-mid',
                    label=f'{stopping_thickness_mm[thickness_idx]} mm'
                )[0]
                color = line.get_color()

                # ------------------------
                # ROI mask
                # ------------------------
                mask = (pos_axis >= roi_min_cm) & (pos_axis <= roi_max_cm)

                x_data = pos_axis[mask]
                y_data = pos_projection[mask]

                if len(x_data) < 5 or np.max(y_data) == 0:
                    print(f" - - Tube {analysis_tubes[position_among_analysis_tubes]}: skipped (bad ROI)")
                    continue

                # ------------------------
                # INITIAL GUESSES (IMPORTANT)
                # ------------------------

                # find two highest peaks in ROI
                mu1 = initial_mu1_guess
                A1 = np.max(y_data)
                sigma0 = (roi_max_cm - roi_min_cm) / 20

                p0 = [
                    A1, mu1, sigma0,
                ]

                # print initial guesses
                print(colored(
                    f" - thickness {stopping_thickness_mm[thickness_idx]} mm -  initial guesses: "
                    f"mu1_0={mu1:.2f}, "
                    f"sigma1_0={sigma0:.2f}, "
                    f"A1_0={A1:.2f}",
                    "yellow"
                ))

                # plot initial guesses
                x_fit_pre = np.linspace(roi_min_cm, roi_max_cm, 800)
                y_fit_pre = single_gaussian(x_fit_pre, *p0)
                axs[thickness_idx // 2, thickness_idx % 2].plot(
                    x_fit_pre,
                    y_fit_pre,
                    '--',
                    color='gray',
                    linewidth=0.5
                )


                # ------------------------
                # FIT
                # ------------------------
                
                try:
                    popt, pcov = curve_fit(
                        single_gaussian,
                        x_data,
                        y_data,
                        p0=p0,
                        maxfev=20000
                    )

                    (A1, mu1, sigma1) = popt

                    print(colored(
                        f" - thickness {stopping_thickness_mm[thickness_idx]} mm: "
                        f"mu1={mu1:.2f}, sigma1={sigma1:.2f}, A1={A1:.2f}",
                        "cyan"
                    ))

                    # print sum of counts in ROI
                    print(colored(
                        f"   Total counts in ROI: {np.sum(y_data)}",
                        "cyan"
                    ))
                    # Print total counts in tube
                    print(colored(
                        f"   Total counts in tube: {np.sum(pos_projection)}",
                        "cyan"
                    ))

                    counts_in_roi_list.append(np.sum(y_data))

                    # ------------------------
                    # Plot fit (ROI only)
                    # ------------------------
                    x_fit = np.linspace(roi_min_cm, roi_max_cm, 800)
                    y_fit = single_gaussian(x_fit, *popt)

                    axs[thickness_idx // 2, thickness_idx % 2].plot(
                        x_fit,
                        y_fit,
                        '--',
                        color='red',
                        label=f'FWHM={2.355*sigma1:.1f}'
                    )

                    resolution_list.append((stopping_thickness_mm[thickness_idx], 2.355 * sigma1))
                    amplitude_fitted_gaussians.append(A1)
                    total_counts_outside_roi.append(np.sum(pos_projection) - np.sum(y_data))

                except RuntimeError:
                    print(colored(f" - {stopping_thickness_mm[thickness_idx]} mm: fit failed", "red"))
                    continue

                # ----------------------------
                # Final plot formatting
                # ----------------------------
                axs[thickness_idx // 2, thickness_idx % 2].set_xlabel('Position (cm)', fontsize=6)
                axs[thickness_idx // 2, thickness_idx % 2].set_ylabel('Counts', fontsize=6)

                axs[thickness_idx // 2, thickness_idx % 2].set_xlim(roi_min_cm - 5, roi_max_cm + 5)

                axs[thickness_idx // 2, thickness_idx % 2].legend(fontsize=6)
                axs[thickness_idx // 2, thickness_idx % 2].legend(
                    fontsize=6,
                    handlelength=3,
                    handletextpad=0.8
                )

            plt.tight_layout(rect=[0, 0, 1, 0.95])
            pdf.savefig()
            plt.close(fig)

            



       

    process()
    