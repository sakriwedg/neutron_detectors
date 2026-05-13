# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 16:48:26 2026
Definition of ILL-PANTHER detector class including detector parameters and datafiles import methods

@author: marchalj
"""
### This generic ndet neutron detector class should be transfered in the main program in \src
import sys
sys.path.append('../../')
from src import ndet_lib
import numpy as np
import datetime
import h5py
import os
from datetime import datetime


def importnxs(filename):
    f=h5py.File(filename,'r')
    # check if file was correctly imported
    if f is None:
        print(f"Error: File '{filename}' could not be imported.")
        sys.exit(1)
    return f

def read_hex_file(file_path):
    hex_values=[]
    with open(file_path, 'r') as file:

        if file is None:
            print(f"Error: File '{file_path}' could not be opened.")
            sys.exit(1)

        for line in file:
            # Strip whitespace and check for comment character
            stripped_line = line.strip()
            if stripped_line.startswith('#'):
                break  # Skip comment lines
            # Split line into parts and process hexadecimal numbers
            parts = stripped_line.split()
            for part in parts:
                if part:  # Ensure part is not empty
                    try:
                        # Convert hexadecimal to integer
                        hex_values.append(int(part, 16))
                        #print(f"Hexadecimal: {part}, Decimal: {hex_value}")
                    except ValueError:
                        print(f"Invalid hexadecimal number: {part}")
    return hex_values

def import_MCC_cluster_dump(self,filename):
    hex_values=read_hex_file(filename)
    if len(hex_values)==0:
        print(f"Error: No hexadecimal values were read from the file '{filename}'.")
        sys.exit(1)

    T=0.013# MCC clock period in microseconds
    sepValue=16383
    dumpArray=np.asarray(hex_values)
    sepIdx=np.where(dumpArray==sepValue)[0]
    eventWidth=np.diff(sepIdx)-1
    multiMax=int(eventWidth.max()/3)
    sepIdxPerMulti=[sepIdx[np.where(eventWidth==(multi+1)*3)] for multi in range(multiMax)]
    ToAPerMulti=[[T*dumpArray[sepIdxPerMulti[multi]+1+n*3] for n in range(multi+1)]  for multi in range(multiMax)]
    TOTPerMulti=[[T*dumpArray[sepIdxPerMulti[multi]+2+n*3] for n in range(multi+1)]  for multi in range(multiMax)]
    ChNbPerMulti=[[dumpArray[sepIdxPerMulti[multi]+3+n*3] for n in range(multi+1)]  for multi in range(multiMax)]
    return ToAPerMulti,TOTPerMulti,ChNbPerMulti,sepIdxPerMulti,multiMax
    
def createPHS(self,ToAPerMulti,TOTPerMulti,ChNbPerMulti,sepIdxPerMulti,multiMax):
    TOTsumPerMulti=[sum(TOTPerMulti[multi]) for multi in range(multiMax)]
    POSPerMulti=[ChNbPerMulti[multi][0] for multi in range(multiMax)] # Should be modified to take TOT max fro example rather than first channel
    xedges=np.arange(self.det_nrows+self.ncols*self.nmods+1)-0.5
    print('#######################################################')
    print(xedges)
    yedges = np.arange(1,300)+0.5
    PHS_per_multi=[]
    for multi in range(multiMax):
        H, xedges, yedges = np.histogram2d(POSPerMulti[multi],1000*TOTsumPerMulti[multi]/13, bins=(xedges, yedges))
        Hsum=sum(H,0)
        PHS_per_multi.append(H.T)
    PHS_all_multi=0
    for multi in range(multiMax):
        PHS_all_multi=PHS_all_multi+PHS_per_multi[multi]
    return PHS_all_multi

def import_QDIV_nomad_list(self, file_numbers, time_file, data_folder, board_number=0):

        print('---- Importing list files from ' + str(file_numbers[0]) + ' to ' + str(file_numbers[-1]))

        # -----------------------------
        # Read time mapping
        # -----------------------------
        time_map = {}

        if not time_file or not os.path.exists(time_file):
            print("---- Warning: No time file provided, all files will be assigned the time of the analysis execution")
            now = datetime.now()
            for i in file_numbers:
                key = f"0{i}"
                time_map[key] = now 

        else:

            print(f"---- Reading time mapping from {time_file}...")
            with open(time_file, "r") as tf:


                for line in tf:
                    line = line.strip()
                    if not line:
                        continue

                    try:
                        fname, tstr = line.split(";")

                        base = os.path.splitext(os.path.basename(fname.strip()))[0]

                        try:
                            time_map[base] = datetime.fromisoformat(tstr.strip())
                        except ValueError:
                            print(f"WARNING: bad time format for {base}, using now()")
                            time_map[base] = datetime.now()

                    except ValueError:
                        continue

  
        # -----------------------------
        # Containers
        # -----------------------------
        all_hist_ab  = []
        all_hist_aa  = []
        all_hist_bb  = []
        all_hist_pos = []
        all_hist_ene = []
        global_time_list = []

  
        # -----------------------------
        # Loop over files
        # -----------------------------

        for i_file in file_numbers:

            progress = (i_file - file_numbers[0]) / (file_numbers[-1] - file_numbers[0] + 1) * 100
            print(f"---- Processing file {i_file} ({progress:.2f}%)...", end="\r", flush=True)

            
            filename = f"{data_folder}/0{i_file}.lst"
            # check if file exists before processing. Otherwise full stop
            if not os.path.exists(filename):
                raise FileNotFoundError(f"Missing file: {filename}")
            else:
                print(f"---- Processing file: {filename}")

            key = os.path.splitext(os.path.basename(filename))[0]

            try:
                # allocate histogram (x_bins, y_bins)
                H_ab  = np.zeros((self.det_ncols, self.n_bins_ene), dtype=np.int32)
                H_aa  = np.zeros((self.det_ncols, self.n_bins_ene), dtype=np.int32)
                H_bb  = np.zeros((self.det_ncols, self.n_bins_ene), dtype=np.int32)
                H_pos = np.zeros((self.det_ncols, self.n_bins_pos), dtype=np.int32)
                H_ene = np.zeros((self.det_ncols,  self.n_bins_pos, self.n_bins_ene), dtype=np.int32)


                # check if file exists before opening
                if not os.path.exists(filename):
                    raise ValueError(f"---- Missing file: {filename}")
                    

                raw = np.fromfile(filename, dtype=np.int32)

                if raw.size % 4 != 0:
                    raise ValueError("---- File corrupted: not multiple of 4 int32 words")

                events = raw.reshape(-1, 4)

                w0 = events[:, 0].view(np.uint32)
                w1 = events[:, 1].view(np.uint32)
                w2 = events[:, 2].view(np.uint32)
                w3 = events[:, 3].view(np.uint32)

                # -----------------------------
                # DECODE buffer 0 (w0)
                # -----------------------------
                # crate    = (w0 >> 28) & 0xF
                board    = (w0 >> 22) & 0x3F
                channel  = (w0 >> 16) & 0x3F
                # rollover = w0 & 0xFFFF

                # -----------------------------
                # timestamp (w1 + rollover)
                # -----------------------------
                # timestamp = w1.astype(np.int64)
                # absolute_time = (rollover.astype(np.int64) << 32) | timestamp

                # -----------------------------
                # A / B / X / Y decoding (w2 + w3)
                # -----------------------------
                A   = w3 & 0xFFFF
                ApB = (w3 >> 16) & 0xFFFF
                # X   = w2 >> 16
                # Y   = w2 & 0xFFFF

                A = A.astype(np.int64)
                ApB = ApB.astype(np.int64)
                channel = channel.astype(np.int64)
                board = board.astype(np.int64)

                n = len(events)

                # run over all events but skip first two
                for i in range(2, n):
                    if i % 100000 == 0:   # adjust frequency
                        print(f"... {i/n:.2%}", end="\r", flush=True)
  
                    if board[i] != board_number:
                        continue  # skip events from other boards if present

                    ch  = 31 - channel[i] if self.channel_inversion else channel[i]
                    aa  = A[i]
                    ab  = ApB[i]
                    bb  = ab - aa

                    if ab > 0:
                        pos = float(aa)/float(ab)*self.pos_max  
                    
                    # compute bins directly
                    if 0 <= ch < self.det_ncols :
                        ab_bin  = int(ab  / self.ene_bin_size)
                        aa_bin  = int(aa  / self.ene_bin_size)
                        bb_bin  = int(bb  / self.ene_bin_size)
                        pos_bin = int(pos / self.pos_bin_size)

                        if 0 <= ab_bin < self.n_bins_ene:
                            H_ab[ch, ab_bin] += 1
                        if 0 <= aa_bin < self.n_bins_ene:
                            H_aa[ch, aa_bin] += 1
                        if 0 <= bb_bin < self.n_bins_ene:
                            H_bb[ch, bb_bin] += 1
                        if 0 <= pos_bin < self.n_bins_pos :
                            H_pos[ch, pos_bin] += 1
                        if 0 <= ab_bin < self.n_bins_ene and 0 <= pos_bin < self.n_bins_pos:
                            H_ene[ch, pos_bin, ab_bin] += 1

                all_hist_ab.append(H_ab.T)    
                all_hist_aa.append(H_aa.T)    
                all_hist_bb.append(H_bb.T)    
                all_hist_pos.append(H_pos.T)  
                all_hist_ene.append(H_ene.T)


                # -----------------------------
                # time handling (datetime)
                # -----------------------------
                if key in time_map:
                    t = time_map[key]
                else:
                    print(f"---- Warning: Missing time for {key}, using now()")
                    t = datetime.now()

                global_time_list.append(t)

            except FileNotFoundError:
                print(f"---- Warning: Missing file: {filename}")
                continue

        # -----------------------------
        # Stack into 3D array
        # -----------------------------

        data_ab  = np.stack(all_hist_ab, axis=0)  # (time, ene, ch)
        data_aa  = np.stack(all_hist_aa, axis=0)  # (time, ene, ch)
        data_bb  = np.stack(all_hist_bb, axis=0)  # (time, ene, ch)
        data_pos = np.stack(all_hist_pos, axis=0) # (time, pos, ch)
        data_ene = np.stack(all_hist_ene, axis=0) # (time, ene, pos, ch)

        # -----------------------------
        # Sort by time
        # -----------------------------
        idx = np.argsort(global_time_list)
        data_ab  = data_ab[idx]
        data_aa  = data_aa[idx]
        data_bb  = data_bb[idx]
        data_pos = data_pos[idx]
        data_ene = data_ene[idx]
        global_time_list = np.array(global_time_list)[idx]

        return data_ab, data_aa, data_bb, data_pos, data_ene, global_time_list 


class PANTHER(ndet_lib.ndet):
    """ Detector class """
    def __init__(self):
        ndet_lib.ndet.__init__(self,nrows=256,ncols=32,nmods=9)
        self.site='ILL'
        self.name='PANTHER'
        self.readout='charge division'
        self.tube_length=2500 # mm
        self.tube_diameter=2.56 # mm
        self.ntubes=self.ncols
        self.det_ntubes=self.det_ncols
        self.detname=self.site+'-'+self.name
        self.aspect_ratio=(self.tube_length/self.nrows)/(self.tube_diameter*self.ncols) # Aspect ratio for a single module
        self.det_aspect_ratio=1/((self.tube_length/self.nrows)/(self.tube_diameter))
        
    def importListMode(self):
        print('### Import a list-mode data file produced by ' + self.detname + ' detector')
        
    def importPHS(self,filename):
        print('### Import a nexus file with PHS data produced by ' + self.detname + ' detector')
        f=importnxs(filename)
        PHSdata=np.squeeze(np.sum(f['/entry0/data/Detector_data'],2)).T
        file_time=f.attrs['file_time'].decode()
        f.close()
        return PHSdata,file_time
  
        
class CSPEC(ndet_lib.ndet):
    """ Detector class """
    def __init__(self):
        ndet_lib.ndet.__init__(self,nrows=1024,ncols=32,nmods=1)
        self.site='ILL'
        self.name='ESS-CSPEC'
        self.readout='charge division'
        self.tube_length=3500 # mm
        self.tube_diameter=25.6 # mm
        self.ntubes=self.ncols
        self.det_ntubes=self.det_ncols
        self.detname=self.site+'-'+self.name
        self.aspect_ratio=(self.tube_length/self.nrows)/(self.tube_diameter*self.ncols) # Aspect ratio for a single module
        self.det_aspect_ratio=((self.tube_length/self.nrows)/(self.tube_diameter))    

        self.n_bins_ene     = 2**8
        self.ene_max        = 2**16

        self.n_bins_pos     = 2**6#10
        self.pos_max        = 2**16

        self.n_bins_Dt      = 2**17 # 19
        self.Dt_max         = 2**21 # 26

        self.coincidence_window = 100 # 1 us

        self.ene_bin_size   = self.ene_max / self.n_bins_ene
        self.pos_bin_size   = self.pos_max / self.n_bins_pos
        self.Dt_bin_size    = self.Dt_max / self.n_bins_Dt

        self.channel_inversion = True # set to True if channel numbering is reversed in data files (tube 0 on the right)
    
    def set_n_bins_pos(self, n_bins_pos):
        self.n_bins_pos = n_bins_pos
        self.pos_bin_size = self.pos_max / self.n_bins_pos

    def importListMode(self, file_numbers, time_file, data_folder, board_number):
        ab, aa, bb, pos, gain, time_list = import_QDIV_nomad_list(self, file_numbers, time_file, data_folder, board_number)
        return ab, aa, bb, pos, gain, time_list
        
    def importPHS(self,filename):
        #print('### Import a nexus file with PHS data produced by ' + self.detname + ' detector')
        f=importnxs(filename)
        PHSdata=np.squeeze(np.sum(f['/entry0/data1/Detector1_data'],2)).T
        file_time=f.attrs['file_time'].decode()
 
        #from datetime import datetime
        # Define the format of the input string
        #format_str = "%Y-%m-%dT%H:%M:%S+%Z"

        # Parse the string to a datetime object
        #file_time = datetime.strptime(file_time_str, format_str)
        #file_time=str(datetime.datetime.today())
        f.close()
        return PHSdata,file_time
    
    def importPOS(self,filename):

        f=importnxs(filename)        
        PHSdata=np.squeeze(np.sum(f['/entry0/data1/Detector1_data'],2)).T
        BM_sum=f['/entry0/monitor1/monitor1data'][0,0]
        file_time=f.attrs['file_time'].decode()
        
        #from datetime import datetime
        # Define the format of the input string
        #format_str = "%Y-%m-%dT%H:%M:%S+%Z"

        # Parse the string to a datetime object
        #file_time = datetime.strptime(file_time_str, format_str)
        #file_time=str(datetime.datetime.today())
        f.close()
        return PHSdata,file_time,BM_sum



class D20(ndet_lib.ndet):
    """ Detector class """
    def __init__(self):
        ndet_lib.ndet.__init__(self,nrows=0,ncols=64,nmods=25)
        self.site='ILL'
        self.name='D20'
        self.readout='indiv. readout'
        #self.tube_length=3500 # mm
        #self.tube_diameter=2.56 # mm
        self.ntubes=self.ncols
        self.det_ntubes=self.det_ncols
        self.detname=self.site+'-'+self.name
        #self.aspect_ratio=(self.tube_length/self.nrows)/(self.tube_diameter*self.ncols) # Aspect ratio for a single module
        #self.det_aspect_ratio=1/((self.tube_length/self.nrows)/(self.tube_diameter))
        
    def import_list_mode(self):
        print('### Import a list-mode data file produced by ' + self.detname + ' detector')
        
    def import_cluster_dump_list_mode(self):
        print('### Import a cluster-dump list-mode data file produced by ' + self.detname + ' detector')
    
    def importPHS(self,filename):
        print('### Import a cluster-dump list-mode data file produced by ' + self.detname + ' detector')
        file_time=str(datetime.datetime.today())
        ToAPerMulti,TOTPerMulti,ChNbPerMulti,sepIdxPerMulti,multiMax=import_MCC_cluster_dump(self,filename)
        PHSdata=createPHS(self,ToAPerMulti,TOTPerMulti,ChNbPerMulti,sepIdxPerMulti,multiMax)
        return PHSdata,file_time
    

class T3(ndet_lib.ndet):
    """ Detector class """
    def __init__(self):
        ndet_lib.ndet.__init__(self,nrows=128,ncols=128,nmods=1)
        self.site='ILL'
        self.name='T3'
        self.readout='indiv. readout'
        #self.tube_length=3500 # mm
        #self.tube_diameter=2.56 # mm
        self.ntubes=self.ncols
        self.det_ntubes=self.det_ncols
        self.detname=self.site+'-'+self.name
        #self.aspect_ratio=(self.tube_length/self.nrows)/(self.tube_diameter*self.ncols) # Aspect ratio for a single module
        #self.det_aspect_ratio=1/((self.tube_length/self.nrows)/(self.tube_diameter))
        
    def import_list_mode(self):
        print('### Import a list-mode data file produced by ' + self.detname + ' detector')
        
    def import_cluster_dump_list_mode(self):
        print('### Import a cluster-dump list-mode data file produced by ' + self.detname + ' detector')

    def importPHS(self,filename):
        print('### Import a cluster-dump list-mode data file produced by ' + self.detname + ' detector')
        file_time=str(datetime.datetime.today())
        ToAPerMulti,TOTPerMulti,ChNbPerMulti,sepIdxPerMulti,multiMax=import_MCC_cluster_dump(self,filename)
        PHSdata=createPHS(self,ToAPerMulti,TOTPerMulti,ChNbPerMulti,sepIdxPerMulti,multiMax)
        return PHSdata,file_time
    