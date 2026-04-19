# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 16:48:26 2026
Definition of a neutron detector class

@author: marchalj
"""


class ndet:
    """ Neutron detector class """
    def __init__(self,nrows,ncols,nmods):
        self.nrows=nrows
        self.ncols=ncols
        self.nmods=nmods
        self.det_nrows=nrows
        self.det_ncols=ncols*nmods


def importnxs(self,filename):
    import h5py
    f=h5py.File(filename,'r')
    return f


