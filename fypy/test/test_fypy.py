# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : test_fypy.py

@Modify Time :  2022/9/23 18:05   

@Author      : Lee    

@Version     : 1.0   

@Description :

'''
import os
import sys
import numpy as np
import datetime

def test_fy4():

    from fypy import fy4pro
    from fypy import readhdf

    l1name = r'D:\DATA\FY4A\L1\20220810\FY4A-_AGRI--_N_DISK_1047E_L1-_FDI-_MULT_NOM_20220810040000_20220810041459_4000M_V0001.HDF'
    geoname =r'D:\DATA\FY4A\L1\20220810\FY4A-_AGRI--_N_DISK_1047E_L1-_GEO-_MULT_NOM_20220810040000_20220810041459_4000M_V0001.HDF'

    fillvalue = 0.0
    data = readhdf(l1name, 'NOMChannel01')
    cal = readhdf(l1name, 'CALChannel01')
    flag = (data<0) | (data>=len(cal))
    data[flag] = 0

    ref = cal[data]
    ref[flag] = fillvalue
    mpro = fy4pro()
    mpro.nom2gll(ref, outname='./data/FY4A-_AGRI--_N_DISK_1047E_L1-_GEO-_MULT_NOM_20220810040000.tif',
                 fillvalue=fillvalue, bbox=(70, 18, 140, 55),)


def test_h8():
    from fypy import hsd2hdf
    outdir = r'D:\DATA\H8\20230115'
    hsdpath = r'D:\DATA\H8'
    nowdate = datetime.datetime.strptime('20230115_0400', '%Y%m%d_%H%M')

    hsd2hdf(outdir, hsdpath, nowdate, SatID='H09')

if __name__ == '__main__':

    # test_fy4()

    test_h8()
