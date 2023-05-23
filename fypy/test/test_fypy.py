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

def test_fy3():
    from fypy.fy3 import FY3Orbit, fy3pro
    from fypy.tools import readhdf

    l1file = r'D:\DATA\FY3D\MERSI\FY3D_MERSI_GBAL_L1_20220718_0725_1000M_MS.HDF'
    geofile = r'D:\DATA\FY3D\MERSI\FY3D_MERSI_GBAL_L1_20220718_0725_GEO1K_MS.HDF'

    mpro = fy3pro()
    band14 = mpro.Calibration(l1file, '/Data/EV_250_Aggr.1KM_RefSB')
    band519 = mpro.Calibration(l1file, '/Data/EV_1KM_RefSB')
    # band2023 = mpro.Calibration(l1file, '/Data/EV_1KM_Emissive')
    # band2425 = mpro.Calibration(l1file, '/Data/EV_250_Aggr.1KM_Emissive')

    # data = []
    # data.append(band14)
    # data.append(band519)
    from osgeo import gdal
    data = band14
    # data = np.vstack((band14, band519))
    # data = np.array(data, dtype=np.uint8)
    # data = readhdf(l1file, '/Data/EV_1KM_RefSB')
    lat = readhdf(geofile, '/Geolocation/Latitude')
    lon = readhdf(geofile, '/Geolocation/Longitude')
    mpro = FY3Orbit(data, lat, lon, dstfile=r'D:\DATA\FY3D\MERSI\test100-1.tif',
                    resolution=0.0025, vmin=0, vmax=10000, resampleAlg=gdal.GRIORA_Bilinear)

if __name__ == '__main__':

    # test_fy4()

    test_fy3()