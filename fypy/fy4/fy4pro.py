# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : fy4pro.py

@Modify Time :  2022/11/10 14:45

@Author      : Lee

@Version     : 1.0

@Description :

'''
import os
import sys
import numpy as np
from osgeo import gdal, osr
from fypy.tools import readhdf_fileinfo

class fy4pro(object) :

    def __init__(self):
        pass


    def nom2gll(self, data, outname=None, bbox=None, subpoint=104.7, resolution=0.04,
                fillvalue=None, dstSRS='EPSG:4326'):
        '''
        将DISK 投影转换成等经纬投影

        Parameters
        ----------
        data : numpy.array
            输入数据
        outname: str, optional
            输出文件名
        bbox : tuple, optional
            output bounds as (minX, minY, maxX, maxY) in target SRS
        subpoint: float
            星下点
        resolution: float
            数据图像的分辨率，单位为degree
        fillvalue : float
            数据填充值
        dstSRS

        Returns
        -------

        '''
        im_data = np.array(data, dtype=np.float32)
        dtype = self._gettype(im_data.dtype)
        if len(im_data.shape) == 3:
            im_bands, im_height, im_width = im_data.shape
        elif len(im_data.shape) == 2:
            im_bands, (im_height, im_width) = 1,im_data.shape
        else:
            im_bands, (im_height, im_width) = 1,im_data.shape

        Driver = gdal.GetDriverByName('MEM')
        memDs = Driver.Create('', im_height, im_width, im_bands, dtype)

        # 写入数据
        if im_bands == 1:
            memDs.GetRasterBand(1).WriteArray(im_data)
        else:
            for i in range(im_bands):
                memDs.GetRasterBand(i+1).WriteArray(im_data[i])

        srs = osr.SpatialReference()
        srs.ImportFromProj4('+proj=geos +h=35785863 +a=6378137.0 +b=6356752.3 +lon_0={subpoint} +no_defs'.format(
            subpoint=subpoint ))
        memDs.SetProjection(srs.ExportToWkt())
        memDs.SetGeoTransform([-5496000, resolution*100*1000, 0, 5496000, 0, -resolution*100*1000])
        if outname is None :
            warpDs = gdal.Warp('', memDs, format='MEM', dstSRS=dstSRS,
                               outputBounds=bbox, xRes=resolution, yRes=resolution,
                               dstNodata=fillvalue)
        else:
            warpDs = gdal.Warp(outname, memDs, format='GTiff', dstSRS=dstSRS,
                               outputBounds=bbox, xRes=resolution, yRes=resolution,
                               srcNodata=fillvalue, dstNodata=fillvalue, creationOptions=["COMPRESS=LZW"])
        if warpDs is None :
            return None

        data_GLL = warpDs.ReadAsArray(0, 0, warpDs.RasterXSize, warpDs.RasterYSize)
        del warpDs

        return data_GLL

    def getL1Data(self, filename, bandID=1):
        '''
        读取FY4 L1数据，并完成辐射定标
        Parameters
        ----------
        filename: str
            L1数据文件名
        bandID : int
            波段索引，从1开始

        Returns
        -------
            numpy.array
            辐射定标转换后的ref或bt
        '''
        if not os.path.isfile(filename) :
            print('文件不存在【%s】' %(filename))
            return None

        import h5py

        # 转换到区域的行列号（考虑去除图像偏移）
        fileinfo = readhdf_fileinfo(filename)
        Begin_Line_Number = fileinfo['Begin Line Number'][0]
        End_Line_Number = fileinfo['End Line Number'][0]
        Begin_Pixel_Number = fileinfo['Begin Pixel Number'][0]
        End_Pixel_Number = fileinfo['End Pixel Number'][0]

        data = np.zeros(shape=(End_Line_Number+1, End_Pixel_Number+1),dtype=np.float32)

        fp = h5py.File(filename, 'r')
        cal = fp['CALChannel%02d' %(bandID)][:]
        dn = fp['NOMChannel%02d' %(bandID)][:]
        fp.close()

        dn[dn>=len(cal)] = 0

        data[Begin_Line_Number:End_Line_Number+1, Begin_Pixel_Number:End_Pixel_Number+1] = cal[dn]

        return data


    def getGEOData(self, filename, sdsname):
        import h5py
        fp = h5py.File(filename, 'r')
        data1 = fp[sdsname][:]
        fp.close()

        # 转换到区域的行列号（考虑去除图像偏移）
        fileinfo = readhdf_fileinfo(filename)
        Begin_Line_Number = fileinfo['Begin Line Number'][0]
        End_Line_Number = fileinfo['End Line Number'][0]
        Begin_Pixel_Number = fileinfo['Begin Pixel Number'][0]
        End_Pixel_Number = fileinfo['End Pixel Number'][0]

        data = np.zeros(shape=(End_Line_Number+1, End_Pixel_Number+1),dtype=np.float32)

        data[Begin_Line_Number:End_Line_Number+1, Begin_Pixel_Number:End_Pixel_Number+1] = data1

        return data

    def _gettype(self, datatype):
        ''' 根据numpy的数据类型，匹配GDAL中的数据类型 '''

        if datatype == np.byte or datatype == np.uint8:
            return gdal.GDT_Byte
        elif datatype == np.uint16 :
            return gdal.GDT_UInt16
        elif datatype == np.int16 :
            return gdal.GDT_Int16
        elif datatype == np.uint32 :
            return gdal.GDT_UInt32
        elif datatype == np.int32 :
            return gdal.GDT_Int32
        elif datatype == np.float32 or datatype.str in ['>f4', '<f4']:
            return gdal.GDT_Float32
        elif datatype == np.float64 or datatype.str in ['>f8', '<f8']:
            return gdal.GDT_Float64
        else:
            return gdal.GDT_Unknown


