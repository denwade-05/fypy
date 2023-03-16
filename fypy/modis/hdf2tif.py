# coding:utf-8
'''
@Project: fypy
-------------------------------------
@File   : hdf2tif.py
-------------------------------------
@Modify Time      @Author    @Version    
--------------    -------    --------
2021/7/20 17:11     Lee        1.0         
-------------------------------------
@Desciption
-------------------------------------

'''
import glob

from osgeo import gdal, osr
import numpy as np
import os
from scipy.optimize import leastsq


RESAM_GDAL = ['AVERAGE', 'BILINEAR', 'CUBIC', 'CUBIC_SPLINE', 'LANCZOS',
              'MODE', 'NEAREST_NEIGHBOR']
SINU_WKT = 'PROJCS["Sinusoidal_Sanson_Flamsteed",GEOGCS["GCS_Unknown",' \
           'DATUM["D_unknown",SPHEROID["Unknown",6371007.181,"inf"]],' \
           'PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]' \
           ',PROJECTION["Sinusoidal"],PARAMETER["central_meridian",0],' \
           'PARAMETER["false_easting",0],PARAMETER["false_northing",0]' \
           ',UNIT["Meter",1]]'


def getResampling(res):
    """Return the GDAL resampling method

       :param str res: the string of resampling method
    """
    if res == 'AVERAGE':
        return gdal.GRA_Average
    elif res == 'BILINEAR' or res == 'BICUBIC':
        return gdal.GRA_Bilinear
    elif res == 'LANCZOS':
        return gdal.GRA_Lanczos
    elif res == 'MODE':
        return gdal.GRA_Mode
    elif res == 'NEAREST_NEIGHBOR':
        return gdal.GRA_NearestNeighbour
    elif res == 'CUBIC_CONVOLUTION' or res == 'CUBIC':
        return gdal.GRA_Cubic
    elif res == 'CUBIC_SPLINE':
        return gdal.GRA_CubicSpline

class ConverModisByGDAL():

    def __init__(self, outname, hdfname, sdsname, resolution=None, outformat="GTiff",
                 epsg=None, wkt=None, resampl='NEAREST_NEIGHBOR', vrt=False):
        """Function for the initialize the object"""
        # Open source dataset

        self.tempfile = None
        self.outname = outname
        self.sdsname = sdsname
        self.resolution = resolution

        # 设置输出文件投影
        if epsg:
            self.dst_srs = osr.SpatialReference()
            self.dst_srs.ImportFromEPSG(int(epsg))
            self.dst_wkt = self.dst_srs.ExportToWkt()
        elif wkt:
            try:
                f = open(wkt)
                self.dst_wkt = f.read()
                f.close()
            except:
                self.dst_wkt = wkt
        else:
            raise Exception('You have to set one of the following option: '
                            '"epsg", "wkt"')
        # error threshold the same value as gdalwarp
        self.maxerror = 0.125
        self.resampling = getResampling(resampl)

        self.driver = gdal.GetDriverByName(outformat)
        self.vrt = vrt
        if self.driver is None:
            raise Exception('Format driver %s not found, pick a supported '
                            'driver.' % outformat)

        if isinstance(hdfname, list):
            self._MultiFiles(hdfname)
        elif isinstance(hdfname, str):
            self._GetSourceInfo(hdfname)
            self._createWarped(self.src_driver)
            # self._reprojectOne(self.src_driver)
        else:
            raise Exception('Type for subset parameter not supported')

    def _GetSourceInfo(self, filename, ):
        self.in_name = filename
        self.src_ds = gdal.Open(self.in_name)
        layers = self.src_ds.GetSubDatasets()

        # 获取sdsname所在的图层栅格索引
        self.src_raster = self.GetLayer(layers)
        if self.src_raster is None :
            raise Exception('dataset[%s] is not in the %s' %(self.sdsname, filename))
        self.src_driver = gdal.Open(self.src_raster)
        self.src_proj = self.src_driver.GetProjection()
        self.src_trans = self.src_driver.GetGeoTransform()

        # print(self.src_proj)
        # print(self.src_trans)

        self.src_meta = self.src_driver.GetMetadata()
        #  列数
        # self.tileColumns = int(self.src_meta["DATACOLUMNS"])
        #  行数
        # self.tileRows = int(self.src_meta["DATAROWS"])

        self.tiledata = self.src_driver.ReadAsArray()
        # if 'scale_factor' in self.src_meta and 'add_offset' in self.src_meta :
        #     print('scale_factor: ', np.float(self.src_meta['scale_factor']) , ' add_offset: ', np.float(self.src_meta['add_offset']))
        #     self.tiledata = self.tiledata * np.float32(self.src_meta['scale_factor']) + np.float32(self.src_meta['add_offset'])

        datasize = self.tiledata.shape
        if len(datasize) == 2 :
            #  列数
            self.tileColumns = int(datasize[1])
            #  行数
            self.tileRows = int(datasize[0])
        elif len(datasize) == 3 :
            #  列数
            self.tileColumns = int(datasize[2])
            #  行数
            self.tileRows = int(datasize[1])

        # self.src_meta = src_driver.GetMetadata()
        band = self.src_driver.GetRasterBand(1)

        if '_FillValue' in list(self.src_meta.keys()):
            self.data_fill_value = self.src_meta['_FillValue']
        elif band.GetNoDataValue():
            self.data_fill_value = band.GetNoDataValue()
        else:
            self.data_fill_value = None
        self.datatype = band.DataType


    def _MultiFiles(self, filelist):

        alltrans = []
        alldata = []
        tileindex = []

        countfile = len(filelist)
        if countfile == 0 :
            return None
        elif countfile == 1 :
            self._GetSourceInfo(filelist[0])
            self._createWarped(self.src_driver)
            # self._reprojectOne(self.src_driver)
        else:

            for filename in filelist :
                basename = os.path.basename(filename)
                namelist = basename.split('.')
                if len(namelist[2]) != 6:
                    raise Exception('The row and col id [%s] error ! ' %(namelist[2]))

                tileindex.append([int(namelist[2][1:3]), int(namelist[2][4:]),])

                self._GetSourceInfo(filename)
                alldata.append(self.tiledata)
                alltrans.append(self.src_trans)

                # print(np.mean(self.tiledata), np.max(self.tiledata))
                # print(self.src_trans)

            self.tileindex = np.array(tileindex)
            alldata = np.array(alldata)
            alltrans = np.array(alltrans)

            rowindmax =  np.max(self.tileindex[:, 1])
            rowindmin =  np.min(self.tileindex[:, 1])
            colindmax =  np.max(self.tileindex[:, 0])
            colindmin =  np.min(self.tileindex[:, 0])


            rowtile = int(rowindmax - rowindmin + 1)
            coltile = int(colindmax - colindmin + 1)

            xtotal = int(coltile * self.tileColumns)
            ytotal = int(rowtile * self.tileRows)

            dtype = self._GetNumpyType(self.datatype)
            self.srcdata = np.full(shape=(ytotal, xtotal), fill_value=self.data_fill_value, dtype=dtype)

            for i in range(countfile) :
                indx = self.tileindex[i][0]
                indy = self.tileindex[i][1]
                offx = (indx - colindmin) * self.tileColumns
                offy = (indy - rowindmin) * self.tileRows

                self.srcdata[offy:offy+self.tileRows, offx:offx+self.tileColumns] = alldata[i]

            self.src_trans = tuple([np.min(alltrans[:, 0]), np.mean(alltrans[:, 1]), np.mean(alltrans[:, 2]),
                                    np.max(alltrans[:, 3]), np.mean(alltrans[:, 4]), np.mean(alltrans[:, 5])])

            name = os.path.basename(self.outname)
            pathdir = os.path.dirname(self.outname)
            self.tempfile = os.path.join(pathdir, str(hash(name)) + '.tif')
            dr = gdal.GetDriverByName("GTiff")
            self.src_driver = dr.Create(self.tempfile, xtotal, ytotal, 1, self._GetGdalType(dtype))

            self.src_driver.SetProjection(self.src_proj)
            self.src_driver.SetGeoTransform(self.src_trans)
            self.src_driver.GetRasterBand(1).WriteArray(self.srcdata)
            # if self.data_fill_value:
            #     self.src_driver.GetRasterBand(1).SetNoDataValue(float(self.data_fill_value))
            #     self.src_driver.GetRasterBand(1).Fill(float(self.data_fill_value))

            self._createWarped(self.src_driver)

    def __del__(self):
        if self.tempfile :
            tempfile = self.tempfile

            del self.src_driver
            if os.path.isfile(tempfile) :
                os.remove(tempfile)
                tempfile = None


    def _createWarped(self, src_driver):
        '''
        Create a warped VRT file to fetch default values for target raster
        dimensions and geotransform

        :param raster: the name of raster, for HDF have to be one subset
        :return:
        '''

        # 投影转换
        tmp_ds = gdal.AutoCreateWarpedVRT(src_driver, self.src_proj,
                                          self.dst_wkt, self.resampling,
                                          self.maxerror)

        if not self.resolution:
            self.dst_xsize = tmp_ds.RasterXSize
            self.dst_ysize = tmp_ds.RasterYSize
            self.dst_trans = tmp_ds.GetGeoTransform()
        else:
            bbox = self._boundingBox(tmp_ds)
            self.dst_xsize = self._calculateRes(bbox[0][0], bbox[1][0],
                                                self.resolution)
            self.dst_ysize = self._calculateRes(bbox[0][1], bbox[1][1],
                                                self.resolution)
            if self.dst_xsize == 0:
                raise Exception('Invalid number of pixel 0 for X size. The '
                                'problem could be in an invalid value of '
                                'resolution')
            elif self.dst_ysize == 0:
                raise Exception('Invalid number of pixel 0 for Y size. The '
                                'problem could be in an invalid value of '
                                'resolution')
            self.dst_trans = [bbox[0][0], self.resolution, 0.0,
                              bbox[1][1], 0.0, -self.resolution]
        del tmp_ds


    def _boundingBox(self, src):
        """Obtain the bounding box of raster in the new coordinate system

           :param src: a GDAL dataset object

           :return: a bounding box value in lists
        """
        src_gtrn = src.GetGeoTransform(can_return_null=True)

        src_bbox_cells = ((0., 0.),
                          (0, src.RasterYSize),
                          (src.RasterXSize, 0),
                          (src.RasterXSize, src.RasterYSize))

        geo_pts_x = []
        geo_pts_y = []
        for x, y in src_bbox_cells:
            x2 = src_gtrn[0] + src_gtrn[1] * x + src_gtrn[2] * y
            y2 = src_gtrn[3] + src_gtrn[4] * x + src_gtrn[5] * y
            geo_pts_x.append(x2)
            geo_pts_y.append(y2)
        return ((min(geo_pts_x), min(geo_pts_y)), (max(geo_pts_x),
                                                   max(geo_pts_y)))


    def _calculateRes(self, minn, maxx, res):
        """Calculate the number of pixel from extent and resolution

           :param float minn: minimum value of extent
           :param float maxx: maximum value of extent
           :param int res: resolution of output raster

           :return: integer number with the number of pixels
        """
        return int(round((maxx - minn) / res))


    def _progressCallback(self, pct, message, user_data):
        """For the progress status"""
        return 1  # 1 to continue, 0 to stop


    def _GetGdalType(self, datatype):
        if datatype == np.byte :
            return gdal.GDT_Byte
        elif datatype == np.uint16 :
            return gdal.GDT_UInt16
        elif datatype == np.int16 :
            return gdal.GDT_Int16
        elif datatype == np.uint32 :
            return gdal.GDT_UInt32
        elif datatype == np.int32 :
            return gdal.GDT_Int32
        elif datatype == np.float32 :
            return gdal.GDT_Float32
        elif datatype == np.float64 :
            return gdal.GDT_Float64
        else:
            return gdal.GDT_Unknown

    def _GetNumpyType(self, datatype):
        if datatype == gdal.GDT_Byte :
            return np.byte
        elif datatype == gdal.GDT_UInt16 :
            return np.uint16
        elif datatype == gdal.GDT_Int16 :
            return np.int16
        elif datatype == gdal.GDT_UInt32 :
            return np.uint32
        elif datatype == gdal.GDT_Int32 :
            return np.int32
        elif datatype == gdal.GDT_Float32 :
            return np.float32
        elif datatype == gdal.GDT_Float64 :
            return np.float64
        else:
            return None

    def GetLayer(self, layers):
        '''
        获取指定的图层的索引名
        :param layers: tuple
        :return: str
        '''

        if self.sdsname:
            for layer in layers :
                l_name = layer[0].split(':')[-1].replace('"','')
                # print(self.sdsname, l_name)
                if self.sdsname == l_name:
                    return layer[0]

        return None

    def savetif(self, outname):

        out_name = outname

        if self.vrt:
            out_name = "{pref}.tif".format(pref=self.outname)

        # 创建Tiff文件
        try:
            dst_ds = self.driver.Create(out_name,
                                        self.dst_xsize, self.dst_ysize,
                                        1, self.datatype)
        except:
            raise Exception('Not possible to create dataset %s' % out_name)

        dst_ds.SetProjection(self.dst_wkt)
        dst_ds.SetGeoTransform(self.dst_trans)

        if self.data_fill_value:
            dst_ds.GetRasterBand(1).SetNoDataValue(float(self.data_fill_value))
            dst_ds.GetRasterBand(1).Fill(float(self.data_fill_value))
        cbk = self._progressCallback


        # value for last parameter of above self._progressCallback
        cbk_user_data = None
        try:
            gdal.ReprojectImage(self.src_driver, dst_ds, self.src_proj,
                                self.dst_wkt, self.resampling, 0,
                                self.maxerror, cbk, cbk_user_data)
            # if not quiet:
            #     print("Layer {name} reprojected".format(name=l))
        except:
            raise Exception('Not possible to reproject dataset '
                            '{name}'.format(name=self.sdsname))
        dst_ds.SetMetadata(self.src_meta)

        del dst_ds

        print('成功输出【%s】' %(out_name))


def hdf2tif(outname, filename, sdsname, resolution=None,
            format='GTiff', epsg=4326, wkt=None, resampl='NEAREST_NEIGHBOR') :

    try:
        mds = ConverModisByGDAL(outname, filename, sdsname, resolution=resolution,
                                outformat=format, epsg=epsg, wkt=wkt, resampl=resampl )
        mds.savetif(outname)
    except BaseException as e :
        print(e)
        return False

    return True



class Hdf2TifFor5min():

    def __init__(self, outName, filename, sdsname, resolution=0.01, fillvalue=-999.0):
        datasets = gdal.Open(filename)
        #打开子数据集
        layers = datasets.GetSubDatasets()

        layer = self.GetLayer(layers, sdsname)
        #打开ndvi数据
        raster = gdal.Open(layer)
        #获取元数据
        Metadata = datasets.GetMetadata()
        data = raster.ReadAsArray()

        data = data * 0.001 * 10 + 0.0
        data[data<0] = fillvalue

        #  获取四个角的维度
        Latitudes = Metadata["GRINGPOINTLATITUDE.1"]
        #  采用", "进行分割
        LatitudesList = Latitudes.split(", ")
        #  获取四个角的经度
        Longitude = Metadata["GRINGPOINTLONGITUDE.1"]
        #  采用", "进行分割
        LongitudeList = Longitude.split(", ")

        # 图像四个角的地理坐标
        GeoCoordinates = np.zeros((4, 2), dtype = "float32")
        GeoCoordinates[0] = np.array([float(LongitudeList[0]),float(LatitudesList[0])])
        GeoCoordinates[1] = np.array([float(LongitudeList[1]),float(LatitudesList[1])])
        GeoCoordinates[2] = np.array([float(LongitudeList[2]),float(LatitudesList[2])])
        GeoCoordinates[3] = np.array([float(LongitudeList[3]),float(LatitudesList[3])])

        #  列数
        # Columns = float(Metadata["DATACOLUMNS"])
        Rows, Columns = data.shape
        #  行数
        # Rows = float(Metadata["DATAROWS"])
        #  图像四个角的图像坐标
        PixelCoordinates = np.array([[0, 0],
                                     [Columns - 1, 0],
                                     [Columns - 1, Rows - 1],
                                     [0, Rows - 1]], dtype = "float32")

        #  计算仿射变换矩阵
        from scipy.optimize import leastsq
        def func(i):
            Transform0, Transform1, Transform2, Transform3, Transform4, Transform5 = i[0], i[1], i[2], i[3], i[4], i[5]
            return [Transform0 + PixelCoordinates[0][0] * Transform1 + PixelCoordinates[0][1] * Transform2 - GeoCoordinates[0][0],
                    Transform3 + PixelCoordinates[0][0] * Transform4 + PixelCoordinates[0][1] * Transform5 - GeoCoordinates[0][1],
                    Transform0 + PixelCoordinates[1][0] * Transform1 + PixelCoordinates[1][1] * Transform2 - GeoCoordinates[1][0],
                    Transform3 + PixelCoordinates[1][0] * Transform4 + PixelCoordinates[1][1] * Transform5 - GeoCoordinates[1][1],
                    Transform0 + PixelCoordinates[2][0] * Transform1 + PixelCoordinates[2][1] * Transform2 - GeoCoordinates[2][0],
                    Transform3 + PixelCoordinates[2][0] * Transform4 + PixelCoordinates[2][1] * Transform5 - GeoCoordinates[2][1],
                    Transform0 + PixelCoordinates[3][0] * Transform1 + PixelCoordinates[3][1] * Transform2 - GeoCoordinates[3][0],
                    Transform3 + PixelCoordinates[3][0] * Transform4 + PixelCoordinates[3][1] * Transform5 - GeoCoordinates[3][1]]
        #  最小二乘法求解
        GeoTransform = leastsq(func,np.asarray((1,1,1,1,1,1)))

        #  获取数据时间
        # date = Metadata["RANGEBEGINNINGDATE"]

        #  第一个子数据集合,也就是NDVI数据
        # DatasetNDVI = datasets.GetSubDatasets()[0][0]
        # RasterNDVI = gdal.Open(DatasetNDVI)
        # NDVI = ndviRaster.ReadAsArray()
        # print(tuple(GeoTransform[0]), data.shape)

        self.array2raster(outName, GeoTransform[0], data)
        print(outName,"Saved successfully!")


        #命名输出完整路径文件名

        #进行几何校正
        # geoData = gdal.Warp(outName, ndviRaster,
        #                     dstSRS = 'EPSG:4326', format = 'GTiff',
        #                     resampleAlg = gdal.GRA_Bilinear)
        # del geoData



    def GetLayer(self, layers, sdsname):
        '''
        获取指定的图层的索引名
        :param layers: tuple
        :return: str
        '''

        if sdsname:
            for layer in layers :
                l_name = layer[0].split(':')[-1].replace('"','')
                # print(self.sdsname, l_name)
                if sdsname == l_name:
                    return layer[0]

        return None

    #  数组保存为tif
    def array2raster(self, outname, GeoTransform, array, fillvalue=-999.0):
        cols = array.shape[1]  # 矩阵列数
        rows = array.shape[0]  # 矩阵行数
        driver = gdal.GetDriverByName('GTiff')
        outRaster = driver.Create(outname, cols, rows, 1, gdal.GDT_Float32)
        # 括号中两个0表示起始像元的行列号从(0,0)开始
        outRaster.SetGeoTransform(tuple(GeoTransform))
        # 获取数据集第一个波段，是从1开始，不是从0开始
        outband = outRaster.GetRasterBand(1)
        outband.WriteArray(array)
        outRasterSRS = osr.SpatialReference()
        # 代码4326表示WGS84坐标
        outRasterSRS.ImportFromEPSG(4326)
        outRaster.SetProjection(outRasterSRS.ExportToWkt())
        outband.FlushCache()

        if fillvalue is not None:
            # for i in range()
            outRaster.GetRasterBand(1).SetNoDataValue(float(fillvalue))




#
# def hdftotif5min(outName, filename, sdsname, resolution=0.01, fillvalue=-999.0):
#     datasets = gdal.Open(filename)
#     #打开子数据集
#     layers = datasets.GetSubDatasets()
#
#     layer = GetLayer(layers, sdsname)
#     #打开ndvi数据
#     raster = gdal.Open(layer)
#     #获取元数据
#     Metadata = datasets.GetMetadata()
#     data = raster.ReadAsArray()
#
#     data = data * 0.001 * 10 + 0.0
#     data[data<0] = fillvalue
#
#     #  获取四个角的维度
#     Latitudes = Metadata["GRINGPOINTLATITUDE.1"]
#     #  采用", "进行分割
#     LatitudesList = Latitudes.split(", ")
#     #  获取四个角的经度
#     Longitude = Metadata["GRINGPOINTLONGITUDE.1"]
#     #  采用", "进行分割
#     LongitudeList = Longitude.split(", ")
#
#     # 图像四个角的地理坐标
#     GeoCoordinates = np.zeros((4, 2), dtype = "float32")
#     GeoCoordinates[0] = np.array([float(LongitudeList[0]),float(LatitudesList[0])])
#     GeoCoordinates[1] = np.array([float(LongitudeList[1]),float(LatitudesList[1])])
#     GeoCoordinates[2] = np.array([float(LongitudeList[2]),float(LatitudesList[2])])
#     GeoCoordinates[3] = np.array([float(LongitudeList[3]),float(LatitudesList[3])])
#
#     #  列数
#     # Columns = float(Metadata["DATACOLUMNS"])
#     Rows, Columns = data.shape
#     #  行数
#     # Rows = float(Metadata["DATAROWS"])
#     #  图像四个角的图像坐标
#     PixelCoordinates = np.array([[0, 0],
#                                  [Columns - 1, 0],
#                                  [Columns - 1, Rows - 1],
#                                  [0, Rows - 1]], dtype = "float32")
#
#     #  计算仿射变换矩阵
#     from scipy.optimize import leastsq
#     def func(i):
#         Transform0, Transform1, Transform2, Transform3, Transform4, Transform5 = i[0], i[1], i[2], i[3], i[4], i[5]
#         return [Transform0 + PixelCoordinates[0][0] * Transform1 + PixelCoordinates[0][1] * Transform2 - GeoCoordinates[0][0],
#                 Transform3 + PixelCoordinates[0][0] * Transform4 + PixelCoordinates[0][1] * Transform5 - GeoCoordinates[0][1],
#                 Transform0 + PixelCoordinates[1][0] * Transform1 + PixelCoordinates[1][1] * Transform2 - GeoCoordinates[1][0],
#                 Transform3 + PixelCoordinates[1][0] * Transform4 + PixelCoordinates[1][1] * Transform5 - GeoCoordinates[1][1],
#                 Transform0 + PixelCoordinates[2][0] * Transform1 + PixelCoordinates[2][1] * Transform2 - GeoCoordinates[2][0],
#                 Transform3 + PixelCoordinates[2][0] * Transform4 + PixelCoordinates[2][1] * Transform5 - GeoCoordinates[2][1],
#                 Transform0 + PixelCoordinates[3][0] * Transform1 + PixelCoordinates[3][1] * Transform2 - GeoCoordinates[3][0],
#                 Transform3 + PixelCoordinates[3][0] * Transform4 + PixelCoordinates[3][1] * Transform5 - GeoCoordinates[3][1]]
#     #  最小二乘法求解
#     GeoTransform = leastsq(func,np.asarray((1,1,1,1,1,1)))
#
#     #  获取数据时间
#     # date = Metadata["RANGEBEGINNINGDATE"]
#
#     #  第一个子数据集合,也就是NDVI数据
#     # DatasetNDVI = datasets.GetSubDatasets()[0][0]
#     # RasterNDVI = gdal.Open(DatasetNDVI)
#     # NDVI = ndviRaster.ReadAsArray()
#     # print(tuple(GeoTransform[0]), data.shape)
#
#     array2raster(outName, GeoTransform[0], data)
#     print(outName,"Saved successfully!")
#
#
#     #命名输出完整路径文件名
#
#     #进行几何校正
#     # geoData = gdal.Warp(outName, ndviRaster,
#     #                     dstSRS = 'EPSG:4326', format = 'GTiff',
#     #                     resampleAlg = gdal.GRA_Bilinear)
#     # del geoData
#
#     return True
#
# def GetLayer(layers, sdsname):
#     '''
#     获取指定的图层的索引名
#     :param layers: tuple
#     :return: str
#     '''
#
#     if sdsname:
#         for layer in layers :
#             l_name = layer[0].split(':')[-1].replace('"','')
#             # print(self.sdsname, l_name)
#             if sdsname == l_name:
#                 return layer[0]
#
#     return None
#
#
#
# #  数组保存为tif
# def array2raster(TifName, GeoTransform, array, fillvalue=-999.0):
#     cols = array.shape[1]  # 矩阵列数
#     rows = array.shape[0]  # 矩阵行数
#     driver = gdal.GetDriverByName('GTiff')
#     outRaster = driver.Create(TifName, cols, rows, 1, gdal.GDT_Float32)
#     # 括号中两个0表示起始像元的行列号从(0,0)开始
#     outRaster.SetGeoTransform(tuple(GeoTransform))
#     # 获取数据集第一个波段，是从1开始，不是从0开始
#     outband = outRaster.GetRasterBand(1)
#     outband.WriteArray(array)
#     outRasterSRS = osr.SpatialReference()
#     # 代码4326表示WGS84坐标
#     outRasterSRS.ImportFromEPSG(4326)
#     outRaster.SetProjection(outRasterSRS.ExportToWkt())
#     outband.FlushCache()
#
#     if fillvalue:
#         outRaster.GetRasterBand(1).SetNoDataValue(float(fillvalue))
#
#
#
#
#
# def mosaicbyshp():
#
#     raster = r'Z:\Data\jilinnongqi\MODIS\NDVI\16D\MYD13A2.A2021121.tif'
#     shp = r'D:\jilin\Arspy\MODIS_hdf2tif\parm\jilincity.shp'
#     # checkshp(shp)
#     # exit()
#
#     outname = r'Z:\Data\jilinnongqi\MODIS\NDVI\jilin.tif'
#     import shapefile
#     r = shapefile.Reader(shp)
#
#     gdal.Warp(outname, raster,
#               format='GTiff',
#               cutlineDSName =shp,
#               # cutlineDSName=None,
#               # cutlineWhere=None,
#               # cutlineSQL=None,
#               # cutlineBlend=None,
#               cropToCutline=True,
#               # dstSRS='EPSG:4326',
#               dstNodata=-3000,
#               # outputBounds=r.bbox,
#               )
#
# def checkshp(shpFile):
#
#     from osgeo import ogr
#
#     # shpFile = 'F:/m2.shp'  # 裁剪矩形
#
#     # # # 注册所有的驱动
#     ogr.RegisterAll()
#
#     # 打开数据
#     ds = ogr.Open(shpFile, 0)
#     if ds is None:
#         print("打开文件【%s】失败！", shpFile)
#         return
#     print("打开文件【%s】成功！", shpFile)
#     # 获取该数据源中的图层个数，一般shp数据图层只有一个，如果是mdb、dxf等图层就会有多个
#     m_layer_count = ds.GetLayerCount()
#     m_layer = ds.GetLayerByIndex(0)
#     if m_layer is None:
#         print("获取第%d个图层失败！\n", 0)
#         return
#     # 对图层进行初始化，如果对图层进行了过滤操作，执行这句后，之前的过滤全部清空
#     m_layer.ResetReading()
#     count = 0
#     m_feature = m_layer.GetNextFeature()
#     while m_feature is not None:
#         o_geometry = m_feature.GetGeometryRef()
#         if not ogr.Geometry.IsValid(o_geometry):
#             print(m_feature.GetFID())
#             count = count + 1
#
#         m_feature = m_layer.GetNextFeature()
#     print("无效多边形共" + str(count) + "个")
#
#
#
# if __name__ == '__main__':
#
#     # mosaicbyshp()
#     # exit()
#
#     pathin = r"Z:\Data\jilinnongqi\MODIS\CLM"
#
#     #  获取文件夹内的文件名
#     hdfNameList = glob.glob(os.path.join(pathin, 'MOD35_L2.A2021138*.hdf'))
#     #
#     # outname = os.path.join(pathin, 'MYD13A2.A2021121.tif')
#     # hdf2tif(outname, hdfNameList, sdsname='1 km 16 days NDVI')
#
#     # exit()
#     for filename in hdfNameList:
#         #  判断当前文件是否为HDF文件
#         hdf2tif(filename+'.tif', filename, sdsname='Cloud_Mask')
#

