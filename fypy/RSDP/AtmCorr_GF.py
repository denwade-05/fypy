# -*- coding:utf-8 -*-
'''
@Project     : lb_toolkits

@File        : AtmCorr_GF.py

@Modify Time :  2022/12/23 16:38   

@Author      : Lee    

@Version     : 1.0   

@Description :

'''
import os
import shutil
import sys
import numpy as np
import datetime
import json
from tqdm import tqdm
from AtmCorr import AtmCorr
from osgeo import gdal
from dateutil.relativedelta import relativedelta
import xml.dom.minidom


class AtmCorr_GF(AtmCorr):

    def __init__(self):
        super(AtmCorr_GF, self).__init__()

    def FLAASH(self, outname, tiffFile, metedata, nowdate, SatID, InstID, InstType=None, fillvalue=65535):

        outdir = os.path.dirname(outname)
        if not os.path.isdir(outdir) :
            os.makedirs(outdir)
            print('成功创建路径【%s】' %(outdir))

        if not os.path.isfile(tiffFile) :
            print('文件不存在【%s】' %(tiffFile))
            return None

        srcdataset = gdal.Open(tiffFile, gdal.GA_ReadOnly)
        if srcdataset == None:
            print("文件无法打开【%s】" %(tiffFile))
            return None

        cols = srcdataset.RasterXSize                  # 栅格矩阵的列数
        rows = srcdataset.RasterYSize                  # 栅格矩阵的行数
        Bands = srcdataset.RasterCount                 # 波段数

        geoTransform1 = srcdataset.GetGeoTransform()   # 获取仿射矩阵信息
        proj1 = srcdataset.GetProjection()             # 获取投影信息

        # 创建输出结果文件
        driver = gdal.GetDriverByName("GTiff")
        dstdataset = driver.Create(outname, cols, rows, 4, gdal.GDT_UInt16,
                                   options=["COMPRESS=LZW", "BigTIFF=YES"])
        dstdataset.SetGeoTransform(geoTransform1)
        dstdataset.SetProjection(proj1)

        #分别读取4个波段
        for BandID in range(0, Bands):
            ReadBand = srcdataset.GetRasterBand(BandID+1)
            outband = dstdataset.GetRasterBand(BandID+1)
            outband.SetNoDataValue(fillvalue)
            #获取对应波段的增益gain和偏移bias
            Gain, Bias = self.RadiometricCalibration(nowdate, SatID, InstID)

            if 'PMS' in InstID :
                if InstType is None :
                    raise Exception('请传入InstType参数【全色：PAN，多光谱：MSS】')

                if 'MSS' in InstType :
                    BandID += 1
                elif 'PAN' in InstType :
                    BandID = BandID
                else:
                    raise Exception('请传入InstType参数【全色：PAN，多光谱：MSS】')

            self.setParam(nowdate, metedata, SatID, InstID, BandID+1)
            #获取大气校正系数
            AtcCofa, AtcCofb, AtcCofc = self.corrCoeff()
            nBlockSize = 2048
            i = 0
            j = 0
            #进度条参数
            XBlockcount = np.ceil(cols/nBlockSize)
            YBlockcount = np.ceil(rows/nBlockSize)
            # print("第%d波段校正："%BandID)
            try:
                with tqdm(total=XBlockcount*YBlockcount, iterable='iterable',
                          desc = '正在进行第%i波段校正'%(BandID+1), mininterval=1) as pbar:
                    while i < rows:
                        while j < cols:
                            # 保存分块大小
                            nXBK = nBlockSize
                            nYBK = nBlockSize

                            # 最后一块
                            if i+nBlockSize>rows:
                                nYBK = rows - i
                            if j+nBlockSize>cols:
                                nXBK=cols - j

                            # 分块读取影像
                            Image = ReadBand.ReadAsArray(j, i, nXBK,nYBK)

                            outImage =np.where(Image>0, Image*Gain[BandID] + Bias[BandID], fillvalue)

                            y = np.where(outImage!=fillvalue, AtcCofa * outImage - AtcCofb, fillvalue)
                            atcImage = np.where(y!=fillvalue, (y / (1 + y * AtcCofc))*10000, fillvalue)

                            outband.WriteArray(atcImage,j,i)
                            j=j+nXBK
                            pbar.update(1)
                        j=0
                        i=i+nYBK
            except KeyboardInterrupt:
                pbar.close()
                raise
            pbar.close()

        dstdataset = None
        srcdataset = None

    def setParam(self, nowdate, metedata, SatID, InstID, BandId):

        # 读取头文件
        dom = xml.dom.minidom.parse(metedata)

        # 太阳和卫星角度信息
        sunz = 90 - float(dom.getElementsByTagName('SolarZenith')[0].firstChild.data)
        suna = float(dom.getElementsByTagName('SolarAzimuth')[0].firstChild.data)
        # s.geometry.view_z = float(dom.getElementsByTagName('SatelliteZenith')[0].firstChild.data)
        # s.geometry.view_a = float(dom.getElementsByTagName('SatelliteAzimuth')[0].firstChild.data)
        satz = 0
        sata = 0
        # 日期
        DateTimeparm = dom.getElementsByTagName('CenterTime')[0].firstChild.data
        DateTime = datetime.datetime.strptime(DateTimeparm, '%Y-%m-%d %H:%M:%S')

        # 中心经纬度
        TopLeftLat     = float(dom.getElementsByTagName('TopLeftLatitude')[0].firstChild.data)
        TopLeftLon     = float(dom.getElementsByTagName('TopLeftLongitude')[0].firstChild.data)
        TopRightLat    = float(dom.getElementsByTagName('TopRightLatitude')[0].firstChild.data)
        TopRightLon    = float(dom.getElementsByTagName('TopRightLongitude')[0].firstChild.data)
        BottomRightLat = float(dom.getElementsByTagName('BottomRightLatitude')[0].firstChild.data)
        BottomRightLon = float(dom.getElementsByTagName('BottomRightLongitude')[0].firstChild.data)
        BottomLeftLat  = float(dom.getElementsByTagName('BottomLeftLatitude')[0].firstChild.data)
        BottomLeftLon  = float(dom.getElementsByTagName('BottomLeftLongitude')[0].firstChild.data)

        ImageCenterLat = (TopLeftLat + TopRightLat + BottomRightLat + BottomLeftLat) / 4

        self.set_geom(nowdate, sunz=sunz, suna=suna, satz=satz, sata=sata)
        self.set_atm(sLatitude=ImageCenterLat, nowdate=nowdate)
        self.set_aer()
        self.set_vis()
        self.set_altitude(dem=0.010)

        script_path = os.path.split(os.path.realpath(__file__))[0]
        #读取辐射校正和大气校正所需参数:增益、偏移和光谱响应函数
        SRFFile = os.path.join(script_path, 'resp', 'GF', "GF.json")
        if not os.path.isfile(SRFFile) :
            raise Exception('高分卫星光谱响应文件不存在【%s】' %(SRFFile))
        SRF = json.load(open(SRFFile))

        minwl = SRF[SatID][InstID]['B%d' %(BandId)]['wl'][0]
        maxwl = SRF[SatID][InstID]['B%d' %(BandId)]['wl'][1]
        response = SRF[SatID][InstID]['B%d' %(BandId)]['SRF']
        self.set_resp(minwl, maxwl, response)

    def RadiometricCalibration(self, nowdate, SatID, InstID):
        ''' 获取高分卫星定标系数 gain和 bias '''

        script_path = os.path.split(os.path.realpath(__file__))[0]
        #读取辐射校正和大气校正所需参数:增益、偏移和光谱响应函数
        CalCoeffFile = os.path.join(script_path, "CalibrationCoefficient.json")
        if not os.path.isfile(CalCoeffFile) :
            raise Exception('高分卫星定标系数文件不存在【%s】' %(CalCoeffFile))
        CalCoeff = json.load(open(CalCoeffFile))

        satid = SatID.replace('-', '')
        instid = InstID.replace('-', '')

        if not satid in CalCoeff :
            raise Exception('请确认【%s】是否在定标系数列表，当前仅支持' %(SatID), CalCoeff.keys())

        if not instid in CalCoeff[satid] :
            raise Exception('请确认【%s】是否在定标系数列表，当前仅支持' %(InstID), CalCoeff[satid].keys())

        while nowdate >= datetime.datetime.strptime('2000', '%Y') :
            stryear = nowdate.strftime('%Y')
            if not stryear in CalCoeff[satid][instid] :
                nowdate -= relativedelta(years=1)
            else:
                gain = CalCoeff[satid][instid][stryear]['gain']
                bias = CalCoeff[satid][instid][stryear]['bias']
                break

        return gain, bias

    def ortho_rectification(self, dstfile, srcfile, rpcfile, demfile=None):
        """
        正射校正
        :param input:输入原始影像
        :param output:输出正射影像
        """

        if not os.path.isfile(rpcfile) :
            raise Exception('RPC文件不存在【%s】' %(rpcfile))

        srcdir = os.path.dirname(srcfile)
        shutil.copy(rpcfile, srcdir)

        dataset = gdal.Open(srcfile, gdal.GA_Update)#读入影像
        rpc = dataset.GetMetadata("RPC")#读入影像，rpc

        if demfile is None :
            dst_ds = gdal.Warp(dstfile, dataset, dstSRS='EPSG:4326',
                               # xRes=resolution,
                               # yRes=resolution,
                               resampleAlg=gdal.GRIORA_Bilinear,
                               rpc=True, #使用RPC模型进行校正
                               warpOptions=['INIT_DEST=NO_DATA'],
                               creationOptions=["COMPRESS=LZW", "BigTIFF=YES"])
        else:
            dst_ds = gdal.Warp(dstfile, dataset, dstSRS='EPSG:4326',
                               # xRes=resolution,
                               # yRes=resolution,
                               # resampleAlg=gdal.GRIORA_Bilinear,
                               rpc=True, #使用RPC模型进行校正
                               transformerOptions=[r'RPC_DEM=%s' %(demfile),
                                                   "RPC_DEMINTERPOLATION=bilinear"],
                               creationOptions=["COMPRESS=LZW", "BigTIFF=YES"])

        dst_ds = None

    def DoesDriverHandleExtension(self, drv, ext):
        exts = drv.GetMetadataItem(gdal.DMD_EXTENSIONS)
        return exts is not None and exts.lower().find(ext.lower()) >= 0


    def GetExtension(self, filename):
        ext = os.path.splitext(filename)[1]
        if ext.startswith('.'):
            ext = ext[1:]
        return ext


    def GetOutputDriversFor(self, filename):
        drv_list = []
        ext = self.GetExtension(filename)
        for i in range(gdal.GetDriverCount()):
            drv = gdal.GetDriver(i)
            if (drv.GetMetadataItem(gdal.DCAP_CREATE) is not None or
                drv.GetMetadataItem(gdal.DCAP_CREATECOPY) is not None) and \
                    drv.GetMetadataItem(gdal.DCAP_RASTER) is not None:
                if ext and self.DoesDriverHandleExtension(drv, ext):
                    drv_list.append(drv.ShortName)
                else:
                    prefix = drv.GetMetadataItem(gdal.DMD_CONNECTION_PREFIX)
                    if prefix is not None and filename.lower().startswith(prefix.lower()):
                        drv_list.append(drv.ShortName)

        # GMT is registered before netCDF for opening reasons, but we want
        # netCDF to be used by default for output.
        if ext.lower() == 'nc' and not drv_list and \
                drv_list[0].upper() == 'GMT' and drv_list[1].upper() == 'NETCDF':
            drv_list = ['NETCDF', 'GMT']

        return drv_list


    def GetOutputDriverFor(self, filename):
        drv_list = self.GetOutputDriversFor(filename)
        ext = self.GetExtension(filename)
        if not drv_list:
            if not ext:
                return 'GTiff'
            else:
                raise Exception("Cannot guess driver for %s" % filename)
        elif len(drv_list) > 1:
            print("Several drivers matching %s extension. Using %s" % (ext if ext else '', drv_list[0]))
        return drv_list[0]


    def Usage(self):
        print('Usage: gdal_pansharpen [--help-general] pan_dataset {spectral_dataset[,band=num]}+ out_dataset')
        print('                       [-of format] [-b band]* [-w weight]*')
        print('                       [-r {nearest,bilinear,cubic,cubicspline,lanczos,average}]')
        print('                       [-threads {ALL_CPUS|number}] [-bitdepth val] [-nodata val]')
        print('                       [-spat_adjust {union,intersection,none,nonewithoutwarning}]')
        print('                       [-verbose_vrt] [-co NAME=VALUE]* [-q]')
        print('')
        print('Create a dataset resulting from a pansharpening operation.')
        return -1


    def gdal_pansharpen(self, argv):

        argv = gdal.GeneralCmdLineProcessor(argv)
        if argv is None:
            return -1

        pan_name = None
        last_name = None
        spectral_ds = []
        spectral_bands = []
        out_name = None
        bands = []
        weights = []
        frmt = None
        creation_options = []
        callback = gdal.TermProgress_nocb
        resampling = None
        spat_adjust = None
        verbose_vrt = False
        num_threads = None
        bitdepth = None
        nodata = None

        i = 1
        argc = len(argv)
        while i < argc:
            if (argv[i] == '-of' or argv[i] == '-f') and i < len(argv) - 1:
                frmt = argv[i + 1]
                i = i + 1
            elif argv[i] == '-r' and i < len(argv) - 1:
                resampling = argv[i + 1]
                i = i + 1
            elif argv[i] == '-spat_adjust' and i < len(argv) - 1:
                spat_adjust = argv[i + 1]
                i = i + 1
            elif argv[i] == '-b' and i < len(argv) - 1:
                bands.append(int(argv[i + 1]))
                i = i + 1
            elif argv[i] == '-w' and i < len(argv) - 1:
                weights.append(float(argv[i + 1]))
                i = i + 1
            elif argv[i] == '-co' and i < len(argv) - 1:
                creation_options.append(argv[i + 1])
                i = i + 1
            elif argv[i] == '-threads' and i < len(argv) - 1:
                num_threads = argv[i + 1]
                i = i + 1
            elif argv[i] == '-bitdepth' and i < len(argv) - 1:
                bitdepth = argv[i + 1]
                i = i + 1
            elif argv[i] == '-nodata' and i < len(argv) - 1:
                nodata = argv[i + 1]
                i = i + 1
            elif argv[i] == '-q':
                callback = None
            elif argv[i] == '-verbose_vrt':
                verbose_vrt = True
            elif argv[i][0] == '-':
                sys.stderr.write('Unrecognized option : %s\n' % argv[i])
                return self.Usage()
            elif pan_name is None:
                pan_name = argv[i]
                pan_ds = gdal.Open(pan_name)
                if pan_ds is None:
                    return 1
            else:
                # print(last_name)
                if last_name is not None:
                    pos = last_name.find(',band=')
                    # print(last_name)
                    if pos > 0:
                        spectral_name = last_name[0:pos]
                        ds = gdal.Open(spectral_name)
                        if ds is None:
                            return 1
                        band_num = int(last_name[pos + len(',band='):])
                        band = ds.GetRasterBand(band_num)
                        spectral_ds.append(ds)
                        spectral_bands.append(band)
                    else:
                        spectral_name = last_name
                        ds = gdal.Open(spectral_name)
                        if ds is None:
                            return 1
                        for j in range(ds.RasterCount):
                            spectral_ds.append(ds)
                            spectral_bands.append(ds.GetRasterBand(j + 1))

                last_name = argv[i]
                # print(last_name)

            i = i + 1

        # print(spectral_name)
        if pan_name is None or not spectral_bands:
            return self.Usage()
        out_name = last_name
        # print(out_name)
        if frmt is None:
            frmt = self.GetOutputDriverFor(out_name)

        if not bands:
            bands = [j + 1 for j in range(len(spectral_bands))]
        else:
            for band in bands:
                if band < 0 or band > len(spectral_bands):
                    print('Invalid band number in -b: %d' % band)
                    return 1

        if weights and len(weights) != len(spectral_bands):
            print('There must be as many -w values specified as input spectral bands')
            return 1

        vrt_xml = """<VRTDataset subClass="VRTPansharpenedDataset">\n"""
        if bands != [j + 1 for j in range(len(spectral_bands))]:
            for i, band in enumerate(bands):
                sband = spectral_bands[band - 1]
                datatype = gdal.GetDataTypeName(sband.DataType)
                colorname = gdal.GetColorInterpretationName(sband.GetColorInterpretation())
                vrt_xml += """  <VRTRasterBand dataType="%s" band="%d" subClass="VRTPansharpenedRasterBand">
          <ColorInterp>%s</ColorInterp>
      </VRTRasterBand>\n""" % (datatype, i + 1, colorname)

        vrt_xml += """  <PansharpeningOptions>\n"""

        if weights:
            vrt_xml += """      <AlgorithmOptions>\n"""
            vrt_xml += """        <Weights>"""
            for i, weight in enumerate(weights):
                if i > 0:
                    vrt_xml += ","
                vrt_xml += "%.16g" % weight
            vrt_xml += "</Weights>\n"
            vrt_xml += """      </AlgorithmOptions>\n"""

        if resampling is not None:
            vrt_xml += '      <Resampling>%s</Resampling>\n' % resampling

        if num_threads is not None:
            vrt_xml += '      <NumThreads>%s</NumThreads>\n' % num_threads

        if bitdepth is not None:
            vrt_xml += '      <BitDepth>%s</BitDepth>\n' % bitdepth

        if nodata is not None:
            vrt_xml += '      <NoData>%s</NoData>\n' % nodata

        if spat_adjust is not None:
            vrt_xml += '      <SpatialExtentAdjustment>%s</SpatialExtentAdjustment>\n' % spat_adjust

        pan_relative = '0'
        if frmt.upper() == 'VRT':
            if not os.path.isabs(pan_name):
                pan_relative = '1'
                pan_name = os.path.relpath(pan_name, os.path.dirname(out_name))

        vrt_xml += """    <PanchroBand>
          <SourceFilename relativeToVRT="%s">%s</SourceFilename>
          <SourceBand>1</SourceBand>
        </PanchroBand>\n""" % (pan_relative, pan_name)

        for i, sband in enumerate(spectral_bands):
            dstband = ''
            for j, band in enumerate(bands):
                if i + 1 == band:
                    dstband = ' dstBand="%d"' % (j + 1)
                    break

            ms_relative = '0'
            ms_name = spectral_ds[i].GetDescription()
            if frmt.upper() == 'VRT':
                if not os.path.isabs(ms_name):
                    ms_relative = '1'
                    ms_name = os.path.relpath(ms_name, os.path.dirname(out_name))

            vrt_xml += """    <SpectralBand%s>
          <SourceFilename relativeToVRT="%s">%s</SourceFilename>
          <SourceBand>%d</SourceBand>
        </SpectralBand>\n""" % (dstband, ms_relative, ms_name, sband.GetBand())

        vrt_xml += """  </PansharpeningOptions>\n"""
        vrt_xml += """</VRTDataset>\n"""

        if frmt.upper() == 'VRT':
            f = gdal.VSIFOpenL(out_name, 'wb')
            if f is None:
                print('Cannot create %s' % out_name)
                return 1
            gdal.VSIFWriteL(vrt_xml, 1, len(vrt_xml), f)
            gdal.VSIFCloseL(f)
            if verbose_vrt:
                vrt_ds = gdal.Open(out_name, gdal.GA_Update)
                vrt_ds.SetMetadata(vrt_ds.GetMetadata())
            else:
                vrt_ds = gdal.Open(out_name)
            if vrt_ds is None:
                return 1

            return 0

        vrt_ds = gdal.Open(vrt_xml)
        out_ds = gdal.GetDriverByName(frmt).CreateCopy(out_name, vrt_ds, 0, creation_options, callback=callback)
        if out_ds is None:
            return 1
        return 0



# # 6s大气校正
# def AtmosphericCorrection(BandId,metedata,config,SatelliteID,SensorID):
#     # global metedata,config,SatelliteID,SensorID
#     #读取头文件
#     dom = xml.dom.minidom.parse(metedata)
#
#     # 6S模型
#     s = SixS()
#
#     # 传感器类型 自定义
#     s.geometry = Geometry.User()
#     s.geometry.solar_z = 90-float(dom.getElementsByTagName('SolarZenith')[0].firstChild.data)
#     s.geometry.solar_a = float(dom.getElementsByTagName('SolarAzimuth')[0].firstChild.data)
#     # s.geometry.view_z = float(dom.getElementsByTagName('SatelliteZenith')[0].firstChild.data)
#     # s.geometry.view_a = float(dom.getElementsByTagName('SatelliteAzimuth')[0].firstChild.data)
#     s.geometry.view_z = 0
#     s.geometry.view_a = 0
#     # 日期
#     DateTimeparm = dom.getElementsByTagName('CenterTime')[0].firstChild.data
#     DateTime = DateTimeparm.split(' ')
#     Date = DateTime[0].split('-')
#     s.geometry.month = int(Date[1])
#     s.geometry.day = int(Date[2])
#
#     # print(s.geometry)
#     # 中心经纬度
#     TopLeftLat = float(dom.getElementsByTagName('TopLeftLatitude')[0].firstChild.data)
#     TopLeftLon = float(dom.getElementsByTagName('TopLeftLongitude')[0].firstChild.data)
#     TopRightLat = float(dom.getElementsByTagName('TopRightLatitude')[0].firstChild.data)
#     TopRightLon = float(dom.getElementsByTagName('TopRightLongitude')[0].firstChild.data)
#     BottomRightLat = float(dom.getElementsByTagName('BottomRightLatitude')[0].firstChild.data)
#     BottomRightLon = float(dom.getElementsByTagName('BottomRightLongitude')[0].firstChild.data)
#     BottomLeftLat = float(dom.getElementsByTagName('BottomLeftLatitude')[0].firstChild.data)
#     BottomLeftLon = float(dom.getElementsByTagName('BottomLeftLongitude')[0].firstChild.data)
#
#     ImageCenterLat = (TopLeftLat + TopRightLat + BottomRightLat + BottomLeftLat) / 4
#
#     # 大气模式类型
#     if ImageCenterLat > -15 and ImageCenterLat < 15:
#         s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.Tropical)
#
#     if ImageCenterLat > 15 and ImageCenterLat < 45:
#         if s.geometry.month > 4 and s.geometry.month < 9:
#             s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.MidlatitudeSummer)
#         else:
#             s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.MidlatitudeWinter)
#
#     if ImageCenterLat > 45 and ImageCenterLat < 60:
#         if s.geometry.month > 4 and s.geometry.month < 9:
#             s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.SubarcticSummer)
#         else:
#             s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.SubarcticWinter)
#
#     # 气溶胶类型大陆
#     s.aero_profile = AtmosProfile.PredefinedType(AeroProfile.Continental)
#
#     # 下垫面类型
#     s.ground_reflectance = GroundReflectance.HomogeneousLambertian(0.36)
#
#     # 550nm气溶胶光学厚度,对应能见度为40km
#     s.aot550 = 0.14497
#
#     # 通过研究去区的范围去求DEM高度。
#     pointUL = dict()
#     pointDR = dict()
#     pointUL["lat"] = max(TopLeftLat,TopRightLat,BottomRightLat,BottomLeftLat)
#     pointUL["lon"] = min(TopLeftLon,TopRightLon,BottomRightLon,BottomLeftLon)
#     pointDR["lat"] = min(TopLeftLat,TopRightLat,BottomRightLat,BottomLeftLat)
#     pointDR["lon"] = max(TopLeftLon,TopRightLon,BottomRightLon,BottomLeftLon)
#     meanDEM = (MeanDEM(pointUL, pointDR)) * 0.001
#
#     # 研究区海拔、卫星传感器轨道高度
#     s.altitudes = Altitudes()
#     s.altitudes.set_target_custom_altitude(meanDEM)
#     s.altitudes.set_sensor_satellite_level()
#
#     # 校正波段（根据波段名称）
#     if BandId == 1:
#         SRFband = config["Parameter"][SatelliteID][SensorID]["SRF"]["1"]
#         s.wavelength = Wavelength(0.450,0.520,SRFband)
#
#     elif BandId == 2:
#         SRFband = config["Parameter"][SatelliteID][SensorID]["SRF"]["2"]
#
#         s.wavelength = Wavelength(0.520,0.590,SRFband)
#
#     elif BandId == 3:
#         SRFband = config["Parameter"][SatelliteID][SensorID]["SRF"]["3"]
#
#         s.wavelength = Wavelength(0.630,0.690,SRFband)
#
#     elif BandId == 4:
#         SRFband = config["Parameter"][SatelliteID][SensorID]["SRF"]["4"]
#         s.wavelength = Wavelength(0.770,0.890,SRFband)
#
#     s.atmos_corr = AtmosCorr.AtmosCorrLambertianFromReflectance(-0.1)
#
#     # 运行6s大气模型
#     s.run()
#     xa = s.outputs.coef_xa
#     xb = s.outputs.coef_xb
#     xc = s.outputs.coef_xc
#     # x = s.outputs.values
#     return (xa, xb, xc)
#
#

