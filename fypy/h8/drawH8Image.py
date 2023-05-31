# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : drawH8Image.py

@Modify Time :  2022/11/10 14:45

@Author      : Lee    

@Version     : 1.0   

@Description :

'''
import os
import sys
import numpy as np
import datetime
from PIL import Image


def drawH8TrueColor(outname, vis006, vis005, vis004):


    r = vis006 * 255
    g = vis005 * 255
    b = vis004 * 255

    r[r<0] = 0
    g[g<0] = 0
    b[b<0] = 0

    r[r>255] = 255
    g[g>255] = 255
    b[b>255] = 255

    rgbArray = np.zeros((r.shape[0], r.shape[1], 4), dtype=np.float64)
    rgbArray[..., 0] = r
    rgbArray[..., 1] = g
    rgbArray[..., 2] = b
    rgbArray[..., 3] = 255

    img = Image.fromarray(rgbArray.astype(np.uint8))
    img.save(outname, quality=90)


    # im = Image.merge('RGB', (Image.fromarray(np.array(r, dtype=np.uint8), "L"),
    #                          Image.fromarray(np.array(g, dtype=np.uint8), "L"),
    #                          Image.fromarray(np.array(b, dtype=np.uint8), "L")))
    # im.save(outname, quality=95)
