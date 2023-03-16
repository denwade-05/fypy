# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : config.py

@Modify Time :  2022/10/27 16:37   

@Author      : Lee    

@Version     : 1.0   

@Description :

'''
import os
import sys
import numpy as np
import datetime

from fypy import parm

exedir = os.path.abspath(list(parm.__path__)[0])

FontTTF = os.path.join(exedir, 'font', 'simsun.ttf')
if not os.path.isfile(FontTTF) :
    raise Exception('字体文件不存在【%s】' %(FontTTF))
