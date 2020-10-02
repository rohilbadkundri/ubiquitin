# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 13:47:09 2017

@author: Rohil
"""
import os

def ensureDirectory(filePath):
    directory = os.path.dirname(filePath)
    if not os.path.exists(directory):
        os.makedirs(directory)