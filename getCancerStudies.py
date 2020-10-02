# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 10:19:22 2017

@author: Rohil
"""

import pandas as pd

cancerStudyDF = pd.read_csv('/Users/Rohil/Documents/Young Dawgs/cBioPortal Cancer Studies/cancerStudyDF.txt', sep = '\t')
tcgaStudiesDEPRECATED = pd.read_csv('/Users/Rohil/Documents/Young Dawgs/cBioPortal Cancer Studies/cancerStudyDF_TCGA.txt')
broadMskccStudiesDF = pd.read_csv('/Users/Rohil/Documents/Young Dawgs/cBioPortal Cancer Studies/cancerStudyDF_Broad_MSKCC.txt', sep = '\t')
otherInstituteStudies = pd.read_csv('/Users/Rohil/Documents/Young Dawgs/cBioPortal Cancer Studies/cancerStudyDF_Others.txt', sep = '\t')
nonEpithelialStudiesDF = pd.read_csv('/Users/Rohil/Documents/Young Dawgs/cBioPortal Cancer Studies/cancerStudyDF_nonEpithelial.txt', sep = '\t')

    
    