# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 21:44:44 2017

@author: Rohil
"""
import pandas as pd
from urllib.request import urlopen
from io import StringIO
from GSEA_ubiquitin_database.getGeneOntologies import ontologyDict
import os
from ubqFunctions.ensureDirectory import ensureDirectory


ontologyFilePath = '/Users/Rohil/Documents/Young Dawgs/GSEA Ubiquitin Genes'
ontologyFileList = os.listdir(ontologyFilePath)

def grabPortalData(pointCancerUrl, cnaCancerUrl, uniqueGenes, instance, cancer_study, run_study):
    
    
    for i in range(0, uniqueGenes.size, 100):
        
        element = uniqueGenes.iloc[i:i+99].tolist()
        element = ','.join(element)
        
        
        try:
        #open the url, read the output, and convert from binary to utf-8
        
            url = "".join((pointCancerUrl + str(element)).split())
            pointPortalData = urlopen(url).read().decode("utf-8")
            #find the first line where 'entrez_gene_id' appears... this is the header column
            #also convert it to a string io object so it can be fed into read_csv
            #pointPortalData = StringIO(pointPortalData[pointPortalData.find('entrez_gene_id'):])
            
            #read the data into a dataframe
            pointDF = pd.read_csv(StringIO(pointPortalData[pointPortalData.find('entrez_gene_id'):]), sep = '\t', error_bad_lines = False)
    
            pointDF.drop(pointDF.columns[[0,3,4,6,9,10,11,12,13,14,15,16,17,18,19,20,21]], axis = 1, inplace = True)
    
            pointDF.rename(columns = {'case_id' : 'sample_id'}, inplace = True)
            
            pointDF['mutation_type'] = pointDF['mutation_type'].str.lower()
            
            pointDF['study_id'] = cancer_study  
    
            pointDF['pointORcna'] = 'point'
    
            pointDF['ontology'] = pointDF['gene_symbol'].map(ontologyDict)
    
            instance.appendPointDF(pointDF)
            
        except:
            ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + run_study + '/grabPortalDataErrors.txt')

            errorFile = open('/Users/Rohil/Documents/Young Dawgs/' + run_study + '/grabPortalDataErrors.txt', 'a')
            errorFile.write(cancer_study + ' returned no point mutation data for ' + element)
        
        
        try:
        
            url = "".join((cnaCancerUrl + str(element)).split())

            cnaPortalData = urlopen(url).read().decode("utf-8")
        
            cnaPortalData = StringIO(cnaPortalData[cnaPortalData.find('GENE_ID'):])
        
            cnaDF = pd.read_csv(cnaPortalData, sep = '\t', error_bad_lines = False)
        
            cnaDF.drop('GENE_ID', axis = 1, inplace = True)

            cnaDF.rename(columns = {'COMMON' : 'gene_symbol'}, inplace = True)

            cnaDF = cnaDF.melt('gene_symbol', var_name = 'sample_id', value_name = 'mutation_type')
            cnaDF['mutation_type'] = cnaDF['mutation_type'].map({-2.0 : 'homozygous_del', -1.0 : 'hemizygous_del', 1.0 : 'gain', 2.0 : 'high_lvl_amplification' })
            cnaDF.dropna(axis = 0, inplace = True)
        
            cnaDF['pointORcna'] = 'cna'
            cnaDF['study_id'] = cancer_study  
        
            cnaDF['ontology'] = cnaDF['gene_symbol'].map(ontologyDict)
        
            instance.appendCnaDF(cnaDF)
        
        except:
            ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + run_study + '/grabPortalDataErrors.txt')

            errorFile = open('/Users/Rohil/Documents/' + run_study + '/grabPortalDataErrors.txt', 'a')
            errorFile.write(cancer_study + ' returned no cna mutation data for ' + element)
