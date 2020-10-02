# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 20:04:45 2017

@author: Rohil

#originally for substrate landscape, now works for any set of genes you'd like 
"""
import pandas as pd
import os
from matplotlib import style
from ubqClasses.RawMutationDataFrameRepository import RawMutationDataFrameRepository
from ubqFunctions.grabPortalData import grabPortalData

style.use('ggplot')

## the shenanigans at the end of the read_csv line are to convert from dataframe to series
#geneList = pd.read_csv('/Users/Rohil/Documents/Young Dawgs/substrate landscape genes/ubqSubstrateCanIsoformGenes.txt', index_col = False, header = None).transpose().iloc[0]

# define your gene list
#geneList = pd.read_excel('/Users/Rohil/Documents/Young Dawgs/ubq gene lists/machinery_family_genes/final_machinery_gene_list.xlsx').gene_symbol
geneList = pd.read_excel('/Users/Rohil/Documents/Young Dawgs/ubq gene lists/myc_genes.xlsx').gene_symbol

ontologyFilePath = '/Users/Rohil/Documents/Young Dawgs/GSEA Ubiquitin Genes'
ontologyFileList = os.listdir(ontologyFilePath)

tcgaStudyDF = pd.read_csv('/Users/Rohil/Documents/Young Dawgs/cBioPortal Cancer Studies/tcgaStudyDF.csv')

colorList = ['yellow', '#A3E4D7', '#5DADE2', '#16A085', '#0B5345',  'blue', '#D2B4DE', '#FCF3CF', '#F7DC6F', '#E67E22', '#922B21', 'red', '#9B59B6', '#633974', '#A04000', 'magenta']


#cutoffNum = 10
#weightCutoffNum = .1
#geneCutoffNum = 5


def getStudyData(studyDF, run_study):
    
    for x in studyDF['study_id']:
        cancer_study = str(x)
        
        if cancer_study == 'coadread_tcga_pub': #because the tag for coad is different - due to cBioPortal error
            pointUrl = 'http://www.cbioportal.org/webservice.do?cmd=getMutationData&genetic_profile_id=%s_mutations&case_set_id=%s_cna_seq&gene_list=' % (cancer_study, cancer_study) #note that only samples with cna and sequencing data are used as to
            cnaUrl = 'http://www.cbioportal.org/webservice.do?cmd=getProfileData&genetic_profile_id=%s_gistic&case_set_id=%s_cna_seq&gene_list=' % (cancer_study, cancer_study)       #prevent sample bias
            grabPortalData(pointUrl, cnaUrl, geneList, rawMutationDataFrameRepository, cancer_study, run_study)
        
        else:
        
            pointUrl = 'http://www.cbioportal.org/webservice.do?cmd=getMutationData&genetic_profile_id=%s_mutations&case_set_id=%s_cnaseq&gene_list=' % (cancer_study, cancer_study) #note that only samples with cna and sequencing data are used as to
            cnaUrl = 'http://www.cbioportal.org/webservice.do?cmd=getProfileData&genetic_profile_id=%s_gistic&case_set_id=%s_cnaseq&gene_list=' % (cancer_study, cancer_study)       #prevent sample bias
            grabPortalData(pointUrl, cnaUrl, geneList, rawMutationDataFrameRepository, cancer_study, run_study)
    
    return rawMutationDataFrameRepository
        
     
def getData(studyDF, study):
        
    try:
        
        rawMutationDataFrameRepository.rawMutationDF = pd.read_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/rawDF.txt')
        rawMutationDataFrameRepository.ontologyDF = pd.read_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/mappedDF.txt')
    
    except:
        
        getStudyData(studyDF, study)   
        rawMutationDataFrameRepository.assimilateRawMutationDF(study)
        rawMutationDataFrameRepository.groupByCancerType(study)
        #rawMutationDataF rameRepository.groupByOntology(study)
        
    finally:
        '''
        rawMutationDataFrameRepository.getMutData(category = 'gene_symbol', study = study, plot = True, cutoff = 0)
        #rawMutationDataFrameRepository.getOntologyMutationData(study, plot = True)
        #rawMutationDataFrameRepository.getOntologyGeneData(study, plot = True)
        #rawMutationDataFrameRepository.getOntologyMutWeightData(study, plot = True)
        rawMutationDataFrameRepository.getMutData(category = 'cancer_type', study = study, plot = True, cutoff = 0)
        rawMutationDataFrameRepository.getGeneData('cancer_type', study, plot = True, cutoff = 0)
        #rawMutationDataFrameRepository.getCancerOntologyData(study, plot = True)
        rawMutationDataFrameRepository.getMutWeightData(category = 'cancer_type', study = study, plot = True, cutoff = 0)
        rawMutationDataFrameRepository.getMutData(category = 'cancer_supertype', study = study, plot = True, cutoff = 0)
        rawMutationDataFrameRepository.getGeneData('cancer_supertype', study, plot = True, cutoff = 0)
        #rawMutationDataFrameRepository.getCancerSupertypeOntologyData(study, plot = True)
        rawMutationDataFrameRepository.getMutWeightData(category = 'cancer_supertype', study = study, plot = True, cutoff = 0)
        rawMutationDataFrameRepository.getMutData(category = 'cancer_site', study = study, plot = True, cutoff = 0)
        rawMutationDataFrameRepository.getGeneData('cancer_site', study, plot = True, cutoff = 0)
        #rawMutationDataFrameRepository.getCancerSiteOntologyData(study, plot = True)
        rawMutationDataFrameRepository.getMutWeightData(category = 'cancer_site', study = study, plot = True, cutoff = 0)
        '''
        pass

    
rawMutationDataFrameRepository = RawMutationDataFrameRepository()

#rawMutationDataFrameRepository.rawMutationDF = pd.read_csv('/Users/Rohil/Documents/Young Dawgs/substrate_landscape_tcga_ubq_site_exact_only/rawDF.txt')

#rawMutationDataFrameRepository.groupByOntology('substrate_landscape_tcga_ubq_site_exact_only')

#getData(tcgaStudyDF, 'substrate_landscape_tcga_ubq_site_exact_only')
getData(tcgaStudyDF, 'myc_landscape')
