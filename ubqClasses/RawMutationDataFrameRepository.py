# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 21:49:50 2017

@author: Rohil
"""
import pandas as pd
import matplotlib.pyplot as plt
import ubqMutationData.ubqFunctions.ubqPlots as uplot
from ubqMutationData.ubqFunctions.mutationFiltering import combineRepeatMutations, getMutationSupertype
from ubqMutationData.ubqFunctions.ensureDirectory import ensureDirectory
from GSEA_ubiquitin_database.getGeneOntologies import ontologyDict
import numpy as np


class RawMutationDataFrameRepository(object):
    
    pointDF = None
    cnaDF = None
    rawMutationDF = None
    ontologyDF = None
    
    
    def __init__(self):
        
        self.pointDF = pd.DataFrame()
        self.cnaDF = pd.DataFrame()
        
    
    def appendPointDF(self, pointDF):
        #if this is the first time pointDF is being created in grabPortalData, intialize the self.pointDF as the pointDF from grabPortalData
        if self.pointDF.empty:
            self.pointDF = pointDF
        
        else:
            self.pointDF = self.pointDF.append(pointDF, ignore_index = True)
        
    
    def appendCnaDF(self, cnaDF):
        
        if self.cnaDF.empty:
            self.cnaDF = cnaDF
        else:
            self.cnaDF = self.cnaDF.append(cnaDF, ignore_index = True)
            
            
    def assimilateRawMutationDF(self, study):
        
        self.rawMutationDF = self.pointDF.append(self.cnaDF, ignore_index = True)
        self.rawMutationDF = self.rawMutationDF[['gene_symbol', 'pointORcna', 'mutation_type', 'amino_acid_change', 'functional_impact_score', 'sample_id', 'study_id', 'ontology']]
        self.rawMutationDF.sort_values('gene_symbol', inplace = True)
        
        self.rawMutationDF = combineRepeatMutations(self.rawMutationDF)

        ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/')
        #self.rawMutationDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/rawDF.txt', index = False)
    
    
    def groupByCancerType(self, study):

        studyToCancerTypeDF = pd.read_csv('/Users/Rohil/Documents/Young Dawgs/cBioPortal Cancer Studies/tcgaStudyDF.csv') # change file based on which studies are being used
        self.rawMutationDF = pd.merge(self.rawMutationDF, studyToCancerTypeDF, how = 'left')
        self.rawMutationDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/rawDF.txt', index = False)
    
    
    def groupByOntology(self, study):
           
        self.ontologyDF = self.rawMutationDF
        
        #drops and resets index
        self.ontologyDF = self.ontologyDF.reset_index(drop = True)
        
        #maps a list of ontologies to each row based on the gene_symbol
        self.ontologyDF['ontology'] = self.ontologyDF['gene_symbol'].map(ontologyDict)
        
        #if any rows do not have a value for 'ontology' set that value to a list containing 'no_ontology_found'
        self.ontologyDF.loc[self.ontologyDF['ontology'].isnull(),['ontology']] = self.ontologyDF.loc[self.ontologyDF['ontology'].isnull(), 'ontology'].apply(lambda x: ['no_ontology_found'])
        
        #get number of ontologies for each row
        dfLength = self.ontologyDF['ontology'].str.len()

        #convert ontology list to numpy array
        dfValues = self.ontologyDF['ontology'].values
        
        #get an array with the indices repeated by the number of ontologies for each row
        dfIndex = np.repeat(self.ontologyDF.index, dfLength)
        
        #expand rows by duplicated index values
        self.ontologyDF = self.ontologyDF.loc[dfIndex]
        
        #add ontology values, flattened from the list
        self.ontologyDF['ontology'] = np.concatenate(dfValues)
        
        #reset the index again, since it currently contains repeat values
        self.ontologyDF = self.ontologyDF.reset_index(drop = True)
        
        ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/')
        
        self.ontologyDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/mappedDF.txt', index = False)        
    
               
    def getMutData(self, category, study, plot, cutoff):
        
        if category == 'ontology':
            
            df = pd.crosstab(self.ontologyDF[category], self.ontologyDF['mutation_type'])
        
        else:
             df = pd.crosstab(self.rawMutationDF[category], self.rawMutationDF['mutation_type'])
             
        df['sum'] = df.sum(axis = 1)
        
        df = getMutationSupertype(df)
        
        df.sort_values(by = ['sum'], ascending = False, inplace = True)
                
        ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/' + category + '/')

        df.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/' + category + '/' + category + '_mutation_data.txt')    
        
        if plot == True:
        
            #take only the top 50 rows... only applies to gene mutation data but doesn't effect other categories
            df = df.head(50).drop('sum', axis = 1)
                        
            uplot.mutPlot(df, study, category, cutoff = 0, weighted = False)
            uplot.mutSupertypeSumPlot(df, study, category)
            uplot.mutSupertypeIndividualPlot(df, study, category)
            
    
    def getGeneData(self, category, study, plot, cutoff):
        
        if category == 'ontology':
            subcategoryList = self.ontologyDF[category].unique()
            
        else:
            subcategoryList = self.rawMutationDF[category].unique()
            
        for subcat in subcategoryList:
            
            if category == 'ontology':
                df = self.ontologyDF.loc[self.ontologyDF[category] == subcat]
                
            else:
                df = self.rawMutationDF.loc[self.rawMutationDF[category] == subcat]
                
            geneDF = pd.crosstab(df['gene_symbol'], df['mutation_type'])
            geneDF['sum'] = geneDF.sum(axis = 1)
            geneDF.sort_values(by = 'sum', ascending = False, inplace = True)
                        
            geneWeight = df.groupby('gene_symbol').mutation_type.value_counts(normalize = True)
            geneWeightDF = geneWeight.unstack().fillna(0.0)
            geneWeightDF = geneWeightDF.reindex(geneDF.index)
            
            ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/' + category + '/')
           
            geneDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/' + category + '/' + category + '_gene_data.txt')
            geneWeightDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/' + category + '/' + category + '_gene_weight_data.txt')
            
            if plot == True:
                
                geneDF.drop('sum', axis = 1, inplace = True)
                
                uplot.genePlot(geneDF, study, category, subcat, cutoff, weighted = False)
                uplot.genePlot(geneWeightDF, study, category, subcat, cutoff, weighted = True)
                
    
    def getMutWeightData(self, category, study, plot, cutoff):
        
        #get a groupby object containing all the mutations grouped by the category
        #then gets the avg mutations per sample or gene
        
        if category == 'ontology':
            
            grouped = self.ontologyDF.groupby(category)
            mutWeight = grouped.mutation_type.value_counts() / grouped.gene_symbol.nunique()
            
        else:
            
            grouped = self.rawMutationDF.groupby(category)
            mutWeight = grouped.mutation_type.value_counts() / grouped.sample_id.nunique()
            
        
        #unstacks groupby object and converts to dataframe to plot
        #fills NaN values with 0... NaN values here mean that the category did not have that type of mutation
        mutWeightDF = mutWeight.unstack().fillna(0.0)
        
        mutWeightDF['sum'] = mutWeightDF.sum(axis = 1)
            
        mutWeightDF.sort_values(by = 'sum', ascending = False, inplace = True)
        
        ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/' + category + '/')
        
        mutWeightDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/' + category + '/' + category + '_mut_weight_data.txt')
        
        if plot == True:
            
            mutWeightDF.drop('sum', axis = 1, inplace = True)
            
            uplot.mutPlot(mutWeightDF, study, category, cutoff, weighted = True)

    
    def getOntologyData(self, category, study, plot): #DEPRECATED... NO LONGER IN USE (TENTATIVE)
        
        for subcat in self.rawMutationDF[category].unique():
            
            df = self.ontologyDF.loc[self.ontologyDF[category] == subcat]
            
            ontologyDF = pd.crosstab(df['ontology'], df['mutation_type'])
            ontologyDF['sum'] = ontologyDF.sum(axis = 1)
            ontologyDF.sort_values(by = 'sum', ascending = False, inplace = True)
            
            ontologyGrouped = df.groupby('ontology')
            ontologyWeight = ontologyGrouped.mutation_type.value_counts() / ontologyGrouped.gene_symbol.nunique()
            ontologyWeightDF = ontologyWeight.unstack().fillna(0.0)
            ontologyWeightDF['sum'] = ontologyWeightDF.sum(axis = 1) 
            ontologyWeightDF.sort_values(by = 'sum', ascending = False, inplace = True)
            
            ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/'+ category + '/')
           
            ontologyDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/' + category + '/' + category + '_ontology_data.txt')
            ontologyWeightDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/' + category + '/' + category + '_ontology_weight_data.txt')
            
            if plot == True:
            
                ontologyDF.drop('sum', axis = 1, inplace = True)
                ontologyWeightDF.drop('sum', axis = 1, inplace = True)
                
                
                uplot.ontologyPlot(ontologyDF, study, category, subcat)
                uplot.ontologyWeightPlot(ontologyWeightDF, study, category, subcat)
            
