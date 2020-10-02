# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 22:51:43 2017

@author: Rohil
"""

class RawMutationDataFrameRepository(object):
    
    pointDF = None
    cnaDF = None
    rawMutationDF = None
    ontologyDF = None
    
    
    def __init__(self):
        
        self.pointDF = pd.DataFrame()
        self.cnaDF = pd.DataFrame()
        
    
    def appendPointDF(self, pointDF):
        
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
        
        self.rawMutationDF = self.pointDF.append(self.cnaDF)
        self.rawMutationDF = self.rawMutationDF[['gene_symbol', 'pointORcna', 'mutation_type', 'amino_acid_change', 'functional_impact_score', 'sample_id', 'study_id', 'ontology']]
        self.rawMutationDF.sort_values('gene_symbol')
        
        self.rawMutationDF = combineRepeatMutations(self.rawMutationDF)

        ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/')
        #self.rawMutationDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/rawDF.txt', index = False)
    
    
    def groupByCancerType(self, study):

        studyToCancerTypeDF = pd.read_csv('/Users/Rohil/Documents/Young Dawgs/cBioPortal Cancer Studies/tcgaStudyDF.txt', sep = '\t') # change file based on which studies are being used
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
    
    
    def getGeneData(self, study, plot):

        rawMutDF = pd.crosstab(self.rawMutationDF['gene_symbol'], self.rawMutationDF['mutation_type'])
        rawMutDF['sum'] = rawMutDF.sum(axis = 1)
        
        rawMutDF = getMutationSupertype(rawMutDF)
        
        rawMutDF.sort_values(by = ['sum'], ascending = False, inplace = True)
                
        ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/gene_symbol/')

        rawMutDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/gene_symbol/geneMutationData.txt')    

        if plot == True:
            
            rawMutDF = rawMutDF.head(50).drop('sum', axis = 1)
            
            mutPlot(rawMutDF, study, 'gene_symbol')
            
            mutSupertypeSumPlot(rawMutDF, study, 'gene_symbol')
            
            mutSupertypeIndividualPlot(rawMutDF, study, 'gene_symbol')

               
    def getOntologyMutationData(self, study, plot):
            
        ontologyMutDF = pd.crosstab(self.ontologyDF['ontology'], self.ontologyDF['mutation_type'])
        
        ontologyMutDF['sum'] = ontologyMutDF.sum(axis = 1)
        
        ontologyMutDF = getMutationSupertype(ontologyMutDF)
        
        ontologyMutDF.sort_values(by = ['sum'], ascending = False, inplace = True)
                
        ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/ontology/')

        ontologyMutDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/ontology/ontologyMutationData.txt')    
        
        if plot == True:
        
            ontologyMutDF.drop('sum', axis = 1, inplace = True)
                        
            mutPlot(ontologyMutDF, study, 'ontology')
            
            mutSupertypeSumPlot(ontologyMutDF, study, 'ontology')
            
            mutSupertypeIndividualPlot(ontologyMutDF, study, 'ontology')
    
    
    def getOntologyGeneData(self, study, plot):
       
        for ontology in self.ontologyDF.ontology.unique():
            
            df = self.ontologyDF.loc[self.ontologyDF['ontology'] == ontology]
            
            ontologyGeneDF = pd.crosstab(df['gene_symbol'], df['mutation_type'])
            ontologyGeneDF['sum'] = ontologyGeneDF.sum(axis = 1)
            ontologyGeneDF.sort_values(by = 'sum', ascending = False, inplace = True)
            
            ontologyGeneWeight = df.groupby('gene_symbol').mutation_type.value_counts(normalize = True) 
            ontologyGeneWeightDF = ontologyGeneWeight.unstack().fillna(0.0)
            ontologyGeneWeightDF.reindex(ontologyGeneDF.index)
            
            ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/ontology/')
            
            ontologyGeneDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/ontology/ontologyGeneData.txt')
            ontologyGeneWeightDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/ontology/ontologyGeneWeightData.txt')
            
            if plot == True:
                
                ontologyGeneDF.drop('sum', axis = 1, inplace = True)
                
                genePlot(ontologyGeneDF, study, 'ontology', ontology)
                
                geneWeightPlot(ontologyGeneWeightDF, study, 'ontology', ontology)
        
# =============================================================================
#         ontologyGeneDF = pd.crosstab(self.ontologyDF['ontology'], self.ontologyDF['gene_symbol']).replace(0, np.nan)
#         ontologyGeneDF['sum'] = ontologyGeneDF.sum(axis = 1)
#         
#         #calculates avg mutations per gene by taking value of sum column for that row and dividing by num of non-NaN occurences on
#         #that line; we subtract one because the sum column adds an extra non-gene occurence
#         ontologyGeneDF['avg_mut_per_gene'] = ontologyGeneDF['sum'] / (ontologyGeneDF.count(axis = 1) - 1)
#     
#         ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/ontology/')
#         
#         ontologyGeneDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/ontology/ontologyGeneData.txt')    
# 
#         if plot == True:
#         
#             ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/ontology/ontology_gene_images/')
#         
#             ontologyGeneDF.sort_values(by = ['avg_mut_per_gene'], ascending = False, inplace = True)
#             
#             #mutWeightPlot(ontologyGeneDF, study, category = 'ontology', yy = 'avg_mut_per_gene')
#             
#             ontologyGeneDF.drop(['sum', 'avg_mut_per_gene'], axis = 1, inplace = True)
#         
#             for i in range(0, ontologyGeneDF.shape[0]):
#             
#                 genePlot(ontologyGeneDF, study, 'ontology', i)
# =============================================================================
                
    
    def getOntologyMutWeightData(self, study, plot):
        
        #get a groupby object containing all the mutations grouped by ontology
        ontologyGrouped = self.ontologyDF.groupby('ontology')
        
        #get avg num of each type of mutation per gene
        ontologyMutWeight = ontologyGrouped.mutation_type.value_counts() / ontologyGrouped.gene_symbol.nunique()
        
        #unstacks groupby object and converts to dataframe to plot
        #fills NaN values with 0... NaN values here mean that ontology did not have that type of mutation
        ontologyMutWeightDF = ontologyMutWeight.unstack().fillna(0.0)
        
        ontologyMutWeightDF['sum'] = ontologyMutWeightDF.sum(axis = 1)
            
        ontologyMutWeightDF.sort_values(by = 'sum', ascending = False, inplace = True)
        
        ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/ontology/')
        
        ontologyMutWeight.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/ontology/ontologyMutWeightData.txt')
        
        if plot == True:
            
            ontologyMutWeightDF.drop('sum', axis = 1, inplace = True)
            
            mutWeightPlot(ontologyMutWeightDF, study, 'ontology', 'avg_mut_per_gene')
        
        
    
    def getCancerMutData(self, study, plot):
         
         cancerMutDF = pd.crosstab(self.rawMutationDF['cancer_type'], self.rawMutationDF['mutation_type'])
         cancerMutDF['sum'] = cancerMutDF.sum(axis = 1)
         
         cancerMutDF = getMutationSupertype(cancerMutDF)
         
         cancerMutDF.sort_values(by = ['sum'], ascending = False, inplace = True)
                  
         ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_type/')
                 
         cancerMutDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_type/cancerMutData.txt')    
                  
         if plot == True:
            
            cancerMutDF.drop('sum', axis = 1, inplace = True)
                        
            mutPlot(cancerMutDF, study, 'cancer_type')
            mutSupertypeSumPlot(cancerMutDF, study, 'cancer_type')
            mutSupertypeIndividualPlot(cancerMutDF, study, 'cancer_type')            
 
         
    def getCancerGeneData(self, study, plot):
        
        for cancer_type in self.rawMutationDF.cancer_type.unique():
            
            df = self.rawMutationDF.loc[self.rawMutationDF['cancer_type'] == cancer_type]
            
            cancerGeneDF = pd.crosstab(df['gene_symbol'], df['mutation_type'])
            cancerGeneDF['sum'] = cancerGeneDF.sum(axis = 1)
            cancerGeneDF.sort_values(by = 'sum', ascending = False, inplace = True)
            
            cancerGeneWeight = df.groupby('gene_symbol').mutation_type.value_counts(normalize = True)
            cancerGeneWeightDF = cancerGeneWeight.unstack().fillna(0.0)
            cancerGeneWeightDF = cancerGeneWeightDF.reindex(cancerGeneDF.index)
            
            ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_type/')
           
            cancerGeneDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_type/cancer_typeGeneData.txt')
            cancerGeneWeightDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_type/cancer_typeGeneWeightData.txt')
            
            if plot == True:
                
                cancerGeneDF.drop('sum', axis = 1, inplace = True)
                
                genePlot(cancerGeneDF, study, 'cancer_type', cancer_type)  
                geneWeightPlot(cancerGeneWeightDF, study, 'cancer_type', cancer_type)
                
    
    def getCancerOntologyData(self, study, plot):
        
        for cancer_type in self.rawMutationDF.cancer_type.unique():
            
            df = self.ontologyDF.loc[self.ontologyDF['cancer_type'] == cancer_type]
            
            cancerOntologyDF = pd.crosstab(df['ontology'], df['mutation_type'])
            cancerOntologyDF['sum'] = cancerOntologyDF.sum(axis = 1)
            cancerOntologyDF.sort_values(by = 'sum', ascending = False, inplace = True)
            
            cancerOntologyGrouped = df.groupby('ontology')
            cancerOntologyWeight = cancerOntologyGrouped.mutation_type.value_counts() / cancerOntologyGrouped.gene_symbol.nunique()
            cancerOntologyWeightDF = cancerOntologyWeight.unstack().fillna(0.0)
            cancerOntologyWeightDF['sum'] = cancerOntologyWeightDF.sum(axis = 1) 
            cancerOntologyWeightDF.sort_values(by = 'sum', ascending = False, inplace = True)
            
            ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_type/')
           
            cancerOntologyDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_type/cancer_typeOntologyData.txt')
            cancerOntologyWeightDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_type/cancer_typeOntologyWeightData.txt')
            
            if plot == True:
            
                cancerOntologyDF.drop('sum', axis = 1, inplace = True)
                cancerOntologyWeightDF.drop('sum', axis = 1, inplace = True)
                
                
                ontologyPlot(cancerOntologyDF, study, 'cancer_type', cancer_type)
                ontologyWeightPlot(cancerOntologyWeightDF, study, 'cancer_type', cancer_type)
                
# =============================================================================
#         cancerGeneDF = pd.crosstab(self.rawMutationDF['cancer_type'], self.rawMutationDF['gene_symbol']).replace(0, np.nan)
#         cancerGeneDF['sum'] = cancerGeneDF.sum(axis = 1)
#         
#         #calculates avg mutations per gene by taking value of sum column for that row and dividing by num of non-NaN occurences on
#         #that line; we subtract one because the sum column adds an extra non-gene occurence
#         cancerGeneDF['avg mut per gene'] = cancerGeneDF['sum'] / (cancerGeneDF.count(axis = 1) - 1)
#     
#         ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_type/')
#         
#         cancerGeneDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_type/cancerGeneData.txt')    
# 
#         if plot == True:
#         
#             ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_type/cancer_type_gene_images/')
#         
#             cancerGeneDF.drop(['sum', 'avg mut per gene'], axis = 1, inplace = True)
#         
#             for i in range(0, cancerGeneDF.shape[0]):
#             
#                 genePlot(cancerGeneDF, study, 'cancer_type', i)
# =============================================================================

    
    
    def getCancerMutWeightData(self, study, plot):
        
        #get a groupby object containing all the mutations grouped by cancer_type
        cancerGrouped = self.rawMutationDF.groupby('cancer_type')
        
        #get avg num of each type of mutation per sample
        cancerMutWeight = cancerGrouped.mutation_type.value_counts() / cancerGrouped.sample_id.nunique()
        
        #unstacks groupby object and converts to dataframe to plot
        #fills NaN values with 0... NaN values here mean that ontology did not have that type of mutation
        cancerMutWeightDF = cancerMutWeight.unstack().fillna(0.0)
        
        cancerMutWeightDF['sum'] = cancerMutWeightDF.sum(axis = 1)
            
        cancerMutWeightDF.sort_values(by = 'sum', ascending = False, inplace = True)
        
        ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_type/')
        
        cancerMutWeightDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_type/cancer_typeMutWeightData.txt')
        
        if plot == True:
            
            cancerMutWeightDF.drop('sum', axis = 1, inplace = True)
            
            mutWeightPlot(cancerMutWeightDF, study, 'cancer_type', 'avg_mut_per_sample')

# =============================================================================
#         cancerMutDF = pd.crosstab(self.rawMutationDF['cancer_type'], self.rawMutationDF['mutation_type'])
#         cancerMutDF['sum'] = cancerMutDF.sum(axis = 1)
#          
#         cancerSampleDF = pd.crosstab(self.rawMutationDF['cancer_type'], self.rawMutationDF['sample_id'])
#         cancerSampleDF.replace(0, np.nan, inplace = True)
#         cancerSampleDF['unique_samples'] = cancerSampleDF.count(axis = 1)
#         
#         ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_type/')
# 
#         cancerSampleDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_type/cancerSampleData.txt')
#                 
#         cancerMutWeightDF = cancerMutDF.join(cancerSampleDF['unique_samples'])
#         
#         cancerMutWeightDF['avg_mut_per_sample'] = cancerMutWeightDF['sum'] / cancerMutWeightDF['unique_samples']
#                 
#         cancerMutWeightDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_type/cancerMutationData.txt')    
# 
#         if plot == True:
#         
#             mutWeightPlot(cancerMutWeightDF, study, 'cancer_type', 'avg_mut_per_sample')
# =============================================================================
  

    def getCancerSupertypeMutData(self, study, plot):
        
         supcancerMutDF = pd.crosstab(self.rawMutationDF['cancer_supertype'], self.rawMutationDF['mutation_type'])
         supcancerMutDF['sum'] = supcancerMutDF.sum(axis = 1)
         
         supcancerMutDF = getMutationSupertype(supcancerMutDF)
         
         supcancerMutDF.sort_values(by = ['sum'], ascending = False, inplace = True)
                  
         ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_supertype/')
                  
         if plot == True:
            
            supcancerMutDF.drop('sum', axis = 1, inplace = True)
                        
            mutPlot(supcancerMutDF, study, category = 'cancer_supertype')
            mutSupertypeSumPlot(supcancerMutDF, study, category = 'cancer_supertype')
            mutSupertypeIndividualPlot(supcancerMutDF, study, 'cancer_supertype')

            
    def getCancerSupertypeGeneData(self, study, plot):
        
        for cancer_supertype in self.rawMutationDF.cancer_supertype.unique():
            
            df = self.rawMutationDF.loc[self.rawMutationDF['cancer_supertype'] == cancer_supertype]
            
            supcancerGeneDF = pd.crosstab(df['gene_symbol'], df['mutation_type'])
            supcancerGeneDF['sum'] = supcancerGeneDF.sum(axis = 1)
            supcancerGeneDF.sort_values(by = 'sum', ascending = False, inplace = True)
                        
            supcancerGeneWeight = df.groupby('gene_symbol').mutation_type.value_counts(normalize = True)
            supcancerGeneWeightDF = supcancerGeneWeight.unstack().fillna(0.0)
            supcancerGeneWeightDF = supcancerGeneWeightDF.reindex(supcancerGeneDF.index)
            
            ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_supertype/')
           
            supcancerGeneDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_supertype/cancer_supertypeGeneData.txt')
            supcancerGeneWeightDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_supertype/cancer_supertypeGeneWeightData.txt')
            
            if plot == True:
                                
                supcancerGeneDF.drop('sum', axis = 1, inplace = True)
                
                genePlot(supcancerGeneDF, study, 'cancer_supertype', cancer_supertype)
                geneWeightPlot(supcancerGeneWeightDF, study, 'cancer_supertype', cancer_supertype)
                
    
    def getCancerSupertypeOntologyData(self, study, plot):
        
        for cancer_supertype in self.rawMutationDF.cancer_supertype.unique():
            
            df = self.ontologyDF.loc[self.ontologyDF['cancer_supertype'] == cancer_supertype]
            
            supcancerOntologyDF = pd.crosstab(df['ontology'], df['mutation_type'])
            supcancerOntologyDF['sum'] = supcancerOntologyDF.sum(axis = 1)
            supcancerOntologyDF.sort_values(by = 'sum', ascending = False, inplace = True)
            
            supcancerOntologyGrouped = df.groupby('ontology')
            supcancerOntologyWeight = supcancerOntologyGrouped.mutation_type.value_counts() / supcancerOntologyGrouped.gene_symbol.nunique()
            supcancerOntologyWeightDF = supcancerOntologyWeight.unstack().fillna(0.0)
            supcancerOntologyWeightDF['sum'] = supcancerOntologyWeightDF.sum(axis = 1) 
            supcancerOntologyWeightDF.sort_values(by = 'sum', ascending = False, inplace = True)
            
            ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_supertype/')
           
            supcancerOntologyDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_supertype/cancer_supertypeOntologyData.txt')
            supcancerOntologyWeightDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_supertype/cancer_supertypeOntologyWeightData.txt')

            if plot == True:
                                
                supcancerOntologyDF.drop('sum', axis = 1, inplace = True)
                supcancerOntologyWeightDF.drop('sum', axis = 1, inplace = True)
                
                ontologyPlot(supcancerOntologyDF, study, 'cancer_supertype', cancer_supertype)
                ontologyWeightPlot(supcancerOntologyWeightDF, study, 'cancer_supertype', cancer_supertype)
                
# =============================================================================
#         supcancerGeneDF = pd.crosstab(self.rawMutationDF['cancer_supertype'], self.rawMutationDF['gene_symbol']).replace(0, np.nan)
#         supcancerGeneDF['sum'] = supcancerGeneDF.sum(axis = 1)
#         
#         #calculates avg mutations per gene by taking value of sum column for that row and dividing by num of non-NaN occurences on
#         #that line; we subtract one because the sum column adds an extra non-gene occurence
#         supcancerGeneDF['avg mut per gene'] = supcancerGeneDF['sum'] / (supcancerGeneDF.count(axis = 1) - 1)
#     
#         ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_supertype/')
#         
#         supcancerGeneDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_supertype/cancerSupertypeGeneData.txt')    
# 
#         if plot == True:
#         
#             ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_supertype/cancer_supertype_gene_images/')
#         
#             supcancerGeneDF.drop(['sum', 'avg mut per gene'], axis = 1, inplace = True)
#         
#             for i in range(0, supcancerGeneDF.shape[0]):
#                 
#                 genePlot(supcancerGeneDF, study, 'cancer_supertype', i)
# =============================================================================
                

    def getCancerSupertypeMutWeightData(self, study, plot):
        
        #get a groupby object containing all the mutations grouped by cancer_supertype
        supcancerGrouped = self.rawMutationDF.groupby('cancer_supertype')
        
        #get avg num of each type of mutation per sample
        supcancerMutWeight = supcancerGrouped.mutation_type.value_counts() / supcancerGrouped.sample_id.nunique()
        
        #unstacks groupby object and converts to dataframe to plot
        #fills NaN values with 0... NaN values here mean that ontology did not have that type of mutation
        supcancerMutWeightDF = supcancerMutWeight.unstack().fillna(0.0)
        
        supcancerMutWeightDF['sum'] = supcancerMutWeightDF.sum(axis = 1)
        supcancerMutWeightDF.sort_values(by = 'sum', ascending = False, inplace = True)
        
        ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_supertype/')
        
        supcancerMutWeightDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_supertype/cancer_supertypeMutWeightData.txt')
        
        if plot == True:
            
            supcancerMutWeightDF.drop('sum', axis = 1, inplace = True)
            mutWeightPlot(supcancerMutWeightDF, study, 'cancer_supertype', 'avg_mut_per_sample')
        
# =============================================================================
#         supcancerMutDF = pd.crosstab(self.rawMutationDF['cancer_supertype'], self.rawMutationDF['mutation_type'])
#         supcancerMutDF['sum'] = supcancerMutDF.sum(axis = 1)
#          
#         supcancerSampleDF = pd.crosstab(self.rawMutationDF['cancer_supertype'], self.rawMutationDF['sample_id'])
#         supcancerSampleDF.replace(0, np.nan, inplace = True)
#         supcancerSampleDF['unique_samples'] = supcancerSampleDF.count(axis = 1)
#         
#         ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_supertype/')
# 
#         supcancerSampleDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_supertype/cancerSupertypeSampleData.txt')
#                 
#         supcancerMutWeightDF = supcancerMutDF.join(supcancerSampleDF['unique_samples'])
#         
#         supcancerMutWeightDF['avg_mut_per_sample'] = supcancerMutWeightDF['sum'] / supcancerMutWeightDF['unique_samples']
#                 
#         supcancerMutWeightDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_supertype/cancerSupertypeMutationData.txt')    
# 
#         if plot == True:
#         
#             mutWeightPlot(supcancerMutWeightDF, study, 'cancer_supertype', 'avg_mut_per_sample')
# =============================================================================
            
    
    def getCancerSiteMutData(self, study, plot):
        
         cancerSiteMutDF = pd.crosstab(self.rawMutationDF['cancer_site'], self.rawMutationDF['mutation_type'])
         cancerSiteMutDF['sum'] = cancerSiteMutDF.sum(axis = 1)
         
         cancerSiteMutDF = getMutationSupertype(cancerSiteMutDF)
         
         cancerSiteMutDF.sort_values(by = ['sum'], ascending = False, inplace = True)
                  
         ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_site/')
                  
         if plot == True:
            
            cancerSiteMutDF.drop('sum', axis = 1, inplace = True)
                        
            mutPlot(cancerSiteMutDF, study, 'cancer_site')
            mutSupertypeSumPlot(cancerSiteMutDF, study, 'cancer_site')
            mutSupertypeIndividualPlot(cancerSiteMutDF, study, 'cancer_site')
            
    
    def getCancerSiteGeneData(self, study, plot):

        for cancer_site in self.rawMutationDF.cancer_site.unique():
            
            df = self.rawMutationDF.loc[self.rawMutationDF['cancer_site'] == cancer_site]
            
            cancerSiteGeneDF = pd.crosstab(df['gene_symbol'], df['mutation_type'])
            cancerSiteGeneDF['sum'] = cancerSiteGeneDF.sum(axis = 1)
            cancerSiteGeneDF.sort_values(by = 'sum', ascending = False, inplace = True)
            
            cancerSiteGeneWeight = df.groupby('gene_symbol').mutation_type.value_counts(normalize = True)
            cancerSiteGeneWeightDF = cancerSiteGeneWeight.unstack().fillna(0.0)
            cancerSiteGeneWeightDF = cancerSiteGeneWeightDF.reindex(cancerSiteGeneDF.index)
            
            ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_site/')
           
            cancerSiteGeneDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_type/cancer_siteGeneData.txt')
            cancerSiteGeneWeightDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_site/cancer_siteGeneWeightData.txt')
            
            if plot == True:
                                
                cancerSiteGeneDF.drop('sum', axis = 1, inplace = True)
                
                genePlot(cancerSiteGeneDF, study, 'cancer_site', cancer_site)
                geneWeightPlot(cancerSiteGeneWeightDF, study, 'cancer_site', cancer_site)
                
    
    def getCancerSiteOntologyData(self, study, plot):
        
        for cancer_site in self.rawMutationDF.cancer_site.unique():
            
            df = self.ontologyDF.loc[self.ontologyDF['cancer_site'] == cancer_site]
            
            cancerSiteOntologyDF = pd.crosstab(df['ontology'], df['mutation_type'])
            cancerSiteOntologyDF['sum'] = cancerSiteOntologyDF.sum(axis = 1)
            cancerSiteOntologyDF.sort_values(by = 'sum', ascending = False, inplace = True)
            
            cancerSiteOntologyGrouped = df.groupby('ontology')
            cancerSiteOntologyWeight = cancerSiteOntologyGrouped.mutation_type.value_counts() / cancerSiteOntologyGrouped.gene_symbol.nunique()
            cancerSiteOntologyWeightDF = cancerSiteOntologyWeight.unstack().fillna(0.0)
            cancerSiteOntologyWeightDF['sum'] = cancerSiteOntologyWeightDF.sum(axis = 1) 
            cancerSiteOntologyWeightDF.sort_values(by = 'sum', ascending = False, inplace = True)
            
            ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_site/')
           
            cancerSiteOntologyDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_site/cancer_siteOntologyData.txt')
            cancerSiteOntologyWeightDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_site/cancer_siteOntologyWeightData.txt')
            
            if plot == True:
                                
                cancerSiteOntologyDF.drop('sum', axis = 1, inplace = True)
                cancerSiteOntologyWeightDF.drop('sum', axis = 1, inplace = True)
                
                ontologyPlot(cancerSiteOntologyDF, study, 'cancer_site', cancer_site)
                ontologyWeightPlot(cancerSiteOntologyWeightDF, study, 'cancer_site', cancer_site)
                
# =============================================================================
#         cancerSiteGeneDF = pd.crosstab(self.rawMutationDF['cancer_site'], self.rawMutationDF['gene_symbol']).replace(0, np.nan)
#         cancerSiteGeneDF['sum'] = cancerSiteGeneDF.sum(axis = 1)
#         
#         #calculates avg mutations per gene by taking value of sum column for that row and dividing by num of non-NaN occurences on
#         #that line; we subtract one because the sum column adds an extra non-gene occurence
#         cancerSiteGeneDF['avg mut per gene'] = cancerSiteGeneDF['sum'] / (cancerSiteGeneDF.count(axis = 1) - 1)
#     
#         ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/')
#         
#         cancerSiteGeneDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_site/cancerSiteGeneData.txt')    
# 
#         if plot == True:
#         
#             ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_site/cancer_site_gene_images/')
#         
#             cancerSiteGeneDF.drop(['sum', 'avg mut per gene'], axis = 1, inplace = True)
#         
#             for i in range(0, cancerSiteGeneDF.shape[0]):
#                 
#                 genePlot(cancerSiteGeneDF, study, 'cancer_site', i)
# =============================================================================
                
    
    def getCancerSiteMutWeightData(self, study, plot):
    
        #get a groupby object containing all the mutations grouped by cancer_site
        cancerSiteGrouped = self.rawMutationDF.groupby('cancer_site')
        
        #get avg num of each type of mutation per sample
        cancerSiteMutWeight = cancerSiteGrouped.mutation_type.value_counts() / cancerSiteGrouped.sample_id.nunique()
        
        #unstacks groupby object and converts to dataframe to plot
        #fills NaN values with 0... NaN values here mean that ontology did not have that type of mutation
        cancerSiteMutWeightDF = cancerSiteMutWeight.unstack().fillna(0.0)
        
        cancerSiteMutWeightDF['sum'] = cancerSiteMutWeightDF.sum(axis = 1)
        cancerSiteMutWeightDF.sort_values(by = 'sum', ascending = False, inplace = True)
        
        ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_site/')
        
        cancerSiteMutWeightDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_site/cancer_siteMutWeightData.txt')
        
        if plot == True:
            
            cancerSiteMutWeightDF.drop('sum', axis = 1, inplace = True)
            mutWeightPlot(cancerSiteMutWeightDF, study, 'cancer_site', 'avg_mut_per_sample')

# =============================================================================
#         cancerSiteMutDF = pd.crosstab(self.rawMutationDF['cancer_site'], self.rawMutationDF['mutation_type'])
#         cancerSiteMutDF['sum'] = cancerSiteMutDF.sum(axis = 1)
#          
#         cancerSiteSampleDF = pd.crosstab(self.rawMutationDF['cancer_site'], self.rawMutationDF['sample_id'])
#         cancerSiteSampleDF.replace(0, np.nan, inplace = True)
#         cancerSiteSampleDF['unique_samples'] = cancerSiteMutDF.count(axis = 1)
#         
#         ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_site/')
# 
#         cancerSiteSampleDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_site/cancerSiteSampleData.txt')
#                 
#         cancerSiteMutWeightDF = cancerSiteMutDF.join(cancerSiteSampleDF['unique_samples'])
#         
#         cancerSiteMutWeightDF['avg_mut_per_sample'] = cancerSiteMutWeightDF['sum'] / cancerSiteMutWeightDF['unique_samples']
#                 
#         cancerSiteMutWeightDF.to_csv('/Users/Rohil/Documents/Young Dawgs/' + study + '/cancer_site/cancerSiteMutationData.txt')    
# 
#         if plot == True:
#         
#             mutWeightPlot(cancerSiteMutWeightDF, study, 'cancer_site', 'avg_mut_per_sample')
# =============================================================================
