# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 19:19:56 2017

@author: Rohil
"""
import pandas as pd
import matplotlib.pyplot as plt
from ubqFunctions.ensureDirectory import ensureDirectory
from matplotlib import style
style.use('ggplot')

#maybe try to key colors to mutation types
colorList = ['yellow', '#A3E4D7', '#5DADE2', '#16A085', '#0B5345',  'blue', '#D2B4DE', '#FCF3CF', '#F7DC6F', '#E67E22', '#922B21', 'red', '#9B59B6', '#633974', '#A04000', 'magenta']

baseSubList = ['missense_mutation', 'nonsense_mutation','nonstop_mutation', 'splice_site', 'translation_start_site', 'nonsynonymous']
smallIndelList = ['frame_shift_ins', 'frame_shift_del', 'in_frame_ins', 'in_frame_del']
chromInstabilityList = ['homozygous_del', 'hemizygous_del', 'gain', 'high_lvl_amplification']

cutoffNum = 10
weightCutoffNum = .1
geneCutoffNum = 5

def mutPlot(df, study, category):
    
    exclude = ['base_substitution', 'small_indel', 'chromosomal_instability']
    
    
    mutPlot = df.loc[:, df.columns.difference(exclude)].drop(df.columns[df.apply(lambda col: col.sum() < cutoffNum)], axis = 1).plot(kind = 'bar',
                       stacked = True, 
                       figsize = (25, 15),
                       title = category + ' vs #mutations',
                       fontsize = 15,
                       color = colorList)
    
    mutPlot.set_ylabel('#mutations', size = 15)
    mutPlot.set_xlabel(category, size = 15)

    mutPlot.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 15})

    mutFig = mutPlot.get_figure()

    ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/'  + category + '/')
    
    mutFig.savefig('/Users/Rohil/Documents/Young Dawgs/' + study + '/'  + category + '/' + category + '_mut_data.png', bbox_inches='tight')

    plt.clf()
    
    

def genePlot(df, study, category, subcat, focus, cutoff):

    genePlot = df.drop(df.columns[df.apply(lambda col: col.sum() < cutoff)], axis = 1)[:50].plot(kind = 'bar', stacked = True, 
                               figsize = (35, 20), 
                               fontsize = 15,
                               color = colorList)

    genePlot.set_title(category + '\n' + i + '\ngene vs #mutations', fontsize = 15)
    genePlot.set_ylabel('#mutations', fontsize = 15)
    genePlot.set_xlabel('gene', fontsize = 15)

    genePlot.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 15})
    
    geneFig = genePlot.get_figure()

    ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/' + category + '/'  + category + '_gene_images/')
    
    geneFig.savefig('/Users/Rohil/Documents/Young Dawgs/' + study + '/' + category + '/'  + category + '_gene_images/' + i + '(GENE).png', bbox_inches='tight')
 
    plt.close()
    

def ontologyPlot(df, study, category, i):

    ontologyPlot = df.drop(df.columns[df.apply(lambda col: col.sum() < cutoffNum)], axis = 1).plot(kind = 'bar', stacked = True,
                                                                                                          figsize = (35, 15),
                                                                                                          title = category + '\n' + i + '\nontology vs # mutations',
                                                                                                          fontsize = 15,
                                                                                                          color = colorList)
    
    ontologyPlot.set_ylabel(i, fontsize = 15)
    ontologyPlot.set_xlabel('#mutations', fontsize = 15)

    ontologyPlot.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 15})
    
    ontologyFig = ontologyPlot.get_figure()
    
    ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/'  + category + '/' + category + '_ontology_images/')
    
    ontologyFig.savefig('/Users/Rohil/Documents/Young Dawgs/' + study + '/'  + category + '/' + category + '_ontology_images/' + i + 'OntologyData.png', bbox_inches='tight')    
    
    plt.clf()


def mutWeightPlot(df, study, category, yy):
                
    #df.sort_values(by = yy, ascending = False, inplace = True)
    
    #add a line to sort by sum here if needed
    
    mutWeightPlot = df.drop(df.columns[df.apply(lambda col: col.sum() < weightCutoffNum)], axis = 1)[:50].plot(kind = 'bar', stacked = True,
                                                                                                          figsize = (35, 15),
                                                                                                          title = category + ' vs ' + yy,
                                                                                                          fontsize = 15,
                                                                                                          color = colorList)
    
    mutWeightPlot.set_ylabel(yy, fontsize = 15)
    mutWeightPlot.set_xlabel(category, fontsize = 15)

    mutWeightPlot.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 15})
    
    mutWeightFig = mutWeightPlot.get_figure()
    
    ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/'  + category + '/')
    
    mutWeightFig.savefig('/Users/Rohil/Documents/Young Dawgs/' + study + '/'  + category + '/' + category + 'MutWeightData.png', bbox_inches='tight')    
    
    plt.clf()
    


def geneWeightPlot(df, study, category, i):
        
    geneWeightPlot = df.drop(df.columns[df.apply(lambda col: col.sum() < weightCutoffNum)], axis = 1)[:30].plot(kind = 'bar', stacked = True,
                             figsize = (25, 15),
                             title = i + '\ngene vs %mutations',
                             fontsize = 15,
                             color = colorList)
    
    geneWeightPlot.set_ylabel('%mutations', fontsize = 15)
    geneWeightPlot.set_xlabel('gene', fontsize = 15)

    geneWeightPlot.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 15})
    
    geneWeightFig = geneWeightPlot.get_figure()
        
    ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/' + category + '/'  + category + '_gene_weight_images/')
        
    geneWeightFig.savefig('/Users/Rohil/Documents/Young Dawgs/' + study + '/' + category + '/'  + category + '_gene_weight_images/' + i + '(gene_weight).png', bbox_inches='tight')
    
    plt.close()
    

def ontologyWeightPlot(df, study, category, i):
    
    ontologyWeightPlot = df.drop(df.columns[df.apply(lambda col: col.sum() < weightCutoffNum)], axis = 1).plot(kind = 'bar', stacked = True,
                             figsize = (35, 25),
                             title = i + '\nontology vs avg mut per gene',
                             fontsize = 15,
                             color = colorList)
    
    ontologyWeightPlot.set_ylabel('%avg mut per gene', fontsize = 15)
    ontologyWeightPlot.set_xlabel('ontology', fontsize = 15)

    ontologyWeightPlot.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 15})
    
    ontologyWeightFig = ontologyWeightPlot.get_figure()
        
    ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/' + category + '/'  + category + '_ontology_weight_images/')
        
    ontologyWeightFig.savefig('/Users/Rohil/Documents/Young Dawgs/' + study + '/' + category + '/'  + category + '_ontology_weight_images/' + i + '(ontology_weight).png', bbox_inches='tight')
    
    plt.close()
