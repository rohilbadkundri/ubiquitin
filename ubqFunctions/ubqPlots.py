# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 08:10:33 2017

@author: Rohil
"""
import pandas as pd
import matplotlib.pyplot as plt
from ubqMutationData.ubqFunctions.ensureDirectory import ensureDirectory
from GSEA_ubiquitin_database.getGeneOntologies import ontologyDict
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

#category refers to cancer_type, ontology, cancer_site, etc
#subcat refers to specific subcategories, e.g. specific cancer site -- blood, brain, etc.
#focus refers to what you are analyzing, e.g. analyzing mutations in each gene within a specific cancer type -- gene would be the focus

exclude = ['base_substitution', 'small_indel', 'chromosomal_instability']

#add kwargs/args for colors
#try to generalize for mutSupertypeSumPlot
def mutPlot(df, study, category, cutoff, weighted):
    
    if weighted == True:
        metric = 'normalized mutation count'
        m = 'normalizedmut'
    else:
        metric = 'mutation count'
        m = 'mut'
    
    
    df = df.loc[:, df.columns.difference(exclude)]
    plot = df.drop(df.columns[df.apply(lambda col: col.sum() < cutoff)], axis = 1).plot(kind = 'bar', stacked = True, 
                               figsize = (35, 20), 
                               fontsize = 15,
                               color = colorList)

    plot.set_title(category + ' vs ' + metric, fontsize = 15)
    plot.set_ylabel(metric, fontsize = 15)
    plot.set_xlabel(category, fontsize = 15)

    plot.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 15})
    
    fig = plot.get_figure()

    ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/' + category + '/')
    
    fig.savefig('/Users/Rohil/Documents/Young Dawgs/' + study + '/' + category + '/' + category + '_' + m + '_image.png', bbox_inches='tight')
 
    plt.close()


def genePlot(df, study, category, subcat, cutoff, weighted):
    
    if weighted == True:
        metric = 'normalized mutation count'
        m = 'weight'
    else:
        metric = 'mutation count'
        m = ''
    
    plot = df.drop(df.columns[df.apply(lambda col: col.sum() < cutoff)], axis = 1)[:50].plot(kind = 'bar', stacked = True, 
                               figsize = (35, 20), 
                               fontsize = 15,
                               color = colorList)

    plot.set_title(category + '\n' + subcat + '\n' + 'gene vs ' + metric, fontsize = 15)
    plot.set_ylabel(metric, fontsize = 15)
    plot.set_xlabel('gene', fontsize = 15)

    plot.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 15})
    
    fig = plot.get_figure()

    ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/' + category + '/'  + category + '_gene' + m + '_images/')
    
    fig.savefig('/Users/Rohil/Documents/Young Dawgs/' + study + '/' + category + '/'  + category + '_gene' + m + '_images/' + subcat + '_gene' + m + '_image.png', bbox_inches='tight')
 
    plt.close()


def mutSupertypeSumPlot(df, study, category):
    
    mutSupertypePlot = df.reset_index().plot(x = category, y = ['base_substitution', 'small_indel', 'chromosomal_instability', 'fusion'],
                                kind = 'bar',
                                stacked = True,
                                figsize = (25,15),
                                title = category + ' vs #mutations',
                                fontsize = 15,
                                color = ['yellow', 'green', 'red', 'pink'])
            
    mutSupertypePlot.set_ylabel('#mutations')
    mutSupertypePlot.set_xlabel(category, size = 15)

    mutSupertypePlot.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 15})

    mutSupertypeFig = mutSupertypePlot.get_figure()

    ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/'  + category + '/')
    
    mutSupertypeFig.savefig('/Users/Rohil/Documents/Young Dawgs/' + study + '/'  + category + '/' + category + 'MutSupertypeSumData.png', bbox_inches='tight')

    plt.clf()
    

def mutSupertypeIndividualPlot(df, study, category):
    
    baseSubCols = df.columns[df.columns.isin(baseSubList)]
    
    smallIndelCols = df.columns[df.columns.isin(smallIndelList)]
    
    chromInstabCols = df.columns[df.columns.isin(chromInstabilityList)]
    
    supertypeDict = {'base_substitution' : baseSubCols, 'small_indel' : smallIndelCols, 'chromosomal_instability' : chromInstabCols}
    
    ensureDirectory('/Users/Rohil/Documents/Young Dawgs/' + study + '/' + category + '/')
    
    for key, value in supertypeDict.items():
    
        if df[value].empty:
            
            break
        
    else:
        
            df['sum'] = df[value].sum()
            df.sort_values(by = ['sum'], ascending = False, inplace = True)
            df.drop('sum', axis = 1, inplace = True)
                    
            supertypePlot = df.reset_index().plot(x = category, y = value, 
                              kind = 'bar',
                              stacked = True, 
                              figsize = (25, 15),
                              title = category + ' vs #' + key + ' mutations',
                              fontsize = 15)
            
            supertypePlot.set_ylabel('#mutations', size = 15)
            supertypePlot.set_xlabel(category, size = 15)
        
            supertypePlot.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 15})
        
            supertypeFig = supertypePlot.get_figure()
        
            supertypeFig.savefig('/Users/Rohil/Documents/Young Dawgs/' + study + '/' + category + '/' + category + '_' +  key +'.png', bbox_inches='tight')
        
            plt.clf()