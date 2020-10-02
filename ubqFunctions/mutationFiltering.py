# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 08:58:48 2017

@author: Rohil
"""
import pandas as pd

baseSubList = ['missense_mutation', 'nonsense_mutation','nonstop_mutation', 'splice_site', 'translation_start_site', 'nonsynonymous']
smallIndelList = ['frame_shift_ins', 'frame_shift_del', 'in_frame_ins', 'in_frame_del']
chromInstabilityList = ['homozygous_del', 'hemizygous_del', 'gain', 'high_lvl_amplification']

def combineRepeatMutations(df):
    
    df['mutation_type'].replace(to_replace = 'nonframeshift[.,_,-]deletion', value = 'in_frame_del', regex = True, inplace = True) #working
    #df['mutation_type'].replace(to_replace = 'splice.*', value = 'splice_site', regex = True, inplace = True) #not working, fixed - need to test #don't need anyways
    #df['mutation_type'].replace(to_replace = 'missense.*', value = 'missense_mutation', regex = True, inplace = True) #not working, fixed - need to test #don't need anyways
    df.loc[df['mutation_type'].str.contains('missense'), 'mutation_type'] = 'missense_mutation' #working
    df.loc[df['mutation_type'].str.contains('splice'), 'mutation_type'] = 'splice_site' #working
    df['mutation_type'].replace(to_replace = 'stop_gain', value= 'nonsense_mutation', inplace = True) #working
    
    return df
    

def getMutationSupertype(df):
   
    baseSubCols = df.columns[df.columns.isin(baseSubList)]
    df['base_substitution'] = df[baseSubCols].sum(axis = 1)

    smallIndelCols = df.columns[df.columns.isin(smallIndelList)]
    df['small_indel'] = df[smallIndelCols].sum(axis = 1)
    
    chromInstabCols = df.columns[df.columns.isin(chromInstabilityList)]
    df['chromosomal_instability'] = df[chromInstabCols].sum(axis = 1)
    
    if 'fusion' not in df:
        
        df['fusion'] = 0
    
    return df