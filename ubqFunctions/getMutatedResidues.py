# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 22:00:39 2018

@author: Rohil
"""
import pandas as pd
import regex as re

def getMutRes(mutDF):
    
    missense = mutDF[mutDF.mutation_type.str.match('missense')][mutDF.amino_acid_change.str.match('^[A-Z]\d{1,4}[A-Z]$')]
    missense['wild_res'] = missense.amino_acid_change.apply(lambda x: x[0])
    missense['mut_res'] = missense.amino_acid_change.apply(lambda x: x[-1])
    missense['mut_pos'] = missense.amino_acid_change.replace(regex=True, to_replace= '\D', value= '')
    #print ('missense: ' + str(missense.shape))
    
    splice1 = mutDF[mutDF.mutation_type.str.match('splice')][mutDF.amino_acid_change.str.match('^X\d{1,4}_splice$')] 
    splice1['wild_res'] = splice1.amino_acid_change.replace(regex=True, to_replace= '[^A-Z]', value= '')
    splice1['mut_res'] = 'splice'
    splice1['mut_pos'] = splice1.amino_acid_change.replace(regex=True, to_replace= '\D', value= '')
    
    splice2 = mutDF[mutDF.mutation_type.str.match('splice')][mutDF.amino_acid_change.str.match('^[A-Z]\d{1,4}[A-Z]$')]
    splice2['wild_res'] = splice2.amino_acid_change.apply(lambda x: x[0])
    splice2['mut_res'] = splice2.amino_acid_change.apply(lambda x: x[-1])
    splice2['mut_pos'] = splice2.amino_acid_change.replace(regex=True, to_replace= '\D', value= '')
    
    splice3 = mutDF[mutDF.mutation_type.str.match('splice')][mutDF.amino_acid_change.str.match('^[A-Z]\d{1,4}=$')]
    splice3['wild_res'] = splice3.amino_acid_change.apply(lambda x: x[0])
    splice3['mut_res'] = splice3.wild_res
    splice3['mut_pos'] = splice3.amino_acid_change.replace(regex=True, to_replace= '\D', value= '')
    #print ('splice: ' + str(splice1.shape[0]+splice2.shape[0]))
    
    splice4 = mutDF[mutDF.mutation_type.str.match('splice')][mutDF.amino_acid_change.str.match('^[A-Z]\d{1,4}$')]
    splice4['wild_res'] = splice4.amino_acid_change.apply(lambda x: x[0])
    splice4['mut_res'] = splice4.wild_res
    splice4['mut_pos'] = splice4.amino_acid_change.replace(regex=True, to_replace= '\D', value= '')
    #print ('splice: ' + str(splice1.shape[0]+splice2.shape[0]))
    
    
    translation = mutDF[mutDF.mutation_type.str.match('translation')]
    translation['wild_res'] = 'M'
    translation['mut_res'] = 'translation'
    translation['mut_pos'] = 1
    #print ('translation: ' + str(translation.shape))
    
    nonstop = mutDF[mutDF.mutation_type.str.match('nonstop')][mutDF.amino_acid_change.str.match('^\*\d{1,4}[A-Z]ext\*')]
    nonstop['wild_res'] = 'stop'
    nonstop['mut_res'] = nonstop.amino_acid_change.apply(lambda x: x.split('ext')[0][-1])
    nonstop['mut_pos'] = nonstop.amino_acid_change.apply(lambda x: re.findall('^\*\d{1,4}[A-Z]ext',x)[0]).str.lstrip('*').str.replace('[A-Za-z]', '')            
    
    #print ('nonstop: ' + str(nonstop.shape))
    
    #inf_del1 = mutDF[mutDF.mutation_type.str.match('in_frame_del')][mutDF.amino_acid_change.str.match('[A-Z]\d{1,4}del')]
    #inf_del1['wild_res'] = inf_del1.amino_acid_change.apply(lambda x: x[0])
    #inf_del1['mut_res'] = 'del'
    #inf_del1['mut_pos'] = inf_del1.amino_acid_change.replace(regex=True, to_replace= '\D', value= '')
    
    #inf_del2 = mutDF[mutDF.mutation_type.str.match('in_frame_del')][mutDF.amino_acid_change.str.match('[A-Z]\d{1,4}_[A-Z]\d{1,4}del')]
    #inf_del2_res = inf_del2.apply(get_inf_del2_res, axis=1)
    #inf_del2_res.columns = ['wild_res', 'mut_res']
    #inf_del2 = pd.concat([inf_del2, inf_del2_res], axis=1)
    #print('inf_del: ' + str(inf_del1.shape[0]+inf_del2.shape[0]))
    
    nonsense = mutDF[mutDF.mutation_type.str.match('nonsense')][mutDF.amino_acid_change.str.match('^[A-Z]\d{1,4}\*$')]
    nonsense['wild_res'] = nonsense.amino_acid_change.apply(lambda x: x[0])
    nonsense['mut_res'] = 'stop'
    nonsense['mut_pos'] = nonsense.amino_acid_change.replace(regex=True, to_replace= '\D', value= '')
    #print('nonsense: ' + str(nonsense.shape))
    
    frameshift1 = mutDF[mutDF.mutation_type.str.match('frame_shift')][mutDF.amino_acid_change.str.match('^[A-Z]\d{1,4}[A-Z]fs\*')]
    frameshift1['wild_res'] = frameshift1.amino_acid_change.apply(lambda x: x[0])
    frameshift1['mut_res'] = 'fs'
    frameshift1['mut_pos'] = frameshift1.amino_acid_change.apply(lambda x: x.split('fs')[0]).replace(regex=True, to_replace= '\D', value= '')
    
    frameshift2 = mutDF[mutDF.mutation_type.str.match('frame_shift')][mutDF.amino_acid_change.str.match('^[A-Z]\d{1,4}\*')]
    frameshift2['wild_res'] = frameshift2.amino_acid_change.apply(lambda x: x[0])
    frameshift2['mut_res'] = 'fs'
    frameshift2['mut_pos'] = frameshift2.amino_acid_change.replace(regex=True, to_replace= '\D', value= '')
    #print ('fs: ' + str(frameshift1.shape[0] + frameshift2.shape[0]))
    
    return (pd.concat([missense, splice1, splice2, splice3, splice4, translation, nonstop, nonsense, frameshift1, frameshift2], axis=0))