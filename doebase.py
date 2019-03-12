#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OptExp (c) University of Manchester 2019

OptExp is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>

Created on Tue Nov 27 16:01:46 2018

@author: Pablo Carbonell
@description: Basic routines common to several doe scripts.
"""
import pandas as pd
import numpy as np

def read_excel(e, s=1, doedf=None):
    """ Reads a DoE sheet """
    if e is None:
        df = doedf
    else:
        df = pd.read_excel(e)
    fact = {}
    partinfo = {}
    offset = 1
    fcol = 0
    r = 0
    while fcol is not None:
        r += 1
        fcol = None
        try:
            fcol = df.iloc[r-1, offset-1]
            factor = int( fcol )
            positional = df.iloc[r-1, offset]
            component = str( df.iloc[r-1, offset+1] )
            part = str( df.iloc[r-1, offset+2] )
            if positional != positional:
                positional = None
        except:
            continue
        if part is None:
            if factor in fact:
                i = len(fact[factor]['levels'])+1
            else:
                i = 1
            part = 'P'+str(factor)+'_'+str(i)

        if factor not in fact:
            fact[factor] = {
                    'positional': positional,
                    'component': component,
                    'levels': []
            }
        if part == 'blank':
            part = None
        fact[factor]['levels'].append(part)
    return fact, partinfo


def uniformData(tree, promoters, positional):
    """ Uniformize the format of the data in order
    to allow polymorphic inputs """
    if type(positional) is bool:
        if positional:
            positional = tree
        else:
            positional = []
    if type(promoters) is not dict:
        promoters = {step:promoters for step in tree}        
    return promoters, positional

def entry( info, columns=['DoE position', 'DoE permutation class', 'DoE designation', 
                   'Part number', 'Type', 'Name', 
                   'Info', 'Comments', 'Enzyme'] ):
    """ Use pandas Series for each entry in the sheet 
    """   
    labels = {
     'origin': ('origin', 'Replication origin'), 
     'resistance': ('resistance', 'Antibiotic resistance'),
     'backbone': ('promoter', 'Backbone promoter'), 
     'gene': ('gene', 'ORF'),
     'promoter': ('promoter', 'promoter')
     }   
    row = pd.Series( None, index=columns )
    if info in labels:
        row[ 'DoE designation' ] = labels[info][0]
        row[ 'Type' ] = labels[info][1]
    return row

def doeTemplate(tree,
                origins=['ori1','ori2'],
                promoters=['prom1', 'prom2', None, None],
                genes = {},
                positional=True):
    """ Generates the DoE sheet:
        - tree: ordered list of reaction steps.
        - origins: list of origins of replication (plasmid copy number)
        - promoters: 
            - list of promoters (Nones at the end for ratio of no-promoters)
            - dictionary with list of promoters for each step
        - genes: dictionary containing the genes for each step
        - positional: gene permutation
            - Boolean: for all genes;
            - List of steps for permutation
    """    
    promoters, positional = uniformData( tree, promoters, positional )
    
    origin = entry('origin')
    resistance = entry('resistance')
    bbpromoter = entry('backbone')
    orf = entry('gene')
    promoter = entry('promoter')
    pos = 1

    doe = pd.DataFrame( )  
    for ori in origins:
        origin.Name = ori
        origin['Part number'] = origin.Name
        origin['DoE position'] = pos
        doe = doe.append( origin, ignore_index=True, sort=False )
    pos += 1
    resistance['Name'] = 'res1'
    resistance['Part number'] = resistance.Name
    resistance['DoE position'] = pos
    doe = doe.append( resistance, ignore_index=True )

    first = True
    for rid in tree:
        pos += 1
        for p in promoters[rid]:
            if first:
                if p is None:
                    continue
                else:
                    bbpromoter['DoE position'] = pos
                    bbpromoter['Name'] = p
                    bbpromoter['Part number'] = bbpromoter.Name
                    doe = doe.append( bbpromoter, ignore_index=True )
            else:        
                promoter['DoE position'] = pos
                promoter.Name = p
                if promoter['Name'] is not None:
                    promoter['Part number'] = promoter.Name
                else:
                    promoter['Part number'] = '-'
                doe = doe.append( promoter, ignore_index=True )
        pos += 1
        if rid not in genes:
            genes[rid] = rid
        for gene in genes[rid]:
            orf['DoE position'] = pos
            orf['Name'] = gene
            orf['Part number'] = orf.Name
            if rid in positional:
                orf['DoE permutation class'] = 1
            else:
                orf['DoE permutation class'] = None
            doe = doe.append( orf, ignore_index=True )
        first = False 
    doe = doe.reindex(origin.index, axis=1) 
    return doe

def promoterList( levels, ratio=0.5 ):
    """ A simple list of promoters based on number of levels.
        Ratio is the percentage of empty promoters respect to full promoters"""
    plist = []
    for i in np.arange( levels ):
        plist.append( 'prom%d' % (i+1,) )
    empty = np.ceil( levels*ratio/(1-ratio) )
    for i in np.arange( empty ):
        plist.append( None )
    return plist

def plasmidList( levels ):
    """ A simple list of plasmids based on the number of levels. """
    plist = []
    for i in np.arange( levels ):
        plist.append( 'pl%d' % (i+1,) )
    return plist
