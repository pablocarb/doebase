#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
synbioParts (c) University of Manchester 2019

synbioParts is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

Created on Fri May 31 13:38:19 2019

@author:  Pablo Carbonell, SYNBIOCHEM
@description: Routines to interface with synbio repositories 
    - parts: list of pointers to SynBiohub, we assume that each part 
    (origin, resistance, promoter) has been already registered in the part registry;
    - pathway: list of enzymes referred by either: 
        - SynbioHub identifier: useful if the part has been previously registered or
        for gene variants (RBS library or mutants)
        - UniProt identifier: a new empty container for this part is created in SynBioHub, 
        which will be later filled by the DNA design step
"""

import os
import re
import sbol
import numpy as np
import pandas as pd
from .OptDes import getDoe

def doeSBOL(pfile='RefParts.csv', gfile='GeneParts.csv', libsize=32, ofile='out.sbol'):
    """
    Perform the DoE and generate the SBOL file from the 
    parts and genes files
    - RefParts.csv: Name, Type, Part
    - GeneParts.csv: Name, Type, Part, Step
        Type: origin, resistance, promoter, gene 
        Step: Enzyme step in the pathway (eventually could be implemented 
        for the other genetic parts)
    """
    diagnostics, cons = getTheDoe(pfile,gfile,libsize)
    doc = getSBOL(pfile,gfile,cons)
    doc.write(ofile)
    
def doeGetSBOL(pfile='RefParts.csv', gfile='GeneParts.csv', libsize=32):
    """
    Perform the DoE and generate the SBOL file from the 
    parts and genes files
    - RefParts.csv: Name, Type, Part
    - GeneParts.csv: Name, Type, Part, Step
        Type: origin, resistance, promoter, gene 
        Step: Enzyme step in the pathway (eventually could be implemented 
        for the other genetic parts)
    """
    diagnostics, cons = getTheDoe(pfile,gfile,libsize)
    doc = getSBOL(pfile,gfile,cons)
    diagnostics['sbol'] = str(doc)
    return diagnostics
            
def _ReadParts(infile='RefParts.csv',registry='https://synbiohub.org'):
    """ A tabular csv file containing columns: Name, Type, Part is read 
    and each part in added to the SBOL doc.
    Part type should is overriden by the ontology definition in the registry """

    tf = pd.read_csv(infile)
    col = []
    doc = sbol.Document()
    for i in tf.index:
        namespace = "http://synbiochem.co.uk"
        sbol.setHomespace( namespace+'_'+str(i) )
        repo = sbol.PartShop(registry)
        name = tf.loc[i,'Name']
        ptype = tf.loc[i,'Type']
        part = tf.loc[i,'Part']
        if type(part) is str:
            print(part)
            repo.pull(part, doc)
    return doc

def _defineParts(doc,parts):
    """ Parts is a dataframe containing columns: Name, Type, Part is read 
    and each part in added to the SBOL doc.
    Part type should is overriden by the ontology definition in the registry """

    sboldef = {'promoter': sbol.SO_PROMOTER, 'gene': sbol.SO_CDS, 
               'origin': 'http://identifiers.org/so/SO:0000296', 
               'resistance': sbol.SO_CDS, 'terminator': sbol.SO_TERMINATOR}

    for i in parts.index:
        name = parts.loc[i,'Name']
        ptype = parts.loc[i,'Type']
        part = parts.loc[i,'Part']
        if ptype in sboldef:
            if ptype == 'promoter':
                promoter = sbol.ComponentDefinition(name)
                promoter.roles = sboldef[ptype]
                promoter.setPropertyValue('http://purl.org/dc/terms/description',part)
                doc.addComponentDefinition(promoter)
                try:
                    terminator = sbol.ComponentDefinition('Ter')
                    terminator.roles = sboldef['terminator']
                    doc.addComponentDefinition(terminator)       
                except:
                    continue
            elif ptype == 'origin':
                origin = sbol.ComponentDefinition(name)
                origin.roles = sboldef[ptype]
                origin.setPropertyValue('http://purl.org/dc/terms/description',part)
                doc.addComponentDefinition(origin)
            elif ptype == 'gene' or ptype == 'resistance':
                origin = sbol.ComponentDefinition(name)
                origin.roles = sboldef[ptype]
                origin.setPropertyValue('http://purl.org/dc/terms/description',part)
                doc.addComponentDefinition(origin)
    return doc


def getTheDoe(pfile='RefParts.csv', gfile='GeneParts.csv',size=32):
    diagnostics = getDoe(pfile,gfile,size)
    names = diagnostics['names']
    M = diagnostics['M']
    fact = diagnostics['fact']
    doemap = {}
    for n in names:
        x = int(re.findall('\d+',n)[0])
        doemap[x] = n
    cons = []
    for i in np.arange(0,M.shape[0]):
        z = 0
        construct = []
        for j in np.arange(1,int(max(fact))+1):
            if j in doemap:
                ### TO DO: find parts and genes and define the construct
                lev = M[i,z]
                z += 1
            else:
                lev = 0
            construct.append( fact[j].levels[lev])
        cons.append(construct)
    cons = np.array(cons)
    return diagnostics, cons

    
def getSBOL(pfile,gfile,cons):
    namespace = "http://synbiochem.co.uk"
    sbol.setHomespace( namespace )
    doc = sbol.Document()
    parts = pd.read_csv(pfile)
    doc = _defineParts(doc, parts)
    print('Parts defined')
    print(doc)
    genes = pd.read_csv(gfile)
    doc = _defineParts(doc, genes)
    print('Genes defined')
    print(doc)
    for row in np.arange(0,cons.shape[0]):
        plasmid = []
        for col in np.arange(0,cons.shape[1]):
            part = cons[row,col]
            if part == '-':
                continue
            component = doc.componentDefinitions[part]
            if sbol.SO_PROMOTER in component.roles and col > 2:
                plasmid.append('Ter')
            plasmid.append(part)
        plasmid.append('Ter')
        # Create a new empty device named `my_device`
        plid = "plasmid%02d" % (row+1,)
        my_device = doc.componentDefinitions.create(plid)

        # Assemble the new device from the promoter, rbs, cds, and terminator from above.
        my_device.assemblePrimaryStructure(plasmid)
        
        # Set the role of the device with the Sequence Ontology term `gene`
        my_device.roles = sbol.SO_PLASMID
        
    return(doc)
    
    
