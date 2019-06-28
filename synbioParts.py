#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
synbioParts (c) University of Manchester 2019

synbioParts is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

Created on Fri May 31 13:38:19 2019

@author:  Pablo Carbonell, SYNBIOCHEM
@description: Routines to interface with synbio repositories 
"""

import os
import re
import sbol
import numpy as np
import pandas as pd
from .doebase import doeTemplate, promoterList, plasmidList, read_excel
from .OptDes import callDoE

def combine(parts, pathway):
    """ This routine combines together the parts and the coding regions of the pathway.
    - parts: list of pointers to SynBiohub, we assume that each part 
    (origin, resistance, promoter) has been already registered in the part registry;
    - pathway: list of enzymes referred by either: 
        - SynbioHub identifier: useful if the part has been previously registered or
        for gene variants (RBS library or mutants)
        - UniProt identifier: a new empty container for this part is created in SynBioHub, 
        which will be later filled by the DNA design step """
        
def ReadParts(infile='RefParts.csv',registry='https://synbiohub.org'):
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

def defineParts(doc,parts):
    """ Parts is a dataframe containing columns: Name, Type, Part is read 
    and each part in added to the SBOL doc.
    Part type should is overriden by the ontology definition in the registry """

    sboldef = {'promoter': sbol.SO_PROMOTER, 'gene': sbol.SO_GENE, 
               'origin': sbol.SO_PLASMID, 'resistance': sbol.SO_GENE}

    for i in parts.index:
        name = parts.loc[i,'Name']
        ptype = parts.loc[i,'Type']
        part = parts.loc[i,'Part']
        if ptype in sboldef:
            if ptype == 'promoter':
                promoter = sbol.ComponentDefinition(name)
                promoter.roles = sbol.SO_PROMOTER
                promoter.setPropertyValue('http://purl.org/dc/terms/description',part)
                doc.addComponentDefinition(promoter)
            elif ptype == 'origin':
                origin = sbol.ComponentDefinition(name)
                origin.roles = sbol.SO_PLASMID
                origin.setPropertyValue('http://purl.org/dc/terms/description',part)
                doc.addComponentDefinition(origin)
            elif ptype == 'gene' or ptype == 'resistance':
                origin = sbol.ComponentDefinition(name)
                origin.roles = sbol.SO_GENE
                origin.setPropertyValue('http://purl.org/dc/terms/description',part)
                doc.addComponentDefinition(origin)
                
    return doc

def defineTemplate(pfile='RefParts.csv', gfile='GeneParts.csv'):
    parts = pd.read_csv(pfile)
    genes = pd.read_csv(gfile)
    prom = []
    ori = []
    for i in parts.index:
        ptype = parts.loc[i,'Type']
        name = parts.loc[i,'Name']
        if ptype == 'promoter':
            prom.append(name)
        elif ptype == 'origin':
            ori.append(name)
    for i in range(0,len(prom)):
        prom.append(None)
    tree = []
    gdict = {}
    for i in genes.index:
        name = genes.loc[i,'Name']
        step = "gene%00d" % (int(genes.loc[i,'Step']),)
        tree.append(step)
        if step not in gdict:
            gdict[step] = []
        gdict[step].append(name)
    doe = doeTemplate(tree, origins=ori, promoters=prom, genes=gdict, positional=False)
    return doe, parts, genes

def getConstruct(doe,diagnostics,parts,genes):
    import pdb
    pdb.set_trace()


def getDoe(pfile='RefParts.csv', gfile='GeneParts.csv',size=32):
    doe,parts,genes = defineTemplate(pfile=pfile, gfile=gfile)
    fact, partinfo = read_excel( None, doedf=doe )
    seed = np.random.randint(10000)
    diagnostics = callDoE(fact, size, seed=seed)
    names = diagnostics['names']
    M = diagnostics['M']
    doemap = {}
    for n in names:
        x = int(re.findall('\d+',n)[0])
        doemap[x] = n
    cons = []
    for i in np.arange(0,M.shape[0]):
        z = 0
        construct = []
        for j in np.arange(1,int(max(doe['DoE position']))+1):
            if j in doemap:
                ### TO DO: find parts and genes and define the construct
                lev = M[j-1,z]
                z += 1
            else:
                lev = 0
            construct.append( fact[j].levels[lev])
        cons.append(construct)
    cons = np.array(cons)
    return doe, diagnostics, cons

    
def getSBOL(pfile,gfile,cons):
    namespace = "http://synbiochem.co.uk"
    sbol.setHomespace( namespace )
    doc = sbol.Document()
    parts = pd.read_csv(pfile)
    doc = defineParts(doc, parts)
    print('Parts defined')
    print(doc)
    genes = pd.read_csv(gfile)
    doc = defineParts(doc, genes)
    print('Genes defined')
    print(doc)
    for row in np.arange(0,cons.shape[0]):
        plasmid = []
        for col in np.arange(0,cons.shape[1]):
            part = cons[row,col]
            if part == '-':
                continue
            component = doc.componentDefinitions[part]
            plasmid.append(part)
        # Create a new empty device named `my_device`
        plid = "plasmid%02d" % (row+1,)
        my_device = doc.componentDefinitions.create(plid)

        # Assemble the new device from the promoter, rbs, cds, and terminator from above.
        my_device.assemblePrimaryStructure(plasmid)
        
        # Set the role of the device with the Sequence Ontology term `gene`
        my_device.roles = sbol.SO_PLASMID
        
    return(doc)
    
    
def example(pfile='RefParts.csv', gfile='GeneParts.csv', ofile='out.sbol'):
#    pwd = os.path.dirname(os.path.realpath(__file__))
#    pfile = os.path.join(pwd,'RefParts.csv')
#    gfile = os.path.join(pwd,'GeneParts.csv')
    doe, diagnostics, cons = getDoe(pfile,gfile,32)
    doc = getSBOL(pfile,gfile,cons)
    doc.write(ofile)
    
# Set default namespace
if __name__ == '__main__':    
    namespace = "http://synbiochem.co.uk"
    sbol.setHomespace( namespace )
    doc = sbol.Document()
    pfile = 'RefParts.csv'
    gfile = 'GeneParts.csv'
    parts = pd.read_csv(pfile)
    doc = defineParts(doc, parts)
    print('Parts defined')
    print(doc)
    genes = pd.read_csv(gfile)
    doc = defineParts(doc, genes)
    print('Genes defined')
    print(doc)
    doe = defineTemplate(parts,genes)
    
