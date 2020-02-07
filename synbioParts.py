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
import requests
import time
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
    
def doeGetSBOL(pfile='RefParts.csv', gfile='GeneParts.csv', 
               gsbol=None, libsize=32,
               getSequences=True,backtranslate=True, codontable='Eecoli.cut'):
    """
    Perform the DoE and generate the SBOL file from the 
    parts and genes files
    - RefParts.csv: Name, Type, Part
        Name: part name
        Type: origin, resistance, promoter, gene
        Part: SynBioHub URI
    - GeneParts.csv: Name, Type, Part, Step, Sequence
        Name: gene name
        Type: gene
        Step: Enzyme step in the pathway (eventually could be implemented 
                                          for the other genetic parts)
        Part: identifier (UniProt, etc) or SynBioHub URI
        Sequence: gene sequence (optional: if not given they are retrieved 
                                 from UniProt or SynBioHub)
    - Gsbol: SBOL file containing optimised versions of the genes (RBS, etc)
    """
    parts = pd.read_csv(pfile)
    genes = pd.read_csv(gfile)
    if gsbol is not None and os.path.exists(gsbol):
       genes = _readGenesSBOL(gsbol, genes)
    diagnostics, cons = getTheDoe(parts,genes,libsize)
    doc = getSBOL(parts,genes,cons,getSequences,backtranslate,codontable)
    diagnostics['sbol'] = doc.writeString()
    return diagnostics

def _readGenesSBOL(gsbol, genes, roles=[sbol.SO_GENE]):
    """ Create a new gene table containing the genes in the SBOL """
    partsdoc = sbol.Document()
    partsdoc.read(gsbol)
    gdict = {}
    for i in genes.index:
         gid = genes.loc[i,'Part']
         if gid not in gdict:
             gdict[gid] = []
         gdict[gid].append( genes.loc[i,:] )
    ngenes = []
    for part in partsdoc.componentDefinitions:
        if len( set(part.roles) & set(roles) ) > 0:
            dispid = part.displayId
            disp = dispid.split('_')
            rbs = disp[-2]
            cds = '_'.join( disp[0:(len(disp)-2)] )
            if cds in gdict:
                for row in gdict[cds]:
                    row = row.copy()
                    row.loc['Name'] = dispid
                    row.loc['Part'] = part.persistentIdentity
                    ngenes.append( row )
    ngenes = pd.DataFrame(ngenes, columns=genes.columns).sort_values(by=['Step','Name'])
    ngenes.index = np.arange(0, ngenes.shape[0])
    return ngenes
            
def _ReadParts(parts,registry='https://synbiohub.org'):
    """ A tabular csv file containing columns: Name, Type, Part is read 
    and each part in added to the SBOL doc.
    Part type should is overriden by the ontology definition in the registry """

    tf = parts
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

def _defineParts(doc,parts,getSequences=True,backtranslate=True,codontable='Eecoli.cut',registry='https://synbiohub.org',local=None):
    """ Parts is a dataframe containing columns: Name, Type, Part is read 
    and each part in added to the SBOL doc.
    Part type should is overriden by the ontology definition in the registry.
    If sequences should be provided in the SBOL doc, check if they are on the parts table
    or select if backtranslate from UniProt and provide codon table:
    (Eecoli.cut, Eyeast.cut, etc) See https://www.ebi.ac.uk/Tools/st/emboss_backtranseq/ """
    
    local = '/mnt/SBC1/data/ibisba/Genie/sbol_in.xml'
    locdoc = None
    if local is not None:
        locdoc = sbol.Document()
        locdoc.read(local)

    sboldef = {'promoter': sbol.SO_PROMOTER, 'gene': sbol.SO_CDS, 
               'origin': 'http://identifiers.org/so/SO:0000296', 
               'resistance': sbol.SO_CDS, 'terminator': sbol.SO_TERMINATOR}
    registry='https://synbiohub.org'
    repo = sbol.PartShop(registry)
    if locdoc is not None:
        for com in locdoc.componentDefinitions:
            if com.name is None:
                com.name = com.displayId
        locdoc.copy('http://liverpool.ac.uk',doc)
        
    for i in parts.index:
        name = parts.loc[i,'Name']
        ptype = parts.loc[i,'Type']
        part = parts.loc[i,'Part']
        if ptype in sboldef:
            if ptype == 'promoter':
                if getSequences:
                    try:
                        repo.pull(part,doc)
                        promoter = doc.getComponentDefinition(part)
                        promoter.name = name
                    except:
                        promoter = sbol.ComponentDefinition(name)
                else:
                    promoter = sbol.ComponentDefinition(name)                    
                promoter.roles = sboldef[ptype]
                promoter.setPropertyValue('http://purl.org/dc/terms/description',part)
                try:
                    if part not in doc.componentDefinitions:
                        promoter = doc.getComponentDefinition(part)
                except:
                    doc.addComponentDefinition(promoter)
                try:
                    terminator = sbol.ComponentDefinition('Ter')
                    terminator.roles = sboldef['terminator']
                    terminator.name = 'Ter'
                    doc.addComponentDefinition(terminator)       
                except:
                    continue
            elif ptype == 'origin':
                if getSequences:
                    try:
                        repo.pull(part,doc)
                        origin = doc.getComponentDefinition(part)
                    except:
                        origin = sbol.ComponentDefinition(name)
                else:
                        origin = sbol.ComponentDefinition(name)                    
                origin.roles = sboldef[ptype]
                origin.setPropertyValue('http://purl.org/dc/terms/description',part)
                origin.name = name
                try:
                    if part not in doc.componentDefinitions:
                        origin = doc.getComponentDefinition(part)
                except:
                    doc.addComponentDefinition(origin)
            elif ptype == 'gene' or ptype == 'resistance':
                gene = None
                # Case 1: where the gene is defined in the sbol 
                if name in doc.componentDefinitions:
                    gene = doc.componentDefinitions[name]
                try:
#                    gene = locdoc.getComponentDefinition(part)
                    gene = doc.componentDefinitions[name]
                except:
                    gene = None
                if gene is None:
                # Case 2: where the gene is in the repo
                    try:
                        repo.pull(part,doc)
                        gene = doc.getComponentDefinition(part)
                    except:
                        pass
                # Case 3: where the sequence is retrieved manually    
                if gene is None: 
                    gene = sbol.ComponentDefinition(name)
                    if getSequences:
                        # Get from the Sequence column if given (assume naseq)
                        if 'Sequence' in parts and not pd.isnan( parts.loc[i,'Sequence'] ):
                            seq = parts.loc[i,'Sequence']
                            seqdef = sbol.Sequence( name, seq, 'https://www.qmul.ac.uk/sbcs/iubmb/misc/naseq.html' )
                            gene.sequence = seqdef
                        else:
                            try:
                                query = 'https://www.uniprot.org/uniprot/' + name + '.fasta'
                                response = requests.get(query)
                                if backtranslate:
                                    resp = requests.post( 'https://www.ebi.ac.uk/Tools/services/rest/emboss_backtranseq/run', 
                                                         {'email': '@'.join(['pablo.carbonell','manchester.ac.uk']), 'title':'OptBioDes SBOL', 
                                                          'codontable': codontable, 'sequence': response.text} )
                                    jobid = resp.text
                                    status = 'RUNNING'
                                    while status == 'RUNNING':
                                        resp1  = requests.get('https://www.ebi.ac.uk/Tools/services/rest/emboss_backtranseq/status/'+jobid)
                                        status = resp1.text
                                        time.sleep(0.1)
                                    if status == 'FINISHED':
                                        resp2  = requests.get('https://www.ebi.ac.uk/Tools/services/rest/emboss_backtranseq/result/'+jobid+'/out')
                                        seq = ''.join( resp2.text.split('\n')[1:] )
                                        seqdef = sbol.Sequence( name, seq, 'https://www.qmul.ac.uk/sbcs/iubmb/misc/naseq.html' )
                                    else:
                                        seq = ''.join( response.text.split('\n')[1:] )                                   
                                        seqdef = sbol.Sequence( name, seq, 'https://www.qmul.ac.uk/sbcs/iupac/AminoAcid/' )
                                else:
                                    seq = ''.join( response.text.split('\n')[1:] )
                                    seqdef = sbol.Sequence( name, seq, 'https://www.qmul.ac.uk/sbcs/iupac/AminoAcid/' )
                                gene.sequence = seqdef
                            # If UniProt/EMBOSS was not succesful, leave the sequence empty
                            except:
                                pass
                gene.roles = sboldef[ptype]
                gene.setPropertyValue('http://purl.org/dc/terms/description',part)
                gene.name = name
                try:
                    doc.addComponentDefinition(gene)
                except:
                    pass
    return doc

def _definePartsSBOL(doc,parts,roles):
    """ Parts are given as an SBOL collection with the desired roles
    """
    partsdoc = sbol.doc()
    partsdoc.read(parts)
    for part in partsdoc.componentDefinitions:
        if len( set(part.roles) & set(roles) ) > 0:
            doc.addComponentDefinition( part )

def _definePartsOld(doc,parts):
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


def getTheDoe(parts, genes,size=32):
    diagnostics = getDoe(parts,genes,size)
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

    
def getSBOL(parts,genes,cons,getSequences=True,backtranslate=True,codontable='Eecoli.cut'):
    """
        parts: df with parts (see doeGetSBOL)
        genes: df with genes (see doeGetSBOL)
        cons: library of constructs generated by the DoE
        getSequences: retrieve sequences if not available in the gene list
        backtranslate: back translate protein sequences if not given
        codontable: only required for back translation
    """
    dic = {}
    for p in parts.index:
        dic[parts.loc[p,'Name']] = parts.loc[p,'Part']
    for g in genes.index:
        dic[genes.loc[g,'Name']] = genes.loc[g,'Part']
    namespace = "http://synbiochem.co.uk"
    sbol.setHomespace( namespace )
    doc = sbol.Document()
    doc = _defineParts(doc, parts)
    print('Parts defined')
    doc = _defineParts( doc, genes, getSequences, backtranslate, codontable )
    print('Genes defined')
    for row in np.arange(0,cons.shape[0]):
        plasmid = []
        for col in np.arange(0,cons.shape[1]):
            name = cons[row,col]
            if name == '-':
                continue
            try:
                component = doc.componentDefinitions[name]
            except:
                part = dic[name]
                component = doc.getComponentDefinition(part)
            if sbol.SO_PROMOTER in component.roles and col > 2:
                plasmid.append(doc.componentDefinitions['Ter'])
            plasmid.append(component)
        plasmid.append(doc.componentDefinitions['Ter'])
        # Create a new empty device named `my_device`
        plid = "plasmid%02d" % (row+1,)
        my_device = doc.componentDefinitions.create(plid)
        
        # Assemble the new device from the promoter, rbs, cds, and terminator from above.
        my_device.assemblePrimaryStructure(plasmid)
        
        # Set the role of the device with the Sequence Ontology term `gene`
        my_device.roles = sbol.SO_PLASMID
    return(doc)
    
    
