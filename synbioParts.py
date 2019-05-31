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

import sbol

doc = sbol.Document()
repo = sbol.PartShop('https://synbiohub.org')
repo.pull('http://synbiohub.org/public/igem/BBa_R0010/1', doc)