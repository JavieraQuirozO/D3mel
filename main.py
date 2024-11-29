#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 19:26:03 2024

@author: javiicolors
"""

from D3mel import D3mel



user = "prueba2"
password = "prueba"
email = "prueba2@gmail.com"
   

#%%

d3m = D3mel(user, password)

#%%
GO = d3m.fbd.Ontology_Terms.GO()

#%%
GO = d3m.get_GO()
GOlist = ["GO:0006974", "GO:0006302", "GO:0006281", "GO:0006301", "GO:0000724", "GO:0006298", "GO:0000729", "GO:0006289", "GO:0006282", "GO:0006977"]

#%%
hijos = d3m.get_descendent(GOlist, i_max=1)

#%%

padres = d3m.get_ascendent(GOlist, i_max=1)
#%%
GOs_byCmp = d3m.get_GOs_byCmp()
gaf = d3m.get_GAF()
genes = list(gaf[gaf["GO ID"].isin(GOlist)]["DB Object ID"])
proyect = d3m.get_proyects()
#%%
dfs = d3m.RNA_by_proyect(genes, proyect[0:2])