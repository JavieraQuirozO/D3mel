#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 19:26:03 2024

@author: javiicolors
"""

from D3mel import D3mel
import random



#%%

d3m = D3mel()

#%%
GO = d3m.get_GO()
GOlist = ["GO:0006793", "GO:0016310", "GO:0006796", "GO:0006470", "GO:0046835"] #metabolismo del fosforo

#%%
hijos_nietos = d3m.get_descendent(GOlist, i_max=2)

#%%
padres_abuelos = d3m.get_ascendent(GOlist, i_max=2)
#%%
elementos = GOlist

for dic in [hijos_nietos, padres_abuelos]:
    for key, values in dic.items():
        elementos.extend(values)
elementos = set(elementos)

#%%
negativos = list(set(GO.nodes()) - elementos)
negativos = random.sample(negativos, len(elementos))

#%%
GOs_byCmp = d3m.get_GOs_byCmp()
gaf = d3m.get_GAF()

#%%
proyects = d3m.get_proyects()
genes_positivos = list(gaf[gaf["GO ID"].isin(GOlist)]["DB Object ID"])
conj_positivo = d3m.RNA_by_proyect(genes_positivos, proyects)

#%%
genes_negativos = list(gaf[gaf["GO ID"].isin(negativos)]["DB Object ID"])
conj_negativos = d3m.RNA_by_proyect(genes_negativos, proyects)