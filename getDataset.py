#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 16:11:54 2026

@author: javiicolors
"""

from collections import Counter
import random
import re

from D3mel import D3mel


GO_PATTERN = re.compile(r"^GO:[0-9]{7}$")


class Get_dataset():

  def __init__(self, DB_Object_Type = "protein"):

    self.d3m = D3mel()
    self.GO = self.d3m.get_GO()
    self.gaf = self.gaf_ext(DB_Object_Type)

  def _go_terms(self, value):
    if isinstance(value, list):
      terms = value
    elif isinstance(value, str):
      terms = [value]
    else:
      return []

    return [
      term for term in terms
      if isinstance(term, str) and GO_PATTERN.match(term)
    ]

  def gaf_ext(self, DB_Object_Type):
    def concatenar(val1, val2):
        terms = []
        if isinstance(val2, list):
            terms.extend(val2)
        terms.append(val1)
        return sorted(set(self._go_terms(terms)))

    gaf = self.d3m.get_GAF()
    gaf = gaf[gaf['DB Object Type'] == DB_Object_Type]
    gos = list(gaf['GO ID'].drop_duplicates())
    prop = self.d3m.get_ascendent(gos)
    gaf['Propagacion'] = gaf['GO ID'].map(prop).fillna(0)
    gaf['GO ID'] = gaf.apply(lambda row: concatenar(row['GO ID'], row['Propagacion']), axis=1)
    return gaf.iloc[:,:-1]


  def get_go_count(self, i_min = 50, i_max = 1000):

    gaf = self.gaf
    resultado = []
    for _, gene_rows in gaf.groupby("DB Object ID"):
        gos = gene_rows["GO ID"].to_list()
        aux = set()
        for go_terms in gos:
          aux.update(self._go_terms(go_terms))
        resultado.extend(list(aux))

    frecuencia = Counter(resultado)
    return {key: count for key, count in frecuencia.items() if i_min <= count <= i_max}


  def positives(self, GOlist):
    if isinstance(GOlist, str):
      return set(self.gaf.loc[self.gaf['GO ID'].apply(lambda x: GOlist in self._go_terms(x)), 'DB Object ID'].tolist())
    elif isinstance(GOlist, list):
      requested = set(GOlist)
      return set(self.gaf.loc[self.gaf['GO ID'].apply(lambda x: bool(requested & set(self._go_terms(x)))), 'DB Object ID'].tolist())
    else:
      print("Ingrese el termino GO en formato string o lista de strings")
      return None


  def _rng(self, seed = None):
    return random.Random(seed) if seed is not None else random


  def RU(self, GOlist, proyects = None, seed = None):
    if proyects == None:
      proyects = self.d3m.get_proyects()

    genes_positivos = self.positives(GOlist)
    if genes_positivos is None:
      return [None, None]

    genes_negativos = set(self.gaf.loc[~self.gaf['DB Object ID'].isin(genes_positivos), 'DB Object ID'].tolist())
    if len(list(genes_positivos)) < len(list(genes_negativos)):
      genes_negativos = self._rng(seed).sample(sorted(genes_negativos), k = len(list(genes_positivos)))
    conj_positivo = self.d3m.RNA_by_proyect(sorted(genes_positivos), proyects)
    conj_negativos = self.d3m.RNA_by_proyect(sorted(genes_negativos), proyects)

    return conj_positivo, conj_negativos

  def SP(self, GOlist, proyects = None, seed = None):
      if isinstance(GOlist, str):
          GOlist = [GOlist]

      genes_positivos = self.positives(GOlist)

      if genes_positivos is None:
        return [None, None]

      if proyects == None:
        proyects = self.d3m.get_proyects()

      padres_dict = self.d3m.get_descendent(GOlist, i_max = 1)
      padres = set(GOlist + [v for values in padres_dict.values() for v in values])
      hijos_dict = self.d3m.get_ascendent(list(padres), i_max=1)
      hijos = set([v for values in hijos_dict.values() for v in values if v and v not in GOlist])

      genes_negativos = set(self.gaf[self.gaf["GO ID"].apply(lambda x: any(go in x for go in hijos))]["DB Object ID"])
      genes_negativos = genes_negativos - genes_positivos
      if len(list(genes_positivos)) < len(list(genes_negativos)):
          genes_negativos = self._rng(seed).sample(sorted(genes_negativos), k = len(list(genes_positivos)))
      conj_positivo = self.d3m.RNA_by_proyect(sorted(genes_positivos), proyects)
      conj_negativos = self.d3m.RNA_by_proyect(sorted(genes_negativos), proyects)


      return conj_positivo, conj_negativos
