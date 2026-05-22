#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from .utilities_position_map import UtilitiesPositionMap
from .utilities_gos import UtilitiesGOs
from .utilities_protein_complex import UtilitiesProteinComplex
from .utilities_rnaseq import UtilitiesRNAseq
from FBD.fbd import FBD
import pandas as pd


#EXPLICAR. ES POR LA VIEJA VERSIÓN DE FBD??
def _empty_optional_dataset(dataset_name):
    if dataset_name == "gene_map_table":
        return pd.DataFrame(columns=["organism_abbreviation", "sequence_loc", "current_symbol"])

    if dataset_name == "dmel_unique_protein_isoforms":
        return pd.DataFrame(columns=["FB_gene_symbol"])

    return None


class DatasetLoader:
    """
    Thin acquisition layer for D3mel.

    It delegates dataset access to the public FBD API and only provides
    empty fallbacks for datasets that are optional in the current library.
    """

    def __init__(self):
        self.fbd = FBD()

    def load(self, dataset_name, optional=False):
        try:
            return self.fbd.download_file(dataset_name)
        except Exception:
            if optional:
                empty = _empty_optional_dataset(dataset_name)
                if empty is not None:
                    return empty
            raise


class Utilities(UtilitiesGOs, UtilitiesPositionMap,
                UtilitiesProteinComplex, UtilitiesRNAseq):

    def __init__(self):
        loader = DatasetLoader()

        go = loader.load("go-basic")
        gaf_raw = loader.load("gene_association")
        rpkm = loader.load("gene_rpkm_report")
        gene_map = loader.load("gene_map_table", optional=True)
        isoforms = loader.load("dmel_unique_protein_isoforms", optional=True)

        UtilitiesProteinComplex.__init__(self)
        UtilitiesGOs.__init__(self, go, gaf_raw)
        UtilitiesPositionMap.__init__(self, gene_map, isoforms)
        UtilitiesRNAseq.__init__(self, rpkm)

        self.gaf = UtilitiesGOs.get_GAF(self)
