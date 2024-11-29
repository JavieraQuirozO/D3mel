# -*- coding: utf-8 -*-

import ast

#%%

class UtilitiesPositionMap():
    
    def __init__(self, fbd):
        self.fbd = fbd
    
    def coding_genes(self):
        proteins = self.fbd.Genes.Unique_protein_isoforms()
        return list(set(proteins['FB_gene_symbol']))
        
    
    def get_location(self):
        map_loc = self.fbd.Genes.Genes_map()

        #%%
        map_loc = map_loc[map_loc['##organism_abbreviation'] == "Dmel"]
        map_loc.dropna(subset=['sequence_loc'], inplace = True) #Archivo oficial y actualizado de las pocisiones

        #%%

        map_loc[['Chr', 'Start', 'End']] = map_loc['sequence_loc'].str.extract(r'(\w+):(\d+)\.\.(\d+).*')
        map_loc.index = map_loc['current_symbol']
        prot = list(set(self.coding_genes()) & set(map_loc.index)) 
        loc_map = map_loc.loc[prot, ['Chr', 'Start', 'End']].copy()
        return loc_map

    def get_order(self, valid_chr = ['2L', '2R', '3L', '3R', 'X']):
        loc_df = self.get_location()
        loc_df = loc_df[loc_df['Chr'].isin(valid_chr)]
        gene_pos = {}
        pos_gene = {}

        grouped = loc_df.groupby('Chr')
        
        for group_name, group_df in grouped:
            sorted_group = group_df.sort_values(by='Start')

            gene_pos_in = {b: i + 1 for i, b in enumerate(sorted_group.index)}
            pos_gene_in = {i + 1: b for i, b in enumerate(sorted_group.index)}

            # Agrega el diccionario interno al diccionario principal
            gene_pos[group_name] = gene_pos_in
            pos_gene[group_name] = pos_gene_in
        
        return gene_pos, pos_gene

 