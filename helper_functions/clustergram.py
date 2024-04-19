import pandas as pd
import numpy as np
import math
import re

import streamlit as st

import seaborn as sns
import matplotlib.pyplot as plt

import textwrap

from helper_functions.session_state import ss

class GeneHandler():
    def genes_used(self, degs, useDEG= None, textgene=None):
        if useDEG is not None:
            get_dict_values = [degs[s].index.to_list() for s in useDEG] # Get the deg list from the selected keys
            get_dict = {k:degs[k].index.to_list() for k in useDEG}
            get_dict = {k:v for k,v in get_dict.items() if len(v) !=0}
            flatten_dict_values = list(set([item for sublist in get_dict_values for item in sublist]))
            gene_final = flatten_dict_values
        
        if textgene is not None:
            genes = textgene.replace(";", ",").replace(" ", ",").replace("\n", ",").split(',')
            # where the user may use multiple delimiters, convert the other delimiters to commas, and then split by comma
            gene_final = [x.upper() for x in set(genes) if x != ""] # remove spaces in between and make sure to capitalise genes
            get_dict = {'user_genes':gene_final}
            get_dict = {k:v for k,v in get_dict.items() if len(v) !=0}
        
        if textgene == "None":
            gene_final = []
            get_dict = {}
        return gene_final, get_dict
    
    def get_gene_vals(self, log_dict_ready, genes_used=None):
        compiled_logFC = pd.DataFrame()
        for k,v in log_dict_ready.items():
            selected_gene_vals = v.loc[~v.index.duplicated(keep='first')]
            selected_gene_vals = selected_gene_vals[selected_gene_vals.index.isin(genes_used)]
            logFC_selected_gene_vals = selected_gene_vals.filter(regex="log2FC", axis=1)
            logFC_selected_gene_vals = logFC_selected_gene_vals.add_prefix(f"{k}_")
            compiled_logFC = pd.concat([compiled_logFC, logFC_selected_gene_vals], axis=1)
        return compiled_logFC

class Clustergram():

    @st.cache_data
    def cluster_plot(_self,
                     compiled_logFC,
                     gene_dict,
                     vminmax = (-2.0, 2.0),
                     cbar_left = 0.96,
                     cbar_bottom = 0.02,
                     cbar_width = 0.15,
                     cbar_height = 0.02,
                     width = 10, 
                     height = 10,
                     dendrogram_r = 0.2,
                     dendrogram_c = 0.12,
                     cluster_cols = True):
        colnames_whitespaced = [i.replace("_"," ") for i in compiled_logFC.columns]
        wrap_colnames = ["\n".join(textwrap.wrap(a, width=30, break_long_words=False)) for a in colnames_whitespaced]

        # drop those FCs with null values
        null_fc = compiled_logFC[compiled_logFC.apply(lambda x: pd.isna(x).any(), axis = 1)].index.to_list()
        reformatted_logFC = compiled_logFC[~compiled_logFC.index.isin(null_fc)]
        reformatted_logFC.columns = wrap_colnames

        dendrogram_c = 0.0 if not cluster_cols else dendrogram_c

        if reformatted_logFC.shape[0] > 3:
            g = sns.clustermap(reformatted_logFC,
                                cmap="vlag",
                                method='average',
                                cbar_pos=(cbar_left, cbar_bottom, cbar_width, cbar_height),
                                center=0, 
                                vmin = vminmax[0],
                                vmax = vminmax[1],
                                z_score=None,
                                col_cluster=cluster_cols,
                                yticklabels=True,
                                figsize=(width, height),
                                dendrogram_ratio=(dendrogram_r, dendrogram_c),
                                linewidths=1, linecolor='white',
                                cbar_kws = {"label": "log2FC", 'orientation':'horizontal', 'ticks':[vminmax[0], 0, vminmax[1]]})

            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=11, rotation=0)
            g.ax_heatmap.set_ylabel("")
            titles = '\n'.join([i for i in gene_dict.keys()])
            g.figure.suptitle(f"Clustergram from \n {titles}", x=0.5, y=1.04, fontsize=14, fontweight='bold')
            for _, spine in g.ax_heatmap.spines.items():
                spine.set_visible(True)
                spine.set_edgecolor("black")
            return g, null_fc
        
        else:
            return None, null_fc

genePP = GeneHandler()
clustergram = Clustergram()