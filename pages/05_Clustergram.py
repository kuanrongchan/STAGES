import pandas as pd
import numpy as np
import math
import re

# from scipy.spatial.distance import pdist, squareform
# from scipy.cluster.hierarchy import dendrogram, linkage

import streamlit as st

import plotly.colors as pc
import matplotlib.pyplot as plt
import seaborn as sns

from helper_functions.session_state import ss
from helper_functions.downloads import file_downloads
from helper_functions.clustergram import genePP, clustergram

ss.initialise_state({'cluster_useDEG': None,
                     'cluster_textgene':'IFIT3;IL15;DDX60;ILK;IGFLR1',
                     'clust_cols': True,
                     'clust_dendror':0.12,
                     'clust_dendroc':0.08,
                     'clust_width': 5,
                     'clust_height': 8,
                     'clust_vminmax':(-2.0, 2.0),
                     'clust_cbarleft':0.90,
                     'clust_cbarbottom':0.05,
                     'clust_cbarwidth':0.15,
                     'clust_cbarheight':0.02,
                     'clust_genelist':None,
                     'clust_genedict':None,
                     'clust_genevals':None,
                     'clust_submit':True,
                     'clustergram_plot':None
                     })

st.header("Clustergram")

try:
    clust_opts = st.sidebar.expander("Clustergram options", expanded=True)
    degs = st.session_state['degs']
    if degs is not None: # If there were DEGs already
        cluster_useDEG = clust_opts.multiselect("Select DEGs to use in clustergram",
                                                help="Leave this field blank if you wish to input a custom set of gene names",
                                                options = list(degs.keys()),
                                                default=st.session_state['cluster_useDEG'])


        if len(cluster_useDEG) == 0: # if no DEGs selected, provide text area for input
            cluster_textgene = clust_opts.text_area("Enter custom list of genes here (if not using DEGs)",
                                                    value = st.session_state['cluster_textgene'],
                                                    placeholder="Enter genes with one of the following delimiters: line breaks, commas, or semicolons")
            ss.save_state({'cluster_textgene':cluster_textgene,
                        'cluster_useDEG':None})
        else: # If DEGs selected, use the cluster_useDEG and revert textgene to None
            ss.save_state({'cluster_textgene': None,
                        'cluster_useDEG':cluster_useDEG})
        
    else:
        cluster_textgene = clust_opts.text_area("Enter custom list of genes here",
                                                value = st.session_state['cluster_textgene'],
                                                placeholder="Enter genes with one of the following delimiters: line breaks, commas, or semicolons")
        cluster_textgene = ss.save_state({'cluster_textgene':cluster_textgene})

    # Customisation settings
    ## Clustering options
    with clust_opts:
        cluster_cols = st.checkbox("Cluster columns", value = st.session_state['clust_cols'], on_change=ss.binaryswitch, args=('clust_cols', ))
        dendrogram_r = st.number_input("Adjust relative row dendrogram length", min_value=0.01, max_value=1.00, step=0.01, value= st.session_state['clust_dendror'])
        if cluster_cols:
            dendrogram_c = st.number_input("Adjust relative column dendrogram height", min_value=0.01, max_value=1.00, step=0.01, value=st.session_state['clust_dendroc'])
            ss.save_state({'clust_dendroc':dendrogram_c})

        ## Size and colourbar options
        c_width = st.number_input("Clustergram width (in inches)", min_value=3, max_value=30, step=1, value=st.session_state['clust_width'])
        c_height = st.number_input("Clustergram height (in inches)", min_value=3, max_value=60, step=1, value=st.session_state['clust_height'])
        c_vminmax = st.slider("Adjust minimum and maximum values of the colour bar", min_value=-10.0, max_value=10.0, step = 0.5, value=st.session_state['clust_vminmax'])
        cbar_left = st.number_input("Adjust colourbar position from left", min_value=0.0, max_value=1.0, step=0.01, value = st.session_state['clust_cbarleft'])
        cbar_bottom = st.number_input("Adjust colourbar position from bottom", min_value=0.0, max_value=1.0, step=0.01, value = st.session_state['clust_cbarbottom'])
        cbar_width = st.number_input("Adjust colourbar width", min_value=0.0, max_value=1.0, step=0.01, value = st.session_state['clust_cbarwidth'])
        cbar_height = st.number_input("Adjust colourbar height", min_value=0.0, max_value=1.0, step=0.01, value = st.session_state['clust_cbarheight'])

        plot_clust = st.checkbox("Plot clustergram", value=st.session_state['clust_submit'], on_change=ss.binaryswitch, args=('clust_submit', ))

    if plot_clust:
        ss.save_state({'clust_dendror':dendrogram_r,
                    'clust_width':c_width,
                    'clust_height':c_height,
                    'clust_vminmax':c_vminmax,
                    'clust_cbarleft':cbar_left,
                    'clust_cbarbottom':cbar_bottom,
                    'clust_cbarwidth':cbar_width,
                    'clust_cbarheight':cbar_height
                    })
        get_genes, gene_dict = genePP.genes_used(degs=degs, useDEG=st.session_state['cluster_useDEG'], textgene=st.session_state['cluster_textgene'])
        gene_vals = genePP.get_gene_vals(st.session_state['log_dict_ready'], get_genes)
        ss.save_state({'clust_genelist':get_genes,
                    'clust_genedict':gene_dict,
                    'clust_genevals':gene_vals})
        
        try:
            get_clustergram = clustergram.cluster_plot(st.session_state['clust_genevals'],
                                                    gene_dict=st.session_state['clust_genedict'],
                                                    vminmax=st.session_state['clust_vminmax'],
                                                    cbar_left=st.session_state['clust_cbarleft'],
                                                    cbar_bottom=st.session_state['clust_cbarbottom'],
                                                    cbar_width=st.session_state['clust_cbarwidth'],
                                                    cbar_height=st.session_state['clust_cbarheight'],
                                                    width=st.session_state['clust_width'],
                                                    height=st.session_state['clust_height'],
                                                    dendrogram_r=st.session_state['clust_dendror'],
                                                    dendrogram_c=st.session_state['clust_dendroc'],
                                                    cluster_cols=st.session_state['clust_cols']
                                                    )
            ss.save_state({'clustergram_plot':get_clustergram})
            st.pyplot(st.session_state['clustergram_plot'])
            file_downloads.create_pdf(st.session_state['clustergram_plot'], "Clustergram", "pyplot")
        
        except ValueError:
            st.error("At least 2 genes must be entered!")

except AttributeError:
    st.error("Perhaps you forgot to run through the pre-processing step?")