import gseapy as gp
import pandas as pd
import numpy as np
import re

from helper_functions.session_state import ss
from helper_functions.uploads import fileuploads
from helper_functions.downloads import file_downloads
from helper_functions.clustergram import genePP
from helper_functions.string import stages_str

import streamlit as st

ss.initialise_state({'string_useDEG':None,
                     'string_textgene':'COL1A2;DCN;IL6;IL8;LIF;MGP;MMP1;MMP2;MMP9',
                     'plot_string':True,
                     'string_plots':None})

st.header("STRING Network Query")
str_netplots, str_other = st.tabs(['STRING Network Plots', "Data"])
string_opts = st.sidebar.expander("STRINGdb Query Options", expanded=True)
degs = st.session_state['degs'] if "degs" in st.session_state else None
if degs is not None: # If there were DEGs already
    string_useDEG = string_opts.multiselect("Select DEGs to use in Enrichr",
                                      help="Leave this field blank if you wish to input a custom set of gene names",
                                      options = list(degs.keys()),
                                      default=st.session_state['string_useDEG'])


    if len(string_useDEG) == 0: # if no DEGs selected, provide text area for input
        string_textgene = string_opts.text_area("Enter custom list of genes here (if not using DEGs)",
                                          value = st.session_state['string_textgene'],
                                          placeholder="Enter genes with one of the following delimiters: line breaks, commas, or semicolons")
        ss.save_state({'string_textgene':string_textgene,
                       'string_useDEG':None})
    else: # If DEGs selected, use the string_useDEG and revert textgene to None
        ss.save_state({'string_textgene': None,
                       'string_useDEG':string_useDEG})
    
else:
    string_textgene = string_opts.text_area("Enter custom list of genes here",
                                            value = st.session_state['string_textgene'],
                                            placeholder="Enter genes with one of the following delimiters: line breaks, commas, or semicolons")
    ss.save_state({'string_textgene':string_textgene})

plot_string = string_opts.checkbox("Plot STRING Network", value=st.session_state['plot_string'], on_change=ss.binaryswitch, args=("plot_string", ))

if plot_string:
    _, gene_dict = genePP.genes_used(degs=degs, useDEG=st.session_state['string_useDEG'], textgene=st.session_state['string_textgene'])

    if len(gene_dict) != 0:
        str_network, tozip = stages_str.string_query(gene_dict)
        ss.save_state({'string_plots':tozip})
        with str_netplots:
            for k, v in str_network.items():
                st.write(f"**{k}**")
                st.image(v, use_column_width=True)
            file_downloads.zip_imgs(tozip)
    else:
        st.warning("Please ensure that there is more than 1 gene from DEGs or genes are manually entered!")