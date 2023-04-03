import gseapy as gp
import pandas as pd
import numpy as np
import re

from helper_functions.session_state import ss
from helper_functions.uploads import fileuploads
from helper_functions.downloads import file_downloads
from helper_functions.clustergram import genePP
from helper_functions.gseapy_functions_stages import prnk

import streamlit as st

ss.initialise_state({'prerank_selected_df_idx':0,
                     'prerank_selected_df':None,
                     'prerank_by':None,
                     'prerank_choose_col_idx':0,
                     'prerank_choose_col':None,
                     'prerank_pthresh':0.05,
                     'prerank_showX':10,
                     'prerank_ht':1000,
                     'plot_prerank':True,
                     'geneset_dict': {"Blood Transcriptomic Modules (BTM)": "accessory_files/BTM.gmt",
                                      "Reactome 2021": "accessory_files/Reactome.gmt",
                                      "Vaccinomics (In-house)": "accessory_files/Vaccinomics.gmt",
                                      "GO Biological Process 2021": "GO_Biological_Process_2021",
                                      "GO Molecular Function 2021": "GO_Molecular_Function_2021",
                                      "GO Cellular Component 2021": "GO_Cellular_Component_2021",
                                      "KEGG 2021 Human": "KEGG_2021_Human",
                                    #   "KEGG 2019 Mouse":"KEGG_2019_Mouse",
                                      "HumanCyc 2016": "HumanCyc_2016"
                                      },
                     'add_geneset_in':None,
                     'add_geneset':None,
                     'geneset_prerank':'Blood Transcriptomic Modules (BTM)',
                     'prerank_res_all':None,
                     'prerank_res_sig':None,
                     'prerank_plots':None}
                     )

st.header("GSEA Preranked Analysis")

# try:
prnk_opts = st.sidebar.expander("GSEA Preranked Analysis Options", expanded=True)

add_geneset_in = prnk_opts.file_uploader("Upload a gene set here (optional)", type="gmt", accept_multiple_files=True, help="Reupload not required if gene set has already been uploaded once", key="add_geneset_prnk")
ss.save_state({'add_geneset_in':add_geneset_in})

prnk_plots_t, prnk_data_t = st.tabs(["Bar plots", "Data"])

if st.session_state['add_geneset_in'] is not None:
    add_geneset = fileuploads.gmt_to_dict(st.session_state['add_geneset_in'])
    ss.save_state({'add_geneset':add_geneset})
    # ss.save_state({'geneset_dict':st.session_state['geneset_dict']|st.session_state['add_geneset']}) # | merges both dictionaries for python 3.10 only
    ss.save_state({'geneset_dict':{**st.session_state['geneset_dict'], **st.session_state['add_geneset']}})
else:
    ss.save_state({'geneset_dict':st.session_state['geneset_dict']})


# Selecting genesets (BTM or reactome) to plot from a list
geneset_opts = list(st.session_state['geneset_dict'].keys())
geneset_prnk = prnk_opts.radio(label='Select a geneset for prerank', options=geneset_opts, index = geneset_opts.index(st.session_state['geneset_prerank']))
ss.save_state({'geneset_prerank':geneset_prnk})

df_opts = list(st.session_state['log_dict_ready'].keys())
prerank_selected_df = prnk_opts.selectbox("Select dataframe to use in GSEA preranked analysis", options = df_opts, index = st.session_state['prerank_selected_df_idx'])
prerank_selected_df_idx = df_opts.index(prerank_selected_df)
ss.save_state({'prerank_selected_df_idx': prerank_selected_df_idx,
            'prerank_selected_df':prerank_selected_df})

prerank_by = prnk.format_cols(st.session_state['log_dict_ready'], st.session_state['comparisons'], selected_df=st.session_state['prerank_selected_df'])
ss.save_state({'prerank_by': prerank_by})

col_opts = list(st.session_state['prerank_by'].keys())
prerank_choose_col = prnk_opts.selectbox("Choose log2 fold-change columns for GSEA preranked analysis", options=col_opts, index=st.session_state['prerank_choose_col_idx'])
prerank_choose_col_idx = col_opts.index(prerank_choose_col)
ss.save_state({'prerank_choose_col_idx':prerank_choose_col_idx,
            'prerank_choose_col':prerank_choose_col})

# prerank_pthresh = prnk_opts.number_input("Choose adjusted p-value threshold of enriched pathways to plot", min_value = 0.00, max_value=1.00, step=0.01, value = st.session_state['prerank_pthresh'])
prerank_showX = prnk_opts.number_input("Display top n pathways from selected p-value cutoff", min_value=1, max_value=100, step=1,
                                    value=st.session_state['prerank_showX'],
                                    help="Show only the top n pathways from a filtered set of pathways")
prerank_ht = prnk_opts.number_input("Bar plot height (in px)", min_value=200, max_value=1600, step=50, value=st.session_state['prerank_ht'])
ss.save_state({
                # 'prerank_pthresh': round(prerank_pthresh,2),
                'prerank_showX':prerank_showX,
                'prerank_ht':prerank_ht})
plot_prerank = prnk_opts.checkbox("Run GSEA Prerank", value = st.session_state['plot_prerank'], on_change=ss.binaryswitch, args=("plot_prerank", ))

if plot_prerank:
    get_col = {st.session_state['prerank_choose_col']:st.session_state['prerank_by'][st.session_state['prerank_choose_col']]}
    get_geneset = st.session_state['geneset_dict'][st.session_state['geneset_prerank']]
    all_res, sig_res = prnk.execute_prerank(col_dict=get_col,
                                            select_dataset=get_geneset,
                                            prerank_pthresh=0.05,
                                            prerank_showX=st.session_state['prerank_showX'])
    ss.save_state({'prerank_res_all':all_res,
                'prerank_res_sig':sig_res})
    sig_plots = prnk.prerank_barplot(st.session_state['prerank_res_sig'],
                                    selected_col = st.session_state['prerank_choose_col'],
                                    select_dataset=st.session_state['geneset_prerank'],
                                    # prerank_pthresh=st.session_state['prerank_pthresh'],
                                    prerank_showX=st.session_state['prerank_showX'],
                                    prerank_ht=st.session_state['prerank_ht'])
    ss.save_state({'prerank_plots':sig_plots})
    prnk_plots_t.plotly_chart(st.session_state['prerank_plots'], theme=None, use_container_width=False)
    with prnk_plots_t:
        file_downloads.create_pdf(st.session_state['prerank_plots'], "Prerank_plots", "plotly")

    with prnk_data_t:
        for k,v in sig_res.items():
            st.write(f"**{k}**")
            st.dataframe(v)
        st.download_button(label="Download GSEA Preranked Results",
                        data=file_downloads.to_excel(st.session_state['prerank_res_all'].values(),sheetnames=st.session_state['prerank_res_all'].keys()),
                        file_name="GSEAPreranked_results.xlsx")