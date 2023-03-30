import gseapy as gp
import pandas as pd
import numpy as np
import re


from helper_functions.session_state import ss
from helper_functions.uploads import fileuploads
from helper_functions.downloads import file_downloads
from helper_functions.clustergram import genePP
from helper_functions.gseapy_functions_stages import enr

import streamlit as st

ss.initialise_state({'add_geneset_in':None,
                     'add_geneset':None,
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
                      'geneset_enr':"Blood Transcriptomic Modules (BTM)",
                      'enr_useDEG':None,
                      'enr_textgene':'COL1A2;DCN;IL6;IL8;LIF;MGP;MMP1;MMP2;MMP9',
                      'enr_pthresh':0.05,
                      'enr_showX':10,
                      'plot_enr':True,
                      'enr_genedict':None,
                      'enr_ht':500,
                      'enr_res_all':None,
                      'enr_res_sig':None,
                      'enrichr_plots':None
                      }
                      )

st.header("Enrichr")
enr_opts = st.sidebar.expander("Enrichr Options", expanded=True)

add_geneset_in = enr_opts.file_uploader("Upload a gene set here (optional)", type="gmt", accept_multiple_files=True, help="Reupload not required if gene set has already been uploaded once")
ss.save_state({'add_geneset_in':add_geneset_in})

enr_plots_t, enr_data_t = st.tabs(["Bar plots", "Data"])

if st.session_state['add_geneset_in'] is not None:
    add_geneset = fileuploads.gmt_to_dict(st.session_state['add_geneset_in'])
    ss.save_state({'add_geneset':add_geneset})
    # ss.save_state({'geneset_dict':st.session_state['geneset_dict']|st.session_state['add_geneset']}) # | merges both dictionaries for python 3.10 only
    ss.save_state({'geneset_dict':{**st.session_state['geneset_dict'], **st.session_state['add_geneset']}})

# Selecting genesets (BTM or reactome) to plot from a list
geneset_opts = list(st.session_state['geneset_dict'].keys())
geneset = enr_opts.radio(label='Select a geneset', options=geneset_opts, index = geneset_opts.index(st.session_state['geneset_enr']))
ss.save_state({'geneset_enr':geneset})

degs = st.session_state['degs'] if "degs" in st.session_state else None
if degs is not None: # If there were DEGs already
    enr_useDEG = enr_opts.multiselect("Select DEGs to use in Enrichr",
                                      help="Leave this field blank if you wish to input a custom set of gene names",
                                      options = list(degs.keys()),
                                      default=st.session_state['enr_useDEG'])
    ss.save_state({'enr_useDEG':enr_useDEG})

    if len(enr_useDEG) == 0: # if no DEGs selected, provide text area for input
        enr_textgene = enr_opts.text_area("Enter custom list of genes here (if not using DEGs)",
                                          value = st.session_state['enr_textgene'],
                                          placeholder="Enter genes with one of the following delimiters: line breaks, commas, or semicolons")
        ss.save_state({'enr_textgene':enr_textgene,
                       'enr_useDEG':None})
    else: # If DEGs selected, use the enr_useDEG and revert textgene to None
        ss.save_state({'enr_textgene': None})
    
else:
    enr_textgene = enr_opts.text_area("Enter custom list of genes here",
                                            value = st.session_state['enr_textgene'],
                                            placeholder="Enter genes with one of the following delimiters: line breaks, commas, or semicolons")
    ss.save_state({'enr_textgene':enr_textgene})


# enr_pthresh = enr_opts.number_input("Choose adjusted p-value threshold of enriched pathways to plot", min_value = 0.00, max_value=1.00, step=0.01, value = st.session_state['enr_pthresh'])
enr_showX = enr_opts.number_input("Display top n pathways from selected p-value cutoff", min_value=1, max_value=100, step=1,
                                     value=st.session_state['enr_showX'],
                                     help="Show only the top n pathways from a filtered set of pathways")
enr_ht = enr_opts.number_input("Bar plot height (in px)", min_value=200, max_value=1600, step=50, value=st.session_state['enr_ht'])
ss.save_state({
            #    'enr_pthresh': round(enr_pthresh,2),
               'enr_showX':enr_showX,
               'enr_ht':enr_ht})
plot_enr = enr_opts.checkbox("Run Enrichr", value = st.session_state['plot_enr'], on_change=ss.binaryswitch, args=("plot_enr", ))

if plot_enr:
    _, gene_dict = genePP.genes_used(degs=degs, useDEG=st.session_state['enr_useDEG'], textgene=st.session_state['enr_textgene'])
    ss.save_state({'enr_genedict':gene_dict})
    get_geneset = st.session_state['geneset_dict'][st.session_state['geneset_enr']]
    if len(gene_dict) == 0:
        st.warning("Please ensure that there is more than 1 gene from DEGs or genes are manually entered!")
    else:
        res_all, res_sig = enr.execute_enrichr(gene_dict=st.session_state['enr_genedict'], select_dataset=get_geneset, enr_pthresh=st.session_state['enr_pthresh'], enr_showX=st.session_state['enr_showX'])
        ss.save_state({'enr_res_all':res_all,
                       'enr_res_sig':res_sig})
        
        if 'bar_fc' in st.session_state:
            enr_plots = enr.enr_barplot(st.session_state['enr_res_sig'],
                                        enr_useDEG=st.session_state['enr_useDEG'],
                                        enr_pthresh=0.05,
                                        deg_fc=st.session_state['bar_fc'],
                                        deg_pval=st.session_state['bar_pval'],
                                        use_corrected_pval=st.session_state['use_corrected_pval'],
                                        select_dataset=st.session_state['geneset_enr'],
                                        # enr_pthresh=st.session_state['enr_pthresh'],
                                        enr_showX=st.session_state['enr_showX'],
                                        enr_ht=st.session_state['enr_ht'])
        else:
            enr_plots = enr.enr_barplot(st.session_state['enr_res_sig'],
                            enr_useDEG=st.session_state['enr_useDEG'],
                            enr_pthresh=0.05,
                            deg_fc=None,
                            deg_pval=None,
                            use_corrected_pval=st.session_state['use_corrected_pval'],
                            select_dataset=st.session_state['geneset_enr'],
                            # enr_pthresh=st.session_state['enr_pthresh'],
                            enr_showX=st.session_state['enr_showX'],
                            enr_ht=st.session_state['enr_ht'])
        ss.save_state({'enrichr_plots':enr_plots})
        with enr_plots_t:
            st.plotly_chart(enr_plots, theme=None, use_container_width=False)
            file_downloads.create_pdf(enr_plots, fn="Enrichr plots", graph_module="plotly")
        
        with enr_data_t:
            for k,v in res_all.items():
                st.write(f"**{k}**")
                st.dataframe(v)
            st.download_button(label="Download Enrichr Results",
                               data=file_downloads.to_excel(st.session_state['enr_res_all'].values(), sheetnames=st.session_state['enr_res_all'].keys()),
                               file_name="Enrichr_results.xlsx")
