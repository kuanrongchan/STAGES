import pandas as pd
import numpy as np
import regex as re

# Stats modules
from statsmodels.stats import multitest
from anndata import AnnData
import decoupler as dc
import scanpy as sc

import streamlit as st
from streamlit_tags import st_tags, st_tags_sidebar

import plotly.colors as pc
import matplotlib.pyplot as plt 

from helper_functions.preprocessing import tested, counts_pp
from helper_functions.session_state import ss
from helper_functions.downloads import file_downloads

st.session_state.update(st.session_state)

ss.initialise_state({'test_fdr':'None',
                     'vthresh':0,
                     'adata':None,
                     'violin1':None,
                     'violin2':None,
                     'comparisons':None,
                     'comp_var': 0,
                     'comparison_options': [],
                     'baseline': None,
                     'comparisonopts_nobaseline':[],
                     'against_baseline':None,
                     'equalvar':True,
                     'use_corrected_pval':False,
                     'submit_comparison':False,
                     'ready':None,
                     'log_dict_ready':None})

try:
    exprdict, metadatadict, anovadict = st.session_state['expr_dict'], st.session_state['meta_dict'], st.session_state['anova_dict']

    adjusted_dfs = {}
    prep_exp = st.sidebar.expander("Pre-processing Options")
    padj_mtds = {"None":None, "Bonferroni":"bonferroni", "Sidak":'sidak', 'Holm-Sidak':'holm-sidak', 'Holm':'holm', 'Simes-Hochberg':'simes-hochberg', 'Hommel':'hommel',
                'Benjamini-Hochberg FDR':'fdr_bh', 'Benjamini-Yekutieli FDR':'fdr_by',
                'Two-stage Benjamini-Hochberg FDR':'fdr_tsbh', 'Two-stage Benjamini-Yekutieli FDR':'fdr_tsbky'}

    # Conditions here should mainly be
    ## 1. Process ANOVA data (not None)
    ## 2. Process Counts + metadata
    ## NOTE: there should not be both anova AND count data

    if anovadict is not None:
        comps = tested.comparison_finder(anovadict)
        test_fdr = prep_exp.selectbox("Select multiple test correction method", options = list(padj_mtds.keys()), index = list(padj_mtds.keys()).index(st.session_state['test_fdr']))
        ss.save_state({'test_fdr':test_fdr})
        test_fdr_match = padj_mtds[test_fdr]

        if st.session_state['test_fdr'] != 0:
                use_corrected_pval = prep_exp.checkbox("Use corrected p-values for subsequent analyses", value=st.session_state['use_corrected_pval'], on_change=ss.binaryswitch, args= ('use_corrected_pval', ))
        else:
            ss.save_state({'use_corrected_pval':False})
        
        if test_fdr_match is not None:
            for k,v in comps.items():
                anova_file = anovadict[k]
                adj_df_per_k = pd.DataFrame()
                for comp in v:
                    comp_df = anova_file.filter(regex = comp, axis=1)
                    pval_col = comp_df.filter(regex = "^pval", axis=1) # need to accommodate other regexes????
                    pval_col = pval_col.sort_values(by=pval_col.columns[0], ascending = True)
                    pval_array = pval_col.iloc[:,0].to_numpy()
                    rej, corrected, alphacSidak, alphacBonf = multitest.multipletests(pval_array,
                                                                                    method=test_fdr_match,
                                                                                    is_sorted = True)
                    corrected_vals = pd.DataFrame(data = {f"adj_pval_{comp}":corrected}, index = pval_col.index)
                    try:
                        adj_df_per_k = pd.concat([adj_df_per_k, comp_df, corrected_vals], axis=1, verify_integrity=True)
                    except ValueError:
                        st.error("Duplicated columns found. Perhaps you already have adjusted p-values? If so, opt for None in multiple test correction and tick the checkbox to use corrected p-values for subsequent analysis.")
                        st.stop()

                adjusted_dfs[k] = adj_df_per_k
            ss.save_state({'test_fdr':test_fdr, 'comparisons':comps, 'ready':adjusted_dfs})
        
        else:
            adjusted_dfs = anovadict
            ss.save_state({'test_fdr':"None", 'comparisons':comps, 'ready': adjusted_dfs})

    elif exprdict is not None and metadatadict is not None: # RNAseq or microarray data
        exprdict = {k:v.groupby(v.index).mean() for k,v in exprdict.items()}
        expr_obj = exprdict[list(exprdict.keys())[0]].T.sort_index(axis=0, ascending=True)
        expr_key = list(exprdict.keys())[0]
        meta_obj = metadatadict[list(metadatadict.keys())[0]].sort_index(axis=0, ascending=True)

        adata = AnnData(expr_obj, obs = meta_obj, dtype='float32') # expr data should be genes in cols, subjects in rows, obs should be the same
        ss.save_state({"adata":adata})

        if st.session_state['file_type'] == "RNAseq Counts": # specifically RNAseq data
            st.header("Count Normalisation")
            split_long_violins = counts_pp.chunks(list(adata.obs_names), chunk_size=12)
            adata = st.session_state['adata']
            # Provide info here
            st.info("The violin plots From the pre-processing options in the sidebar, select a threshold value that separates the violin plots into two (narrowest point).")
            bef, aft = st.tabs(['Before pre-processing', 'After pre-processing'])
            
            # Create violin plot
            max_y = counts_pp.violin_maxy(adata)
            violin_thresh = prep_exp.slider(label = "Select violin plot threshold", min_value = 0, max_value= max_y,
                                    help = "This function will retain genes whose counts are above the threshold.",
                                    value=st.session_state['vthresh'])
            ss.save_state({'vthresh':violin_thresh})

            violin1, vaxes = st.cache_data(counts_pp.multiviolin)(adata, split_long_violins=split_long_violins)
            if len(split_long_violins) == 1:
                _ = vaxes.axhline(y=st.session_state['vthresh'], color ='red', linewidth=1, linestyle='--')
            else:
                _ = [ax.axhline(y = st.session_state['vthresh'], color = 'red', linewidth = 1, linestyle = '--') for ax in vaxes]
            ss.save_state({'violin1':violin1})
            bef.pyplot(st.session_state['violin1'])

            adata = st.session_state['adata']
            adata_vars = list(adata.obs.columns)
            comp_var = prep_exp.selectbox(label="Select variable to use for comparison", options = adata_vars, index= st.session_state['comp_var'])
            comparison_options = adata.obs[comp_var].unique()
            ss.save_state({'comp_var': adata_vars.index(comp_var), 'comparison_options':comparison_options})
            baseline = prep_exp.multiselect(label=f"Select groups within **{adata_vars[st.session_state['comp_var']]}** as the baseline for comparison ie. choose A where B vs A", options = st.session_state['comparison_options'], default = st.session_state['baseline'])
            ss.save_state({'baseline':baseline})
            comparisonopts_nobaseline = [i for i in comparison_options if i not in st.session_state['baseline']]
            ss.save_state({'comparisonopts_nobaseline':comparisonopts_nobaseline})
            against_baseline = prep_exp.multiselect(label=f"Select groups within **{adata_vars[st.session_state['comp_var']]}** to compare against baseline ie choose B where B vs A", options = st.session_state['comparisonopts_nobaseline'], default = st.session_state['against_baseline'])
            ss.save_state({'against_baseline':against_baseline})
            equalvar = prep_exp.checkbox(label="Assume equal population variance", value = st.session_state['equalvar'], on_change=ss.binaryswitch, args = ('equalvar', ))
            test_fdr = prep_exp.selectbox("Select multiple test correction method", options = list(padj_mtds.keys()), index = list(padj_mtds.keys()).index(st.session_state['test_fdr']))
            ss.save_state({'test_fdr':test_fdr})
            if st.session_state['test_fdr'] != 0:
                use_corrected_pval = prep_exp.checkbox("Use corrected p-values for subsequent analyses", value=st.session_state['use_corrected_pval'], on_change=ss.binaryswitch, args= ('use_corrected_pval', ))
            else:
                ss.save_state({'use_corrected_pval':False})
            submit_comparison = prep_exp.checkbox("Selection complete", value = st.session_state['submit_comparison'], on_change=ss.binaryswitch, args=('submit_comparison', ))

            if submit_comparison:
                dc.mask_features(adata, log=True, thr=st.session_state['vthresh'])
                sc.pp.filter_genes(adata, min_cells=adata.shape[0])
                violin2, vaxes2 = counts_pp.multiviolin(adata, split_long_violins=split_long_violins)
                ss.save_state({'violin2':violin2, 'adata':adata})
                aft.pyplot(st.session_state['violin2'])

                ratios = counts_pp.ratio(adata, comp_var = adata_vars[st.session_state['comp_var']],
                                        baseline= st.session_state['baseline'],
                                        against_baseline= st.session_state['against_baseline'],
                                        is_log=False)
                ttest = counts_pp.pval_scipy(adata, comp_var = adata_vars[st.session_state['comp_var']],
                                            baseline= st.session_state['baseline'],
                                            against_baseline= st.session_state['against_baseline'],
                                            equalvar= st.session_state['equalvar'])
                
                comps = tested.comparison_finder({expr_key:ttest})
                test_fdr_match = padj_mtds[st.session_state['test_fdr']]
                adj_df_per_k = pd.DataFrame()
                if test_fdr_match is not None:
                    corrected_pval_df = pd.DataFrame()
                    for c in ttest.columns:
                        pval_col = ttest.loc[:,c].sort_values(ascending = True)
                        pval_array = pval_col.to_numpy()
                        rej, corrected, alphacSidak, alphacBonf = multitest.multipletests(pval_array,
                                                                                        method=test_fdr_match,
                                                                                        is_sorted = True)
                        corrected_vals = pd.DataFrame(data = {f"adj_{c}":corrected}, index = pval_col.index)
                        corrected_pval_df = pd.concat([corrected_pval_df, corrected_vals], axis=1)
                    adj_df_per_k = pd.concat([adj_df_per_k, ratios, ttest, corrected_pval_df], axis=1)
                else:
                    adj_df_per_k = pd.concat([ratios, ttest], axis=1)

                sort_by_comparison = pd.concat([adj_df_per_k.filter(regex=comp, axis=1) for comp in comps[expr_key]], axis=1)
                ss.save_state({'ready': {expr_key:sort_by_comparison}, 'comparisons':comps})

        else:
            adata = st.session_state['adata']
            adata_vars = list(adata.obs.columns)
            comp_var = prep_exp.selectbox(label="Select variable to use for comparison", options = adata_vars, index= st.session_state['comp_var'])
            comparison_options = adata.obs[comp_var].unique()
            ss.save_state({'comp_var': adata_vars.index(comp_var), 'comparison_options':comparison_options})
            baseline = prep_exp.multiselect(label=f"Select the groups within **{adata_vars[st.session_state['comp_var']]}** as the baseline for comparison ie. choose A where B vs A", options = st.session_state['comparison_options'], default = st.session_state['baseline'])
            ss.save_state({'baseline':baseline})
            comparisonopts_nobaseline = [i for i in comparison_options if i not in st.session_state['baseline']]
            ss.save_state({'comparisonopts_nobaseline':comparisonopts_nobaseline})
            against_baseline = prep_exp.multiselect(label=f"Select the groups within **{adata_vars[st.session_state['comp_var']]}** to compare against baseline ie choose B where B vs A", options = st.session_state['comparisonopts_nobaseline'], default = st.session_state['against_baseline'])
            ss.save_state({'against_baseline':against_baseline})
            equalvar = prep_exp.checkbox(label="Assume equal population variance", value = st.session_state['equalvar'], on_change=ss.binaryswitch, args = ('equalvar', ))
            test_fdr = prep_exp.selectbox("Select multiple test correction method", options = list(padj_mtds.keys()), index = list(padj_mtds.keys()).index(st.session_state['test_fdr']))
            ss.save_state({'test_fdr':test_fdr})
            if st.session_state['test_fdr'] != 0:
                use_corrected_pval = prep_exp.checkbox("Use corrected p-values for subsequent analyses", value=st.session_state['use_corrected_pval'], on_change=ss.binaryswitch, args= ('use_corrected_pval', ))
            else:
                ss.save_state({'use_corrected_pval':False})
            submit_comparison = prep_exp.checkbox("Selection complete", on_change=ss.binaryswitch, args=('submit_comparison', ))
            
            if submit_comparison:
                ratios = counts_pp.ratio(adata, comp_var = adata_vars[st.session_state['comp_var']],
                                        baseline= st.session_state['baseline'],
                                        against_baseline= st.session_state['against_baseline'],
                                        is_log = True)
                ttest = counts_pp.pval_scipy(adata, comp_var = adata_vars[st.session_state['comp_var']],
                                            baseline= st.session_state['baseline'],
                                            against_baseline= st.session_state['against_baseline'],
                                            equalvar= st.session_state['equalvar'])
                
                comps = tested.comparison_finder({expr_key:ttest})
                test_fdr_match = padj_mtds[st.session_state['test_fdr']]
                adj_df_per_k = pd.DataFrame()
                if test_fdr_match is not None:
                    corrected_pval_df = pd.DataFrame()
                    for c in ttest.columns:
                        pval_col = ttest.loc[:,c].sort_values(ascending = True)
                        pval_array = pval_col.to_numpy()
                        rej, corrected, alphacSidak, alphacBonf = multitest.multipletests(pval_array,
                                                                                        method=test_fdr_match,
                                                                                        is_sorted = True)
                        corrected_vals = pd.DataFrame(data = {f"adj_{c}":corrected}, index = pval_col.index)
                        corrected_pval_df = pd.concat([corrected_pval_df, corrected_vals], axis=1)
                    adj_df_per_k = pd.concat([adj_df_per_k, ratios, ttest, corrected_pval_df], axis=1)
                else:
                    adj_df_per_k = pd.concat([ratios, ttest], axis=1)

                sort_by_comparison = pd.concat([adj_df_per_k.filter(regex=comp, axis=1) for comp in comps[expr_key]], axis=1)
                ss.save_state({'ready': {expr_key:sort_by_comparison}, 'comparisons':comps})

    if st.session_state['ready'] is not None:
        log_dict = tested.log_transform(st.session_state['ready'], comparison_dict=st.session_state['comparisons'], use_corrected_pval=st.session_state['use_corrected_pval'])
        ss.save_state({'log_dict_ready':log_dict})
        st.header("Pre-processed data (ratios and p-values)")
        for k,v in st.session_state['ready'].items():
            st.subheader(k)
            st.dataframe(v)
        st.download_button(label="Download Processed Data", data=file_downloads.to_excel(st.session_state['ready'].values(),
                                                                                         sheetnames=st.session_state['ready'].keys()), file_name="processed_data.xlsx")
        

except KeyError:
    st.error("Perhaps you forgot to upload a dataset or use the demo data?")