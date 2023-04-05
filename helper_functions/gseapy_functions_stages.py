import plotly.graph_objects as go
import plotly.figure_factory as ff
from plotly.subplots import make_subplots
import textwrap

import numpy as np
import pandas as pd
import gseapy as gp
from gseapy.biomart import Biomart

import streamlit as st

from datetime import datetime, timedelta
import pytz

class Enrichr_STAGES():
    '''
    Class to run enrichr functions for STAGES
    '''

    def background_enrichr(_self, filename):
        start_date = datetime(2023, 3, 23, tzinfo=pytz.timezone("Asia/Singapore"))
        today = datetime.now(tz=pytz.timezone("Asia/Singapore"))
        interval = timedelta(days = 90)
        if today >= start_date + interval:
            bm = Biomart()
            query = bm.query(dataset='hsapiens_gene_ensembl', #  e.g. 'hsapiens_gene_ensembl'
                             attributes=['ensembl_gene_id', 'external_gene_name', 'entrezgene_id'],
                             filename="accessory_files/hsapiens_gene_ensembl.txt")
            query = query.dropna(subset=["entrezgene_id"], inplace=True)
            start_date = today
        else:
            query = pd.read_csv(filename, sep="\t",index_col=0)
            query = query.dropna(subset=["entrezgene_id"], inplace=True)
        return query
    
    @st.cache_data
    def execute_enrichr(_self, gene_dict, select_dataset, enr_pthresh=0.05, enr_showall=True, enr_showX=10):
        '''
        Parameters
        ----------
        gene_dict: dict | keys containing deg keys or user_genes and values containing list of genes to use
        select_dataset: str | one of the gene sets for enrichr
        enr_pthresh: float | pvalue to filter pathways by
        enr_showX: int | number of pathways to display from filtered dataset
        '''
        enr_significant, enr_all = {}, {}
        non_zero = {k:gene_dict[k] for k in gene_dict.keys() if len(gene_dict[k]) !=0} # in case deg sets with 0 genes were selected
        for k,v in non_zero.items():
            enr = gp.enrichr(
                gene_list=v,
                gene_sets=select_dataset,
                outdir=None,
                no_plot=True,
                cutoff=0.5,
                background=_self.background_enrichr(filename="accessory_files/hsapiens_gene_ensembl.txt")
            )

            # Sort values by adjusted p-value
            data = enr.results
            data.set_index("Term", inplace=True)
            data = data.sort_values(by=['Adjusted P-value'])
            # Drop the unimportant variables
            data_truncated = data.drop("Gene_set", axis=1)
            # Filter data by adjusted p-value and rank
            data_sig = data_truncated[(data_truncated['Adjusted P-value'] < enr_pthresh)]
            # Calculate -logP
            data_sig['-logadjP'] = np.log10(data_sig['Adjusted P-value']) * (-1)

            # Sort and return
            if not enr_showall:
                enr_significant[k] = data_sig.sort_values(by = "-logadjP", ascending = True).tail(enr_showX)
            else:
                enr_significant[k] = data_sig.sort_values(by = "-logadjP", ascending = True)
            enr_all[k] = data
        return enr_all, enr_significant
    
    @st.cache_data
    def enr_barplot(_self, enr_significant, enr_useDEG=None, deg_fc=1.30, deg_pval=0.05, use_corrected_pval=True, select_dataset="BTM", enr_pthresh=0.05, enr_showall=True, enr_showX=10, enr_ht=500):
        title_fmt = f"All significant Enrichr {select_dataset} pathways" if enr_showall else f"Top {enr_showX} Enrichr {select_dataset} pathways"
        if enr_useDEG is not None: # which implies it will require DEGs
            fig = make_subplots(rows=len(enr_significant), cols=1, subplot_titles=list(enr_significant.keys()),
                                x_title="-log10 (adjusted p-value)", shared_xaxes=True,
                                vertical_spacing=0.05) # make a subplot regardless
                
            i = 1
            for k,v in enr_significant.items(): # need to rethink this part as it creates a subplot for every comparison
                wrap_enr_sig = ["<br>".join(textwrap.wrap(a, 30)) for a in v.index]
                wrap_genes = ["<br>".join(textwrap.wrap(a, 80)) for a in v.Genes]
                marker_clr = "#EF553B" if "UP" in k else "#636EFA"
                fig.add_trace(go.Bar(x=v['-logadjP'], y=wrap_enr_sig,
                                     customdata=wrap_genes,
                                     orientation='h', marker_color=marker_clr,
                                     name="",
                                     hovertemplate="<b>%{y}</b><br>-logadjP: %{x}<br>Enriched genes: %{customdata}"),
                                     row = i, col = 1)
                i += 1
            fig.update_yaxes(title="Term", tickmode='linear', tick0=0, dtick=0, automargin=True)
            fig.update_layout(title=f"{title_fmt}<br>DEGs based on |FC| > {deg_fc}, {'adjusted p-value' if use_corrected_pval else 'p-value'} < {deg_pval}", title_x=0.5,
                              showlegend=False,
                              yaxis={'tickmode': 'linear'},
                              font=dict(family='Arial', size=14),
                              width = 750, height = enr_ht)
            

        else:
            fig = go.Figure()
            user_df = enr_significant['user_genes']
            wrap_usergenes = ["<br>".join(textwrap.wrap(a, 50)) for a in user_df.index]
            fig.add_trace(go.Bar(x=user_df['-logadjP'], y=wrap_usergenes,
                                 customdata=user_df['Genes'],
                                 orientation='h',
                                 marker_color="#4FC04F",
                                 name="",
                                 hovertemplate="<b>%{y}</b><br>-logadjP: %{x}<br>Enriched genes: %{customdata}")
                                 )
            fig.update_yaxes(title="Term", tickmode='linear', tick0=0, dtick=0, automargin=True)
            fig.update_layout(title=f"{title_fmt} for user-input genes", title_x=0.5,
                            showlegend=False,
                            xaxis_title = "-log10 (adjusted p-value)",
                            yaxis={'tickmode': 'linear'},
                            font=dict(family='Arial', size=14),
                            width = 750, height = enr_ht)
        return fig

class Prerank_STAGES():
    '''
    Class to run GSEA prerank functions for STAGES
    '''
    def format_cols(self, log_dict_ready, comparisons, selected_df):
        '''
        Parameters
        ----------
        log_dict_ready: dict | log2FC, -log10(p/adjP), p/adjP
        comparisons: dict | keys containing file name, values containing list of comparisons
        '''
        col_storage = {}
        comps = comparisons[selected_df]
        to_parse = log_dict_ready[selected_df]
        for comp in comps:
            logfc_comp = to_parse.loc[:,f"log2FC_{comp}"]
            logfc_comp = logfc_comp.reset_index()
            logfc_comp.columns = [0,1]
            logfc_comp = logfc_comp.sort_values(by=1, ascending=False)
            col_storage[f"{selected_df}_{comp}"] = logfc_comp
        return col_storage

    @st.cache_data
    def execute_prerank(_self, col_dict, select_dataset, prerank_pthresh=0.05, prerank_showX=10):
        prerank_results_dict = {}
        prerank_all_out = {}
        prerank_sig_out = {}
        for key, data in col_dict.items():
            running = gp.prerank(rnk=data,
                                gene_sets=select_dataset,
                                permutation_num=200,  # reduce number to speed up testing
                                outdir=None,
                                seed=123,
                                no_plot=True)
            prerank_results_dict[key] = running

        for key, result in prerank_results_dict.items():
            results = result.res2d
            prerank_all_out[key] = results

            ranked = results.loc[:,['Term', 'Lead_genes', 'NES', 'FDR q-val']]

            ranked = ranked.set_index("Term")
            ranked.index = [i.replace('"', "") for i in ranked.index]

            pos_nes = ranked[(ranked["NES"] > 0) & (ranked["FDR q-val"] < prerank_pthresh)]
            neg_nes = ranked[(ranked["NES"] < 0) & (ranked["FDR q-val"] < prerank_pthresh)]
            neg_nes["negative NES"] = neg_nes["NES"] * -1

            pos_nes_sort = pos_nes.sort_values(by=['NES'], ascending=True).tail(prerank_showX)
            pos_nes_sort.reset_index(inplace=True) # Compatibility issues with python 3.7, where the names argument was not valid at pandas 1.3.5
            pos_nes_sort = pos_nes_sort.rename(columns = {'index':'Term'}) # Subsequently due to this version error, have to manually rename the index column to term
            pos_nes_sort['direction'] = "positive"

            neg_nes_sort = neg_nes.sort_values(by=['negative NES'], ascending=True).tail(prerank_showX)
            neg_nes_sort.reset_index(inplace=True) # Compatibility issues with python 3.7, where the names argument was not valid at pandas 1.3.5
            neg_nes_sort = neg_nes_sort.rename(columns = {'index':'Term'}) # Subsequently due to this version error, have to manually rename the index column to term
            neg_nes_sort['direction'] = "negative"

            prerank_sig_out[f'Positive_enrichment_{key}'] = pos_nes_sort
            prerank_sig_out[f'Negative_enrichment_{key}'] = neg_nes_sort
        return prerank_all_out, prerank_sig_out

    @st.cache_data
    def prerank_barplot(_self, prerank_sig, selected_col, prerank_pthresh=0.05, select_dataset=None, prerank_showX=10, prerank_ht = 1000):
        fig = make_subplots(rows=2, cols=1, shared_xaxes=True, x_title = "|NES|", vertical_spacing=0.05, subplot_titles=list(prerank_sig.keys()))
        i = 1
        for k,v in prerank_sig.items():
            wrap_prnk_sig = ["<br>".join(textwrap.wrap(a, 30)) for a in v.Term]
            wrap_leadgenes = ["<br>".join(textwrap.wrap(a, 80)) for a in v.Lead_genes]
            if "Positive" in k:
                pos = fig.add_trace(go.Bar(x=v['NES'], y=wrap_prnk_sig,
                                        customdata=np.stack((v['Term'], v['FDR q-val'], wrap_leadgenes), axis=-1),
                                        orientation='h', marker_color="#EF553B",
                                        name="",
                                        hovertemplate="<b>%{customdata[0]}</b><br>" +
                                                      "|NES|: %{x}<br>"+ 
                                                      "FDR q-value: %{customdata[1]:.3f}<br>"+ 
                                                      "Leading edge genes: %{customdata[2]}",
                                        hoverlabel={'align':'left'}),
                                        row = i, col = 1)
                i += 1
            else:
                neg = fig.add_trace(go.Bar(x=v['negative NES'], y=wrap_prnk_sig,
                                        customdata=np.stack((v['Term'], v['FDR q-val'], wrap_leadgenes), axis=-1),
                                        orientation='h', marker_color="#636EFA",
                                        name="",
                                        hovertemplate="<b>%{customdata[0]}</b><br>" +
                                                      "|NES|: %{x}<br>"+ 
                                                      "FDR q-value: %{customdata[1]:.3f}<br>"+ 
                                                      "Leading edge genes: %{customdata[2]}",
                                        hoverlabel={'align':'left'}),
                                        row = i, col = 1)
                i += 1
        fig.update_yaxes(title="Term", tickmode='linear', tick0=0, dtick=0, automargin=True)
        fig.update_layout(title=f"Top {prerank_showX} GSEA Preranked {select_dataset} Pathways<br>({selected_col})", title_x=0.5,
                         showlegend=False,
                         yaxis={'tickmode': 'linear'},
                         font=dict(family='Arial', size=14),
                         width = 750, height = prerank_ht,
                         margin={'t':120}
                         )
        return fig

enr = Enrichr_STAGES()
prnk = Prerank_STAGES()