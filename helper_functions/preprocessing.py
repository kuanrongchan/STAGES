import pandas as pd
import numpy as np
import time
import math
import regex as re

import streamlit as st
from streamlit_tags import st_tags, st_tags_sidebar

import plotly.colors as pc
import matplotlib.pyplot as plt
import matplotlib

import decoupler as dc
# import pingouin as pg
from scipy import stats

class FC_class():
    def comparison_finder(self, cleandict):
        comparison_regex = r"(ratio|p[\.\-value]*)[_\-\s\.](.*[_\-\s\.]vs[_\-\s\.].*)"
        comparison_dict = {}
        for k,v in cleandict.items():
            comparison = list(dict.fromkeys([re.match(pattern=comparison_regex, string=i, flags=re.I).group(2) for i in v.columns if re.match(pattern=comparison_regex, string=i, flags=re.I) is not None]))
            comparison_dict[k] = comparison
        return comparison_dict

    def log_transform(self, cleandict, comparison_dict, use_corrected_pval=False):
        p_regex = "^(adj_pval)[_\-\s\.]" if use_corrected_pval else "^(p[\.\-value]*)[_\-\s\.]"
        neg_log_pval = "neg_log_adj_pval" if use_corrected_pval else "neg_log_pval"
        log_dict = {}
        for k,v in cleandict.items():
            new_combined = pd.DataFrame()
            comps_per_df = comparison_dict[k] # list of comparisons
            for comp in comps_per_df:
                ratio_pattern = f"ratio[_\-\s\.]{comp}"
                ratio_cols = [col for col in v.columns if re.match(ratio_pattern, col, flags=re.I)]
                ratio = v.loc[:,ratio_cols]

                pval_pattern = f"{p_regex}{comp}"
                pval_cols = [col for col in v.columns if re.match(pval_pattern, col, flags=re.I)]
                pval = v.loc[:,pval_cols]
                comp_df = pd.concat([np.log2(ratio), np.log10(pval)*(-1), pval], axis=1)
                comp_df.columns = [f"log2FC_{comp}", f"{neg_log_pval}_{comp}", pval.columns[0]]
                new_combined = pd.concat([new_combined, comp_df], axis=1)
            log_dict[k] = new_combined
        return log_dict


class RNAseq():

    def chunks(self, list_a, chunk_size):
        return [list_a[i:i + chunk_size] for i in range(0, len(list_a), chunk_size)]

    def violin_maxy(self, adata):
        log10adata = np.log1p(adata.to_df())
        maxy = math.ceil(max(log10adata.max(axis=0)))
        return maxy
    

    def multiviolin(self, _adata, split_long_violins):
        '''
        Parameters
        ----------
        adata: AnnData object containing counts and metadata
        split_long_violins: list | chunked list
        vthresh: int | threshold to draw the line and filter genes that are above this value
        '''
        unit_height = 3
        violin1, axes = plt.subplots(figsize = (10, len(split_long_violins) * unit_height), nrows=len(split_long_violins), ncols=1, sharey = True, constrained_layout=True)
        if len(split_long_violins) == 1:
            dc.plot_violins(_adata,
                            log=True,
                            ax=axes,
                            color="#00ABFD")
        else:
            for a, ax in zip(split_long_violins, axes):
                dc.plot_violins(_adata[a,:],
                                log = True,
                                ax = ax,
                                color = "#00ABFD")
        violin1.suptitle("Log1p counts per sample")
        return violin1, axes

    def ratio(self, adata, comp_var, baseline, against_baseline, is_log=False) -> pd.DataFrame:
        ratio_df = pd.DataFrame()
        if not is_log:
            for base in baseline:
                avg_base = adata[adata.obs[comp_var] == base].to_df().mean(axis=0) # Filter with obs first, then to_df() because otherwise it's just an array without index/columns, lastly mean
                for comp in against_baseline:
                    avg_comp = adata[adata.obs[comp_var] == comp].to_df().mean(axis=0)
                    ratio = pd.DataFrame(avg_comp / avg_base, columns = [f'ratio_{comp}_vs_{base}'])
                    ratio_df = pd.concat([ratio_df, ratio], axis=1)
        else:
            for base in baseline:
                basedf = adata[adata.obs[comp_var] == base].to_df()
                unlogbase = 2**basedf
                avg_base = unlogbase.mean(axis=0)
                for comp in against_baseline:
                    compdf = adata[adata.obs[comp_var] == comp].to_df()
                    unlogcomp = 2**compdf
                    avg_comp = unlogcomp.mean(axis=0)
                    ratio = pd.DataFrame(avg_comp / avg_base, columns = [f'ratio_{comp}_vs_{base}'])
                    ratio_df = pd.concat([ratio_df, ratio], axis=1)
        return ratio_df

    def pval_scipy(self, adata, comp_var, baseline, against_baseline, equalvar=False) -> pd.DataFrame:
        pval_df = pd.DataFrame()
        for base in baseline:
            base_ct = adata[adata.obs[comp_var] == base].to_df() # Filter with obs first, then to_df() because otherwise it's just an array without index/columns
            for comp in against_baseline:
                comp_ct = adata[adata.obs[comp_var] == comp].to_df()
                genes_pval = pd.DataFrame()
                for i in base_ct.columns:
                    T, p = stats.ttest_ind(base_ct.loc[:,i], comp_ct.loc[:,i], equal_var=equalvar, nan_policy='omit')
                    genes_pval = pd.concat([genes_pval, pd.DataFrame(data={f'pval_{comp}_vs_{base}':p}, index=[i])], axis=0)
                pval_df = pd.concat([pval_df, genes_pval], axis=1)
        return pval_df


tested = FC_class()
counts_pp = RNAseq()