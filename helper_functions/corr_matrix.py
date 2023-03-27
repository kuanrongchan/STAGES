import pandas as pd
import numpy as np
import plotly.express as px
import phik
from scipy import stats
import textwrap

class Correlation():
    def corr_matrix(self, log_dict, method):
        corr_symbols = {'spearman': 'ρ', 'pearson': 'r', 'kendall':'τ', 'phik':'𝜙k'}
        updated_df_list = []
        for k, v in log_dict.items():
            filtered = v.filter(regex="log2FC", axis=1)
            df_new = filtered.add_prefix(f"{k}_")
            df_new = df_new.loc[~df_new.index.duplicated(keep='first')]
            updated_df_list.append(df_new)

        concat_fc = pd.concat(updated_df_list, axis=1)
        concat_fc = concat_fc.dropna()
        if method != 'phik':
            concat_corr = concat_fc.corr(method=method)
        else:
            concat_corr = concat_fc.phik_matrix()
        
        # plot
        mask = np.triu(np.ones_like(concat_corr, dtype=bool))
        df_mask = concat_corr.mask(mask)
        whitespace_cols = [a.replace("_", " ") for a in df_mask.columns]
        long_tick_names = ["<br>".join(textwrap.wrap(a, 30)) for a in whitespace_cols]
        z_raw = df_mask.to_numpy().round(3)
        corr_matrix = px.imshow(z_raw, x = long_tick_names,
                                y = long_tick_names,
                                color_continuous_scale = px.colors.diverging.RdBu_r,
                                color_continuous_midpoint=0.0,
                                text_auto = True,
                                zmin=-1,
                                zmax=1,
                                labels = {"color":f"{method} {corr_symbols[method]}"})
        
        corr_matrix.update_layout(
            title_text="Comparison correlation matrix (based on log2FC values)",
            title_x=0.5,
            coloraxis_colorbar_x=0.8,
            xaxis_showgrid=False,
            yaxis_showgrid=False,
            xaxis_zeroline=False,
            yaxis_zeroline=False,
            yaxis_autorange='reversed',
            template='plotly_white',
            font=dict(
                family='Arial', size=14)
        )

        # NaN values are not handled automatically and are displayed in the figure
        # So we need to get rid of the text manually
        # for i in range(len(corr_matrix.layout.annotations)):
            # corr_matrix.layout.annotations[i].font.size = 14
            # if corr_matrix.layout.annotations[i].text == '0.000':
            #     corr_matrix.layout.annotations[i].text = ""

        return corr_matrix
    
cmatrix = Correlation()