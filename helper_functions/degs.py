# Housekeeping
from helper_functions.session_state import ss
import regex as re

# Stats and data wrangling
import pandas as pd
import numpy as np
import math

# Plotting modules
import matplotlib.pyplot as plt
import plotly.colors as pc
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.figure_factory as ff

import streamlit as st

class PreDEGs():
    '''
    This class aims to provide users with some sensing of their expression dataset before allowing them to filter by fold change and threshold.
    '''
    def nclrs(self, comparison_dict):
        combined_comparisons = set([s for v in comparison_dict.values() for s in v])
        n_comps = len(combined_comparisons)
        plotly_clrs = pc.qualitative.Plotly
        if n_comps > 10:
            colors = pc.sample_colorscale(plotly_clrs, [n/(n_comps -1) for n in range(n_comps)], colortype='tuple')
        else:
            colors = plotly_clrs[0:n_comps]
        return colors

    @st.cache_data
    def volcano(_self,
    user_log, comparison_dict,
    reset=False, xaxes = (0.0, 0.0), yaxes = 0.0,
    interactive_volcano = False,
    use_corrected_pval = False
    ):
        plt.style.use("ggplot")
        top10annotation, bottom10annotation = [], []
        
        unlist_comparisons = sorted(list(set([item for sublist in comparison_dict.values() for item in sublist])), reverse=True)
        colorlist = _self.nclrs(comparison_dict = comparison_dict)
        legend_dict = {}
        for a, c in zip(unlist_comparisons, colorlist):
            legend_dict[a.replace("_"," ").replace("-", " ")] = c

        if len(user_log) == 1:
            volcano1 = go.Figure()
            fig, ax = plt.subplots()
            highest_y = 0.0
            for k, df in user_log.items():
                comps = comparison_dict[k] # a list of comparisons made for each dataframe that the user uploads
                for i, tp in enumerate(comps):
                    complabels = tp.replace("_", " ").replace("-", " ")
                    hex_clr = legend_dict[complabels]

                    #### selecting the required FC and pval for plotting
                    pval_name = f'neg_log_adj_pval_{tp}' if use_corrected_pval else f'neg_log_pval_{tp}'
                    fc_name = f'log2FC_{tp}'
                    mini_df = df[[fc_name, pval_name]]
                    max_y = math.ceil(mini_df[pval_name].max()) + 1
                    highest_y = max_y if max_y > highest_y else highest_y
                    
                    if xaxes != (0.0, 0.0) and yaxes != (0.0):
                        user_filter = mini_df[(mini_df[pval_name] <= yaxes) & (mini_df[fc_name].between(xaxes[0], xaxes[1],inclusive='both'))]
                    elif xaxes == (0.0, 0.0) and yaxes != (0.0):
                        user_filter = mini_df[(mini_df[pval_name] <= yaxes)]

                    elif xaxes != (0.0, 0.0) and yaxes == (0.0):
                        user_filter = mini_df[(mini_df[fc_name].between(xaxes[0], xaxes[1], inclusive='both'))]
                    else:
                        user_filter = mini_df
                    
                    top_10 = user_filter.sort_values(by=fc_name, ascending=False).head(10)
                    bottom_10 = user_filter.sort_values(by=fc_name, ascending=False).tail(10)
                    
                    bottom10annotation.append(
                        bottom_10.rename(columns={fc_name: "log2FC", pval_name: "neg_log_pval"}))
                    top10annotation.append(
                        top_10.rename(columns={fc_name: "log2FC", pval_name: "neg_log_pval"}))

                    ax.grid(visible=True, which="major", axis="both", alpha=0.5)
                    plt.scatter(user_filter[fc_name], user_filter[pval_name], alpha=0.8, label=complabels, c = [hex_clr])
                    plt.title("Volcano plot across comparisons", loc='center')
                    plt.xlabel('log2(Fold-change)')
                    plt.ylabel('-log10(p-value)')
                    plt.axhline(y=0, color='r', linestyle='dashed')
                    plt.axvline(x=0, linestyle='dashed')

                    leg = ax.legend(bbox_to_anchor=(1, 1), loc='upper left', frameon=False, borderaxespad=0.3)

                    if xaxes != (0.0,0.0):
                        plt.xlim([xaxes[0], xaxes[1]])
                    else:
                        pass

                    if yaxes != (0.0):
                        plt.ylim(-0.50, yaxes)
                    else:
                        plt.ylim(-0.50, highest_y)

                    if interactive_volcano:
                        volcano1.add_trace(go.Scatter(x=user_filter[fc_name], y=user_filter[pval_name],
                                                    mode='markers',
                                                    name = complabels, hovertext=list(user_filter.index),
                                                    marker=dict(color=hex_clr, size=8), opacity=0.9,
                                                    legendgroup=tp
                                                    )
                                        )
                annotationconcat_top = pd.concat(top10annotation, axis=0)
                annotationconcat_top = annotationconcat_top.sort_values(by=["log2FC"], ascending=False).head(10)

                annotationconcat_bottom = pd.concat(bottom10annotation, axis=0)
                annotationconcat_bottom = annotationconcat_bottom.sort_values(by=["log2FC"], ascending=True).head(10)

                # in case I forget for annotation:
                # 1. get the top 10 of each comparison and rename the annotation cols to log2FC and negative_log_pval (generic)
                # 2. concat all the top 10s for each comparison along axis=0 and sort by highest log2FC
                # 3. the result is a long df of 2 cols (log2FC, neglogpval), where my xy annotation coordinates is x = log2FC (1st col), y = neglogpval (2nd col)

                for i in range(len(annotationconcat_top)):
                    plt.annotate(text=annotationconcat_top.index[i], xy=(annotationconcat_top.iloc[i, 0], annotationconcat_top.iloc[i, 1]),
                                xytext=(0, 3), horizontalalignment='center', textcoords='offset points',
                                fontsize=7)  # top 10
                for i in range(len(annotationconcat_bottom)):
                    plt.annotate(text=annotationconcat_bottom.index[i],
                                xy=(annotationconcat_bottom.iloc[i, 0], annotationconcat_bottom.iloc[i, 1]),
                                xytext=(0, 3), horizontalalignment='center', textcoords='offset points',
                                fontsize=7)  # bottom 10
            names = set()
            volcano1.for_each_trace(lambda trace:trace.update(showlegend=False)if trace.name in names else names.add(trace.name))

            volcano1.update_layout(showlegend=True,
                                title="Interactive volcano across comparisons",
                                legend_title_text="Comparisons",
                                font=dict(family='Arial', size=14),
                                xaxis_title="log2(Fold-change)",yaxis_title="-log10(p-value)")

            if xaxes != (0.0,0.0):
                volcano1.update_xaxes(range=[xaxes[0], xaxes[1]])
            else:
                pass
            
        else:
            i = 1
            if len(user_log) % 2 == 0:
                nrows = math.ceil(len(user_log) / 2)
                extras = nrows*2 - len(user_log)
                volcano1 = make_subplots(rows=nrows, cols=2, subplot_titles=(list(user_log.keys())),
                                        x_title="log2(Fold-Change)", y_title="-log10(p-value)", shared_xaxes=True, shared_yaxes=True)
                v_row, v_col = 1, 1
                j = 1
                fig, axs = plt.subplots(nrows=nrows, ncols=2, sharex=True, sharey = True, figsize=(8, 7))

            else:
                nrows = math.ceil(len(user_log) / 3)
                extras = nrows*3 - len(user_log)
                volcano1 = make_subplots(rows=nrows, cols=3, subplot_titles=(list(user_log.keys())),
                                        x_title="log2(Fold-Change)", y_title="-log10(p-value)", shared_xaxes=True, shared_yaxes=True)
                v_row, v_col = 1, 1
                j = 1
                fig, axs = plt.subplots(nrows=nrows, ncols=3, sharex=True, sharey=True, figsize=(8,7))

            min_x, max_x, max_y = 0,0,0
            for k, df in user_log.items():
                comps = comparison_dict[k] # a list of comparisons made for each dataframe that the user uploads
                for i, tp in enumerate(comps):
                    complabels = tp.replace("_", " ").replace("-", " ")
                    hex_clr = legend_dict[complabels]
                    #### selecting the required FC and pval for plotting
                    pval_name = f'negative_log_pval_{tp}'
                    fc_name = f'log2FC_{tp}'
                    mini_df = df[[fc_name, pval_name]]

                    min_x = math.floor(mini_df[fc_name].min()) if math.floor(mini_df[fc_name].min()) < min_x else min_x
                    max_x = math.ceil(mini_df[fc_name].max()) if math.ceil(mini_df[fc_name].max()) > max_x else max_x
                    max_y = math.ceil(mini_df[pval_name].max()) if math.ceil(mini_df[pval_name].max()) > max_y else max_y

                    if xaxes != (0.0, 0.0) and yaxes != (0.0):
                        user_filter = mini_df[(mini_df[pval_name] <= yaxes) & (mini_df[fc_name].between(xaxes[0], xaxes[1],inclusive='both'))]
                    elif xaxes == (0.0, 0.0) and yaxes != (0.0):
                        user_filter = mini_df[(mini_df[pval_name] <= yaxes)]

                    elif xaxes != (0.0, 0.0) and yaxes == (0.0):
                        user_filter = mini_df[(mini_df[fc_name].between(xaxes[0], xaxes[1], inclusive='both'))]
                    else:
                        user_filter = mini_df

                    top_10 = user_filter.sort_values(by=fc_name, ascending=False).head(10)
                    bottom_10 = user_filter.sort_values(by=fc_name, ascending=True).head(10)

                    bottom10annotation.append(
                        bottom_10.rename(columns={fc_name: "log2FC", pval_name: "negative_log_pval"}))
                    top10annotation.append(
                        top_10.rename(columns={fc_name: "log2FC", pval_name: "negative_log_pval"}))

                    ax = plt.subplot(nrows, 2, j) if len(user_log) % 2 == 0 else plt.subplot(nrows, 3, j)
                    ax.grid(visible=True, which="major", axis="both", alpha=0.5)
                    ax.scatter(user_filter[fc_name], user_filter[pval_name], alpha=0.9, label = complabels, c = [hex_clr])
                    ax.axhline(y=0, color='r', linestyle='dashed')
                    ax.axvline(x=0, linestyle='dashed')
                    ax.set_title(f"{k}", fontdict={'fontsize':10})

                    if xaxes != (0.0,0.0):
                        ax.set_xlim([xaxes[0], xaxes[1]])
                    else:
                        ax.set_xlim([min_x, max_x])

                    if yaxes != (0.0):
                        ax.set_ylim([-1.0, yaxes])
                    else:
                        ax.set_ylim([-1.0, max_y])

                    if interactive_volcano:
                        volcano1.add_trace(go.Scatter(x=user_filter[fc_name], y=user_filter[pval_name],
                                                    mode='markers',
                                                    name=complabels, hovertext=list(user_filter.index),
                                                    marker=dict(color=hex_clr, size=8, opacity=0.9), legendgroup=str(i)),
                                        row=v_row, col=v_col)
                        i += 1

                annotationconcat_top = pd.concat(top10annotation, axis=0)
                annotationconcat_top = annotationconcat_top.sort_values(by=["log2FC"], ascending=False).head(10)

                annotationconcat_bottom = pd.concat(bottom10annotation, axis=0)
                annotationconcat_bottom = annotationconcat_bottom.sort_values(by=["log2FC"], ascending=True).head(10)

                for i in range(len(annotationconcat_top)):
                    plt.annotate(text=annotationconcat_top.index[i], xy=(annotationconcat_top.iloc[i, 0], annotationconcat_top.iloc[i, 1]),
                                xytext=(0, 3), horizontalalignment='center', textcoords='offset points',
                                fontsize=7)  # top 10
                for i in range(len(annotationconcat_bottom)):
                    plt.annotate(text=annotationconcat_bottom.index[i],
                                xy=(annotationconcat_bottom.iloc[i, 0], annotationconcat_bottom.iloc[i, 1]),
                                xytext=(0, 3), horizontalalignment='center', textcoords='offset points',
                                fontsize=7)  # bottom 10

                top10annotation.clear()
                bottom10annotation.clear()

                i = 1
                j += 1
                v_col += 1

                if (len(user_log) % 2 == 0) and v_col > 2:
                    v_col = 1
                    v_row += 1
                if (len(user_log) % 2 != 0) and v_col > 3:
                    v_col = 1
                    v_row += 1
            
            handles = {}
            for a in fig.axes:
                hn, lb = a.get_legend_handles_labels()
                for h,l in zip(hn, lb):
                    handles[l] = h
    
            leg = plt.legend(handles = handles.values(), labels = handles.keys(), bbox_to_anchor=(1.1, 1), loc='upper left', frameon=False, borderaxespad=0.3)

            if extras == 1:
                axs[nrows-1, 2].remove()
            elif extras == 2:
                axs[nrows-1, 2].remove()
                axs[nrows-1, 1].remove()
            else:
                pass
            
            fig.add_subplot(111, frame_on=False)
            plt.grid(visible=False)
            plt.tick_params(labelcolor="none", bottom=False, left=False)
            fig.suptitle("Volcano plot across comparisons", fontsize=14)
            plt.xlabel("log2(Fold-change)")
            plt.ylabel("-log10(p-value)")

            if xaxes != (0.0,0.0):
                volcano1.update_xaxes(range=[xaxes[0], xaxes[1]])
            else:
                pass
            ## removing duplicate legends
            names = set()
            volcano1.for_each_trace(lambda trace: trace.update(showlegend=False) if (trace.name in names) else names.add(trace.name))

            volcano1.update_layout(showlegend=True,
                                title="Interactive volcano across comparisons", title_x=0.5,
                                legend_title_text="Comparisons",
                                font=dict(family='Arial', size=14)
                                )
        return fig, volcano1
    
    @st.cache_data
    def deg_cdf(_self, ready_dict, comparison_dict, pval=0.05, markermode='lines', use_corrected_pval=False):
        FC_step = [i/10 for i in range(0, 31, 1)]
        fig = go.Figure()
        for k,v in ready_dict.items():
            comparisons = comparison_dict[k]
            for comp in comparisons:
                df = v.filter(regex=comp, axis=1)
                # Get pvals
                p = [i for i in df.columns if re.search("^pval", i)][0] if not use_corrected_pval else [i for i in df.columns if re.search("^adj_pval", i)][0]
                p = df.loc[:,p]
                ratio = df.filter(regex="^ratio", axis=1)
                fc = pd.DataFrame(index = ratio.index, columns=[f'FC_{comp}'])
                for r in range(len(ratio.index)):
                    value = ratio.iloc[r,0]
                    fc.iloc[r, 0] = value if value >=1 else (-1)/value

                p_fc = pd.concat([p, fc], axis=1)
                n_p = [len(p_fc[(p_fc.iloc[:,0] < pval) & (p_fc.iloc[:,1] > a)]) for a in FC_step]
                n_n = [len(p_fc[(p_fc.iloc[:,0] < pval) & (p_fc.iloc[:,1] < -a)]) for a in FC_step]
                n_total = [x + y for x, y in zip(n_p, n_n)]
                fig.add_trace(go.Scatter(
                    x=FC_step, y=n_total, name=comp.replace("_"," "),
                    hovertemplate=f"p-value: {pval}<br>FC: %{{x}}<br>number of DEGs: %{{y}}",
                    mode=markermode))

        fig.update_layout(title="DE cutoff plot of DEG numbers at various points",
                        title_x = 0.5,
                        xaxis=dict(title="Fold change"),
                        yaxis=dict(title="Number of DEGs"),
                        template="ggplot2",
                        legend=dict(title="Comparison"),
                        width=700, height=500)
        return fig
    

class DEGs():
    '''
    This class will provide the output for bar plots and data containing DEGs.
    '''

    @st.cache_data
    def degs(_self, log_ready_dict, comparison_dict, pval_cutoff=0.0, fc_cutoff=0.0, u_width = 800, u_height=600, use_corrected_pval=False):
        log2fc_cutoff = np.log2(fc_cutoff)
        p_format = "adjusted p-value" if use_corrected_pval else "p-value"
        ####################################### Filter DF by Pvals and FC #################################################
        proportions, deg_dict = {}, {}
        for k,v in log_ready_dict.items():
            comps = comparison_dict[k]
            up, down = [], []
            for cmp in comps:
                pval_name = f"adj_pval_{cmp}" if use_corrected_pval else f"pval_{cmp}"
                logfc_name = f"log2FC_{cmp}"
                filtered = v[[pval_name, logfc_name]]
                degs = filtered[(filtered[pval_name] < pval_cutoff) & ((filtered[logfc_name] > log2fc_cutoff)|(filtered[logfc_name] < -log2fc_cutoff))]
                deg_dict[f"{k}_{cmp}"] = degs
                # calculating proportion of DEGs for pie chart
                upreg_deg = degs[degs.loc[:, logfc_name] > 0]
                up.append(len(upreg_deg))
                proportions[f"{k}_upcount"] = up
                proportions[f"UP_{k}_{cmp}"] = upreg_deg

                downreg_deg = degs[degs.loc[:, logfc_name] < 0]
                down.append(len(downreg_deg))
                proportions[f"{k}_downcount"] = down
                proportions[f"DOWN_{k}_{cmp}"] = downreg_deg
                
        if len(log_ready_dict) == 1:
            stacked1 = go.Figure()
        else:
            stacked_row = 1
            stacked_col = 1
            if len(log_ready_dict) % 2 == 0:
                nrows = math.ceil(len(log_ready_dict) / 2)
                stacked1 = make_subplots(rows=nrows, cols=2, subplot_titles=(list(log_ready_dict.keys())),
                                        y_title='Number of DEGs',
                                        vertical_spacing = 0.5, shared_yaxes=True)
            else:
                nrows = math.ceil(len(log_ready_dict) / 3)
                stacked1 = make_subplots(rows=nrows, cols=3, subplot_titles=(list(log_ready_dict.keys())),
                                        y_title='Number of DEGs', 
                                        vertical_spacing=0.5, horizontal_spacing=0.02,
                                        shared_yaxes=True)
            
        for k,v in log_ready_dict.items():
            if len(log_ready_dict) == 1:
                # Stacked Bar
                stacked1.add_trace(
                    go.Bar(x=comps, y=proportions[f'{k}_downcount'], name="Downregulated", marker_color="#636EFA"))
                stacked1.add_trace(
                    go.Bar(x=comps, y=proportions[f'{k}_upcount'], name="Upregulated", marker_color="#EF553B"))
            else:
                # Stacked Bar
                stacked1.add_trace(
                    go.Bar(x=comps, y=proportions[f'{k}_downcount'], name="Downregulated", marker_color="#636EFA",
                            legendgroup="A"),
                    row=stacked_row, col=stacked_col)
                stacked1.add_trace(
                    go.Bar(x=comps, y=proportions[f'{k}_upcount'], name="Upregulated", marker_color="#EF553B",
                            legendgroup="B"),
                    row=stacked_row, col=stacked_col)

                stacked_col += 1
                if len(log_ready_dict) % 2 == 0 and stacked_col > 2:
                    stacked_col = 1
                    stacked_row += 1
                elif len(log_ready_dict) % 2 != 0 and stacked_col > 3:
                    stacked_col = 1
                    stacked_row += 1

                ## removing duplicate legends
                names = set()
                stacked1.for_each_trace(lambda trace:trace.update(showlegend=False) if (trace.name in names) else names.add(trace.name))
        
        stacked1.update_layout(showlegend=True, barmode='stack',
                            title=f"Number of DEGs across comparisons<br>(FC {fc_cutoff}; {p_format} {pval_cutoff})",
                            title_x=0.5,
                            legend_title_text='DEGs:',
                            font=dict(family='Arial', size=14), width=u_width, height=u_height)
        stacked1.update_xaxes(automargin=True)

        proportions = {k:v for k,v in proportions.items() if type(v) != list}
        return stacked1, proportions


preDE = PreDEGs()
DE = DEGs()