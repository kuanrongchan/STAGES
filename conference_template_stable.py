#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import time
import math
import csv
import re
# import base64
from io import BytesIO

import gseapy as gp
from scipy.spatial.distance import pdist, squareform
from scipy import stats

import streamlit as st
from streamlit_tags import st_tags, st_tags_sidebar

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.figure_factory as ff
import plotly.express as px
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar


# Stable version 1

################################################ for df download #######################################################
def convert_df(df):
    return df.to_csv().encode('utf-8')


def to_excel(df):
    output = BytesIO()
    writer = pd.ExcelWriter(output, engine='xlsxwriter')
    for d, i in zip(df, range(len(df))):
        d.to_excel(writer, sheet_name=f'Sheet {i + 1}')
    writer.save()
    processed_data = output.getvalue()
    return processed_data


# def get_table_download_link(df): keeping just in case download button fails
#     """Generates a link allowing the data in a given panda dataframe to be downloaded
#     in:  dataframe
#     out: href string
#     """
#     val = to_excel(df)
#     b64 = base64.b64encode(val)
#     return f'<a href="data:application/octet-stream;base64,{b64.decode()}" download="deg_files.xlsx">Download DEGs as Excel File</a>' # decode b'abc' => abc

st.title("Omics Dashboard \U0001F4CA")

################################################# Documentation ########################################################
documentation = st.sidebar.checkbox("Read the Docs", value=False, key='documentation')
################################################# File Uploader ########################################################
df_query = st.sidebar.file_uploader(
    'Upload your .csv/.xlsx file. A demo dataset will be uploaded if no files are uploaded',
    accept_multiple_files=True)

df_dict = {}
df_names = []

if len(df_query) != 0:
    for d in df_query:
        head, sep, tail = str(d.name).partition(".")
        if tail == 'csv':
            data = st.experimental_memo(pd.read_csv)(d, index_col=0)
            df_dict[head] = data
            df_names.append(head)

        elif tail == 'xlsx':
            x = st.experimental_memo(pd.read_excel)(d, index_col=0, sheet_name=None, engine='openpyxl')
            selected_sheet = st.sidebar.selectbox(label="* Select which sheet to read in", options=x.keys())
            data = x[selected_sheet]
            df_dict[head] = data
            df_names.append(head)
else:
    testdata = st.experimental_memo(pd.read_csv)("demo_dataframe_corrected.csv", index_col=0)
    testname = "Demo"
    df_dict[testname] = testdata
    df_names.append(testname)

if st.sidebar.checkbox("Show uploaded/demo dataframe"):
    for k, v in df_dict.items():
        st.write(f"**{k} dataframe**", v)

######## Important ########
deg_dict = {}      ########
proportions = {}   ########
###########################


################################################# Read the Docs #######################################################
def read_docs():
    st.subheader("Static and Temporal Analysis of Gene Expression Studies (STAGES) documentation")
    st.image("https://user-images.githubusercontent.com/91276553/141084154-7d84695a-b220-43c5-bd41-08a38fd0ec70.png",
             width=None)
    st.markdown(
        '''
    STAGES is a multi-app that integrates data visualisation and pathway analysis for static and temporal gene expression studies. STAGES is an open source and community funded web tool for creating beautiful charts from gene expression datasets. The multi-page web app built using Streamlit, which currently allows users to analyse an omics data in distinct stages:

    1. Plot interactive volcano plots
    2. Filter data for differentially expressed genes (Users can apply their preferred fold-change and p-value cut-offs to identify DEG number and identities)
    3. Build customised clustergrams based on identified up-regulated DEGs (UP) or down-regulated DEGs (DOWN)
    4. Build customised clustergrams based on user-selected gene list
    5. Perform Enrichr analysis based on DEGs or user-selected gene list
    6. Plot interactive correlation matrix comparing across different time-points or experimental conditions

    ## Getting started

    To use the app, you will need one comparison file which should minimally contain:
    
    1. Gene names on the first column
    2. Ratio values (relative transcript expression versus control or baseline)
    3. Adjusted p-value (or p-value)

    For the app to be able to recognise your ratio and p-values, please label:

    1. Ratio as ratio_X_vs_Y
    2. Adjusted p-values (or p-values) as pval_X_vs_Y,

    where X and Y are the comparison variables. 

    Some examples of labelling "X" include: ratio_virus_vs_ctrl, ratio_drugA_vs_placebo, ratio_hr6_vs_0, ratio_day1_vs_day0. 
    
    Some examples of labelling "Y" include: pval_virus_vs_ctrl, pval_drugA_vs_placebo, pval_hr6_vs_0, pval_day1_vs_day0. 

    For multiple comparisons to be made within the same graph, simply insert more comparison columns (e.g. ratio_A_vs_Y, pval_A_vs_Y, ratio_B_vs_Y, pval_B_vs_Y ...), but please ensure that  "Y" is consistently present in all comparisons. Also, ensure that no icons or symbols used for labelling "X" and "Y." If you have other column statistics, it is not necessary to remove them.

    To perform multiple comparisons for time-course experiments, you can choose to upload multiple .csv or .xls files. But please do ensure that the header columns are labelled the same way (meaning that the data has to measured at same time-points for the different experimental conditions)

    Demo examples are provided. You can try out the demo examples to familiarise yourself with the apps before uploading your dataset

    ## Data safety and security
    The data you upload is safe and is never stored anywhere.

    ## Contributors
    These apps are jointly made by myself (Kuan Rong Chan), Clara Koh, Justin Ooi and Gabrielle Lee from Duke-NUS, Department of Emerging Infectious Diseases. I am also thankful for Eugenia Ong and Ayesa Syenina from VIREMICS for their constructive feedback. These apps are now free for everyone to use, but for a limited period of time as we are constantly upgrading the apps. For more details on what we do, feel free to visit us at [omicsdiary.com](https://omicsdiary.com/).

        ''')
if documentation:
    read_docs()

################################################## Grabs Timepoint ####################################################
# Is a variable to allow the following functions to render either timepoint or comparison in plot titles
is_tp = 1


def timepoints(df_dict):
    pattern_days = "[-_\s]{1}(min[0-9]*[.]?[0-9]+)|([0-9]*[.]?[0-9]+min)[-_\s]?|(hr[0-9]*[.]?[0-9]+)|([0-9]*[.]?[0-9]+hr?)|(D[ay]*[0-9]*[.]?[0-9]+)|([0-9]*[.]?[0-9]+D[ay]*)"
    tpcomps = []
    for k, df in df_dict.items():
        extract_days = df.columns.str.extract(pattern_days, flags=re.IGNORECASE)
        for num in range(6):  # 6 capture groups (0-5 but range including 5 has to +1)
            temp = extract_days[num].unique().tolist()
            for days in temp:
                if days not in tpcomps:
                    tpcomps.append(days)
        tpcomps_nangone = [x for x in tpcomps if x == x]
        if len(tpcomps_nangone) == 0:
            global is_tp
            is_tp = 0

            pattern_comps = "(?<!adj)[-_\s]\w*\d*[-_\s]vs[-_\s]\w*\d*"
            findcomps = re.findall(pattern_comps, str(df.columns.tolist()), flags=re.I)
            for c in findcomps:
                cln = c.replace("_", "", 1)
                if cln not in tpcomps:
                    tpcomps.append(cln)
            tpcomps_nangone = [x for x in tpcomps if str(x) != 'nan']
    return tpcomps_nangone


######################################### Check for log2FC columns only ################################################
def check_log(df_dict):
    df_dict2 = {}
    for k, df in df_dict.items():
        df = df.apply(pd.to_numeric, errors='coerce')
        fc_regex = "log2Fold[-_\s]?Changes?|log2FC"
        cols = df.columns.tolist()
        string_cols = ' '.join(cols)
        check_fc = re.findall(fc_regex, string_cols, flags=re.I)

        if check_fc:
            d_pvals = df.filter(regex=re.compile("p[-_\s]?valu?e?", re.I))
            d_neglogpval = np.log10(d_pvals) * (-1)
            d_neglogpval.columns = d_neglogpval.columns.str.replace(
                "p[-_\s]?valu?e?", "negative_log_pval", regex=True,
                flags=re.IGNORECASE)
            df2 = pd.concat([df, d_neglogpval], axis=1)
            df_dict2[f"{k}"] = df2

        else:
            d_ratios = df.filter(regex=re.compile("ratio", re.I))
            d_log2FC = np.log2(d_ratios)
            d_log2FC.columns = d_log2FC.columns.str.replace("ratio", "log2FC", regex=True, flags=re.IGNORECASE)

            d_pvals = df.filter(regex=re.compile("p[-_\s]?valu?e?", re.I))
            d_neglogpval = np.log10(d_pvals) * (-1)
            d_neglogpval.columns = d_neglogpval.columns.str.replace(
                "p[-_\s]?valu?e?", "negative_log_pval", regex=True,
                flags=re.IGNORECASE)
            df2 = pd.concat([df, d_log2FC, d_neglogpval], axis=1)
            df_dict2[f"{k}"] = df2
    return df_dict2


######################################### Colours to use for interactive volcano plots ################################
def n_colors(list_of_days):
    len_days = len(list_of_days)
    colors = px.colors.qualitative.Plotly[0:len_days + 1]
    return colors


###################################################### Volcano Plot ###################################################
def volcano(df_dict, list_of_days, colorlist):
    st.subheader("Volcano plot analysis")
    vol_expand = st.sidebar.expander("Expand for volcano plot", expanded=False)
    interactive_volcano = vol_expand.checkbox(label="Show interactive volcano plot", value=False,
                                              help="Facilitates gene name display on hover. This may cause lag")
    if is_tp == 1:
        tp_or_comp = "time-points (vs baseline)"
    else:
        tp_or_comp = "comparisons"

    if len(df_dict) == 1:
        volcano1 = go.Figure()
        fig = plt.figure()
        for k, df in df_dict.items():
            for tp, col in zip(list_of_days, colorlist):
                complabels = tp.replace("_", " ").replace("-", " ")
                #### selecting the required FC and pval for plotting
                FC_col_name = list([col for col in df.columns if tp in col if "log2FC" in col])
                fold_changes = df[FC_col_name[0]]
                pval_col_name = list([col for col in df.columns if tp in col if "negative_log_pval" in col])
                # edited to negative instead of neg so we can pass all the dfs through data formatter
                # (to include the log2ratio and -log pval)
                pvals = df[pval_col_name[0]]
                plt.grid(b=True, which="major", axis="both", alpha=0.3)
                plt.scatter(fold_changes, pvals, alpha=0.7, label=complabels)
                plt.title(f"Volcano plot across {tp_or_comp}", loc='center')
                plt.legend(bbox_to_anchor=(1.04, 1), loc='upper left')
                plt.xlabel('log2(Fold-change)')
                plt.ylabel('-log10(p-value)')
                plt.axhline(y=0, color='r', linestyle='dashed')
                plt.axvline(x=0, linestyle='dashed')

                if interactive_volcano:
                    volcano1.add_trace(go.Scatter(x=fold_changes, y=pvals,
                                                  mode='markers',
                                                  name=complabels, hovertext=list(df.index),
                                                  line=dict(color=col)
                                                  )
                                       )
        volcano1.update_layout(showlegend=True,
                               title=f"Interactive volcano across {tp_or_comp}",
                               legend_title_text="Timepoint",
                               font=dict(family='Arial', size=14),
                               xaxis_title="log2(Fold-change)", yaxis_title="-log10(p-value)")
    else:
        if len(df_dict) % 2 == 0:
            nrows = math.ceil(len(df_dict) / 2)
            volcano1 = make_subplots(rows=nrows, cols=2, subplot_titles=(list(df_dict.keys())),
                                     x_title="log2(Fold-Change)", y_title="-log10(p-value)")
            v_row, v_col = 1, 1
            j = 1
            fig, axes = plt.subplots(nrows=nrows, ncols=2, sharex=True, sharey=True)
        elif math.ceil(len(df_dict) % 3) == 0:
            nrows = math.ceil(len(df_dict) / 2)
            volcano1 = make_subplots(rows=nrows, cols=2, subplot_titles=(list(df_dict.keys())),
                                     x_title="log2(Fold-Change)", y_title="-log10(p-value)")
            v_row, v_col = 1, 1
            j = 1
            fig, axes = plt.subplots(nrows=nrows, ncols=3, sharex=True, sharey=True)

        for k, df, in df_dict.items():
            for tp, clr in zip(list_of_days, colorlist):
                complabels = tp.replace("_", " ").replace("-", " ")
                #### selecting the required FC and pval for plotting
                FC_col_name = list([col for col in df.columns if tp in col if "log2FC" in col])
                fold_changes = df[FC_col_name[0]]
                pval_col_name = list([col for col in df.columns if tp in col if "negative_log_pval" in col])
                # edited to negative instead of neg so we can pass all the dfs through data formatter
                # (to include the log2ratio and -log pval)
                pvals = df[pval_col_name[0]]
                if len(df_dict) % 2 == 0:
                    plt.subplot(nrows, 2, j)
                elif len(df_dict) % 3 == 0:
                    plt.subplot(nrows, 3, j)
                plt.grid(b=True, which="major", axis="both", alpha=0.3)
                plt.scatter(fold_changes, pvals, alpha=0.7, label=complabels)
                plt.axhline(y=0, color='r', linestyle='dashed')
                plt.axvline(x=0, linestyle='dashed')

                if interactive_volcano:
                    volcano1.add_trace(go.Scatter(x=fold_changes, y=pvals,
                                                  mode='markers',
                                                  name=complabels, hovertext=list(df.index),
                                                  line=dict(color=clr)
                                                  ),
                                       row=v_row, col=v_col
                                       )
            j += 1
            v_col += 1
            if (len(df_dict) % 2 == 0) and j > 2:
                j = 1
            if (len(df_dict) % 3 == 0) and j > 3:
                j = 1

            if (len(df_dict) % 2 == 0) and v_col > 2:
                v_col = 1
                v_row += 1
            if (len(df_dict) % 3 == 0) and v_col > 3:
                v_col = 1
                v_row += 1

        plt.legend(bbox_to_anchor=(1.04, 1), loc='upper left')
        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor="none", bottom=False, left=False)
        plt.title(f"Volcano plot across {tp_or_comp}", loc='center')
        plt.xlabel('log2(Fold-change)')
        plt.ylabel('-log10(p-value)')
        plt.tight_layout(h_pad=1.0)

        ## removing duplicate legends
        names = set()
        volcano1.for_each_trace(
            lambda trace:
            trace.update(showlegend=False)
            if (trace.name in names) else names.add(trace.name))

        volcano1.update_layout(showlegend=True,
                               title=f"Interactive volcano across {tp_or_comp}",
                               legend_title_text="Timepoint",
                               font=dict(family='Arial', size=14)
                               )
    if interactive_volcano:
        ivolcano = st.success("Plot complete!")
        time.sleep(0.25)
        ivolcano.empty()
        st.pyplot(fig, use_container_width=True)
        st.plotly_chart(volcano1, use_container_width=True)
    else:
        svolcano = st.success("Plot complete!")
        time.sleep(0.25)
        svolcano.empty()
        st.pyplot(fig)


############################################# Stacked DEGs ############################################################
def degs(df_dict, list_of_days, colorlist):
    st.subheader('Differential expression analysis')
    stack_expand = st.sidebar.expander("Expand for stacked bar chart", expanded=False)
    with stack_expand:
        pval_cutoff = st.slider(label=f"Adjust p-value cutoff here", value=0.05, min_value=0.0, max_value=1.0,
                                step=0.01)
        fc_cutoff = st.slider(label="Adjust fold-change cutoff here ", min_value=0.0, max_value=5.0, value=2.0,
                              step=0.1)
        log2fc_cutoff = np.log2(fc_cutoff)

    if is_tp == 1:
        tp_or_comp = "time"
    else:
        tp_or_comp = "comparisons"
    ####################################### Filter DF by Pvals and FC #################################################

    pval_regex = "p[-_\s]?valu?e?"
    fc_regex = "log2Fold[-_\s]?Changes?|log2FC"

    new_dfs = {}

    def organised_df():
        df_list = [df_dict[x] for x in df_names]
        for name, df in zip(df_names, df_list):
            pvalcols = df.filter(regex=re.compile(pval_regex, flags=re.I))
            pvalcols.columns = pvalcols.columns.str.replace(pval_regex, "pval", regex=True)
            fccols = df.filter(regex=re.compile(fc_regex, flags=re.I))
            fccols.columns = fccols.columns.str.replace(fc_regex, "log2FC", regex=True)
            newdf = pd.concat([pvalcols, fccols], axis=1)
            new_dfs[name] = newdf

    organised_df()
    newdf_list = [new_dfs[x] for x in df_names]  # only contains adj pval and log2FC

    upreg_deg_count = []
    downreg_deg_count = []

    if len(df_dict) == 1:
        stacked1 = go.Figure()
        for df, name in zip(newdf_list, df_names):
            for tp, clrs in zip(list_of_days, colorlist):
                pval_name_pertp = list([col for col in df.columns if f"pval_{tp}" in col])
                pval_col_pertp = df[pval_name_pertp[0]]
                pval_filter_pertp = (pval_col_pertp[(pval_col_pertp <= pval_cutoff)]).to_frame()

                fc_name_pertp = list([col for col in df.columns if "log2FC_{}".format(tp) in col])
                fc_col_pertp = df[fc_name_pertp[0]]
                fc_filter_pertp = (
                fc_col_pertp[((fc_col_pertp >= log2fc_cutoff) | (fc_col_pertp < -(log2fc_cutoff)))]).to_frame()

                pval_fc_pertp = pval_filter_pertp.merge(fc_filter_pertp, how='inner', left_index=True, right_index=True)
                deg_dict[f"{name}_{tp}"] = pval_fc_pertp

                # calculating proportion of DEGs for pie chart
                upreg_deg = pval_fc_pertp[pval_fc_pertp.iloc[:, 1] > 0]
                upreg_deg_count.append(len(upreg_deg))
                proportions["upcount"] = upreg_deg_count
                proportions[f"UP_{name}_{tp}"] = upreg_deg

                downreg_deg = pval_fc_pertp[pval_fc_pertp.iloc[:, 1] < 0]
                downreg_deg_count.append(len(downreg_deg))
                proportions["downcount"] = downreg_deg_count
                proportions[f"DOWN_{name}_{tp}"] = downreg_deg

                # Stacked Bar

            stacked1.add_trace(
                go.Bar(x=list_of_days, y=proportions['downcount'], name="Downregulated", marker_color="#636EFA"))
            stacked1.add_trace(
                go.Bar(x=list_of_days, y=proportions['upcount'], name="Upregulated", marker_color="#EF553B"))
            upreg_deg_count.clear()
            downreg_deg_count.clear()

        stacked1.update_layout(showlegend=True, barmode='stack',
                               title=f"Number of DEGs across {tp_or_comp} (based on selected cutoffs)",
                               xaxis_title=tp_or_comp.title(), yaxis_title="Number of DEGs",
                               legend_title_text='DEGs:',
                               font=dict(
                                   family='Arial', size=14))
    else:
        if len(df_dict) % 2 == 0:
            nrows = math.ceil(len(df_dict) / 2)
            stacked1 = make_subplots(rows=nrows, cols=2, subplot_titles=(list(df_dict.keys())),
                                     x_title=tp_or_comp.title(), y_title='Number of DEGs', shared_yaxes=True)
        elif math.ceil(len(df_dict)) % 3 == 0:
            nrows = math.ceil(len(df_dict) / 3)
            stacked1 = make_subplots(rows=row_no, cols=3, subplot_titles=(list(df_dict.keys())),
                                     x_title=tp_or_comp.title(), y_title='Number of DEGs')
        stacked_row = 1
        stacked_col = 1
        for df, name in zip(newdf_list, df_names):
            for tp, clrs in zip(list_of_days, colorlist):
                pval_name_pertp = list([col for col in df.columns if f"pval_{tp}" in col])
                pval_col_pertp = df[pval_name_pertp[0]]
                pval_filter_pertp = (pval_col_pertp[(pval_col_pertp <= pval_cutoff)]).to_frame()

                fc_name_pertp = list([col for col in df.columns if "log2FC_{}".format(tp) in col])
                fc_col_pertp = df[fc_name_pertp[0]]
                fc_filter_pertp = (
                fc_col_pertp[((fc_col_pertp >= log2fc_cutoff) | (fc_col_pertp < -(log2fc_cutoff)))]).to_frame()

                pval_fc_pertp = pval_filter_pertp.merge(fc_filter_pertp, how='inner', left_index=True, right_index=True)
                deg_dict[f"{name}_{tp}"] = pval_fc_pertp

                # calculating proportion of DEGs for pie chart
                upreg_deg = pval_fc_pertp[pval_fc_pertp.iloc[:, 1] > 0]
                upreg_deg_count.append(len(upreg_deg))
                proportions["upcount"] = upreg_deg_count
                proportions[f"UP_{name}_{tp}"] = upreg_deg

                downreg_deg = pval_fc_pertp[pval_fc_pertp.iloc[:, 1] < 0]
                downreg_deg_count.append(len(downreg_deg))
                proportions["downcount"] = downreg_deg_count
                proportions[f"DOWN{name}_{tp}"] = downreg_deg

            # Stacked Bar
            stacked1.add_trace(
                go.Bar(x=list_of_days, y=proportions["downcount"], name="Downregulated", marker_color="#636EFA",
                       legendgroup="A"),
                row=stacked_row, col=stacked_col)
            stacked1.add_trace(
                go.Bar(x=list_of_days, y=proportions["upcount"], name="Upregulated", marker_color="#EF553B",
                       legendgroup="B"),
                row=stacked_row, col=stacked_col)
            upreg_deg_count.clear()
            downreg_deg_count.clear()

            stacked_col += 1
            if len(df_dict) % 2 == 0 and stacked_col > 2:
                stacked_col = 1
                stacked_row += 1
            elif math.ceil(len(df_dict) % 3) == 0 and stacked_col > 3:
                stacked_col = 1
                stacked_row += 1

        ## removing duplicate legends
        names = set()
        stacked1.for_each_trace(
            lambda trace:
            trace.update(showlegend=False)
            if (trace.name in names) else names.add(trace.name))

        stacked1.update_layout(showlegend=True, barmode='stack',
                               title=f"Number of DEGs across {tp_or_comp} (based on selected cutoffs)",
                               legend_title_text='DEGs:',
                               font=dict(
                                   family='Arial', size=14))

    stack_success = st.success("Plot complete!")
    time.sleep(0.25)
    stack_success.empty()
    st.plotly_chart(stacked1, use_container_width=True)
    deg_to_dl = [deg_dict[f"{name}_{tp}"] for name in df_names for tp in list_of_days]
    with st.expander("Expand to view and download DEGs per timepoint/comparison"):
        for k, v in deg_dict.items():
            st.write(f"**{k}**", v)
        st.download_button(label="Download DEGs", data=to_excel(deg_to_dl), file_name="DEGs.xlsx")


############################################### Extract DEGs from deg_dict #############################################
def deg_cluster(proportions, log_dfx):
    deglist = [] # first list to add the selected DEGs
    remove_dupes = [] # second list to remove duplicate genes
    temp = [] # third list to add log-filtered datasets to be concatenated
    proportion_keys = list(proportions.keys())
    proportion_keys.remove("upcount")
    proportion_keys.remove("downcount")
    
    select_deg_dicts = postdeg.multiselect("Select DEGs to plot", options=sorted(proportion_keys, key=str.casefold))
    f_width = postdeg.slider("Change clustergram width (in inches)", min_value=5, max_value=20,
                             step=1, value=10)
    f_height = postdeg.slider("Change clustergram height (in inches)", min_value=5, max_value=50,
                              step=1, value=10)
    
    for l in select_deg_dicts:
        degs = proportions[l].index.tolist()
        deglist.append(degs)
    flattened = [val for sublist in deglist for val in sublist]
    for f in flattened:
        if f not in remove_dupes:
            remove_dupes.append(f)

    for k, v in log_dfx.items():
        fc_regex = "log2Fold[-_\s]?Changes?|log2FC"
        log_filter = v.filter(regex=re.compile(fc_regex, flags=re.I))
        log_filter = log_filter.add_prefix(f"{k}_")
        log_filter = log_filter.loc[~log_filter.index.duplicated(keep='first')]
        temp.append(log_filter)
    filter_on = pd.concat(temp, axis=1)

    if len(remove_dupes) > 1:
        specific_cluster = filter_on.loc[remove_dupes]

        # plot
        g = sns.clustermap(specific_cluster, cmap="vlag", method='average', figsize=(f_width, f_height),
                           cbar_pos=(0.01, 0.1, 0.03, 0.1), center=0,
                           col_cluster=True, yticklabels=True)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=11)
        st.pyplot(g)
    else:
        st.warning("Please choose more than 1 DEG")


####################################################### Clustergram #################################################
def clustergram(dfx):
    st.subheader("Gene clustergram")

    dfx_keys = list(dfx.keys())
    # dfx_keys.insert(0, "all")
    fc_regex = "log2Fold[-_\s]?Changes?|log2FC"
    temp = []

    # user input
    with clust_expand:
        select_df = st.multiselect("Select dataframes to use", options=dfx_keys)
        all_df = st.checkbox("All dataframes", value=False)
        gene_list = st.text_area(label="Input list of genes here",
                                 help="Please use one of the following delimiters: line breaks, commas, or semicolons")
        g_width = clust_expand.slider("Change clustergram width (in inches)", min_value=5, max_value=20,
                                      step=1, value=10, key='reg1')
        g_height = clust_expand.slider("Change clustergram height (in inches)", min_value=5, max_value=50,
                                       step=1, value=10, key='reg2')

    if all_df:
        for k, v in dfx.items():
            log_filter = v.filter(regex=re.compile(fc_regex, flags=re.I))
            log_filter = log_filter.add_prefix(f"{k}_")
            log_filter = log_filter.loc[~log_filter.index.duplicated(keep='first')]
            temp.append(log_filter)
        filter_on = pd.concat(temp, axis=1)
    else:
        if len(select_df) != 0:
            for d in select_df:
                used = dfx[d]
                log_filter = used.filter(regex=re.compile(fc_regex, flags=re.I))
                log_filter = log_filter.add_prefix(f"{d}_")
                log_filter = log_filter.loc[~log_filter.index.duplicated(keep='first')]
                temp.append(log_filter)
            filter_on = pd.concat(temp, axis=1)
        else:
            st.stop()

    if len(gene_list) == 1:
        clust_expand.warning(
            "Please enter more than one gene in the list and/or select more than one column to plot the clustergram.")

    elif len(gene_list) > 1:
        genes = gene_list.replace(";", ",").replace(" ", ",").replace("\n", ",").split(',')
        # where the user may use multiple delimiters
        remove_dupes = []
        for g in genes:
            if g not in remove_dupes:
                remove_dupes.append(g)

        gene_final = [x.upper() for x in remove_dupes if x != ""]
        specific_cluster = filter_on.loc[gene_final]

        # clustergram
        g = sns.clustermap(specific_cluster, cmap="vlag",
                           method='average', figsize=(g_width, g_height),
                           cbar_pos=(0.01, 0.1, 0.03, 0.15),
                           center=0, col_cluster=True, yticklabels=True)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=11)
        st.pyplot(g)

# ############################################### Enrichr ##############################################################
degs_but_manual = 0 # If is 0, means user is using DEGs for Enrichr, if is 1, user chooses to add genes manually

def select_enrichr_dataset():
    geneset_dict = {
        "Blood Transcriptomic Modules (BTM)": "BTM.gmt",
        "Reactome": "Reactome.gmt",
        "Vaccinomics (Inhouse)": "Vaccinomics.gmt", "GO Cellular Component 2021": "GO_Cellular_Component_2021",
        "GO Biological Process 2021": "GO_Biological_Process_2021",
        "GO Molecular Function 2021":"GO_Molecular_Function_2021"
    }

    # Selecting genesets (BTM or reactome) to plot from a list
    geneset = enrichr_exp.radio(label='# Select a geneset :', help='Choose BTM or Reactome geneset',
                                options=geneset_dict.keys())
    return geneset_dict[geneset]


def genes_used(premade_dict=None):
    deglist = []
    remove_dupes = []
    temp = []
    if premade_dict is not None:
        choose_genetype = enrichr_exp.selectbox("Select whether to use DEGs or manually add gene list",
                                                options=["DEGs", "Add manually"])
        if choose_genetype == "DEGs":
            proportion_keys = list(premade_dict.keys())
            proportion_keys.remove("upcount")
            proportion_keys.remove("downcount")
            uplist = []
            downlist = []
            upkeys = [x for x in proportion_keys if re.search("UP", x)]
            select_upDEGs = enrichr_exp.multiselect("Select upregulated DEGs to use (optional)", options=upkeys,
                                                    help="Either upregulated or downregulated DEGs must be selected if using the DEGs option")

            for s in select_upDEGs:
                ups = premade_dict[s].index.tolist()
                uplist.append(ups)

            downkeys = [x for x in proportion_keys if re.search("DOWN", x)]
            select_downDEGs = enrichr_exp.multiselect("Select downregulated DEGs to use (optional)", options=downkeys,
                                                      help="Either upregulated or downregulated DEGs must be selected if using the DEGs option")
            for s in select_downDEGs:
                downs = premade_dict[s].index.tolist()
                downlist.append(downs)

            gene_final = [uplist, downlist]


        elif choose_genetype == "Add manually":
            gene_in = enrichr_exp.text_area(label="Input list of at least 3 genes here",
                                            help="Please use one of the following delimiters:line breaks, commas, or semicolons")
            genes = gene_in.replace(";", ",").replace(" ", ",").replace("\n", ",").split(',')
            # where the user may use multiple delimiters
            for g in genes:
                if g not in remove_dupes:
                    remove_dupes.append(g)
            gene_final = [x.upper() for x in remove_dupes if x != ""]
            global degs_but_manual
            degs_but_manual = 1
    else:
        gene_in = enrichr_exp.text_area(label="Input list of at least 3 genes here",
                                        help="Please use one of the following delimiters:line breaks, commas, or semicolons")

        genes = gene_in.replace(";", ",").replace(" ", ",").replace("\n", ",").split(',')
        # where the user may use multiple delimiters
        for g in genes:
            if g not in remove_dupes:
                remove_dupes.append(g)
        gene_final = [x.upper() for x in remove_dupes if x != ""]

    return gene_final


def execute_enrichr(genelist, select_dataset, use_degs=False):
    st.info("Expand the plot to view all of the terms.")
    enrichr_results_exp = st.expander("Expand for enrichr dataframe", expanded=False)
    if not use_degs:
        select_enrichr = gp.enrichr(gene_list=genelist,
                                    gene_sets=select_dataset,
                                    outdir=None,
                                    no_plot=True,
                                    cutoff=0.5  # test dataset, use lower value from range(0,1)
                                    )
        # Sort values by adjusted p-value
        data = select_enrichr.results
        data.set_index("Term", inplace=True)
        data = data.sort_values(by=['Adjusted P-value'])
        # Drop the unimportant variables
        data_trunc = data.drop("Gene_set", axis=1)
        # Filter data by adjusted p-value and rank
        data_sig = data_trunc[(data_trunc['Adjusted P-value'] < 0.05)]
        # Calculate -logP
        data_sig['-logP'] = np.log10(data_sig['Adjusted P-value']) * (-1)
        prettify_geneset = select_dataset.replace(".gmt", "")
        enrichr_results_exp.write(f"**EnrichR analysis using {prettify_geneset}**")
        enrichr_results_exp.write(data_trunc)
        enrichr_download = [data_trunc]
        enrichr_results_exp.download_button(label="Download Enrichr dataframe", data=to_excel(enrichr_download),
                                            file_name="enrichr_analysis.xlsx")

        # Plot bar graph
        toplot = data_sig.sort_values("-logP", ascending=True).tail(10)
        fig = go.Figure(data=go.Bar(x=toplot['-logP'], y=toplot.index, orientation='h', marker_color="#4FC04F"))
        fig.update_layout(title="Enrichr analysis of query genes", title_x=0.5, yaxis={'tickmode':'linear'})
        st.plotly_chart(fig, use_container_width=True)


    else:
        ups = genelist[0]
        downs = genelist[1]
        flattenedup = [val for sublist in ups for val in sublist]
        flatteneddown = [val for sublist in downs for val in sublist]
        remove_dupes_up, remove_dupes_down = [], []
        data_up_trunc, data_down_trunc = None, None
        for f in flattenedup:
            if f not in remove_dupes_up:
                remove_dupes_up.append(f)
        for f in flatteneddown:
            if f not in remove_dupes_down:
                remove_dupes_down.append(f)

        if len(ups) != 0:
            gene_final_up = sorted(remove_dupes_up)
            enrichr_up = gp.enrichr(gene_list=gene_final_up,
                                    gene_sets=select_dataset,
                                    outdir=None,
                                    no_plot=True,
                                    cutoff=0.5  # test dataset, use lower value from range(0,1)
                                    )
            # Sort values by adjusted p-value
            data_up = enrichr_up.results
            data_up.set_index("Term", inplace=True)
            data_up = data_up.sort_values(by=['Adjusted P-value'])
            # Drop the unimportant variables
            data_up_trunc = data_up.drop("Gene_set", axis=1)
            # Filter data by adjusted p-value and rank
            data_up_sig = data_up_trunc[(data_up_trunc['Adjusted P-value'] < 0.05)]
            # Calculate -logP
            data_up_sig['-logP'] = np.log10(data_up_sig['Adjusted P-value']) * (-1)

        if len(downs) != 0:
            gene_final_down = sorted(remove_dupes_down)
            enrichr_down = gp.enrichr(gene_list=gene_final_down,
                                      gene_sets=select_dataset,
                                      outdir=None,
                                      no_plot=True,
                                      cutoff=0.5  # test dataset, use lower value from range(0,1)
                                      )
            # Sort values by adjusted p-value
            data_down = enrichr_down.results
            data_down.set_index("Term", inplace=True)
            data_down = data_down.sort_values(by=['Adjusted P-value'])
            # Drop the unimportant variables
            data_down_trunc = data_down.drop("Gene_set", axis=1)
            # Filter data by adjusted p-value and rank
            data_down_sig = data_down_trunc[(data_down_trunc['Adjusted P-value'] < 0.05)]
            # Calculate -logP
            data_down_sig['-logP'] = np.log10(data_down_sig['Adjusted P-value']) * (-1)

        prettify_geneset = select_dataset.replace(".gmt", "").replace("_", " ")
        enrichr_results_exp.write(f"**EnrichR analysis using {prettify_geneset}**")

        if data_up_trunc is not None and data_down_trunc is not None:
            toplot_up = data_up_sig.sort_values("-logP", ascending=True).tail(10)
            toplot_down = data_down_sig.sort_values("-logP", ascending=True).tail(10)
            with enrichr_results_exp:
                st.write("Enrichr with upregulated DEGs")
                st.dataframe(data_up_trunc)
                st.write("Enrichr with downregulated DEGs")
                st.dataframe(data_down_trunc)
                enrichr_download = [data_up_trunc, data_down_trunc]
                st.download_button(label="Download Enrichr dataframe", data=to_excel(enrichr_download),
                                   file_name="enrichr_updownDEGs_analysis.xlsx")

            fig = make_subplots(rows=2, cols=1, subplot_titles=["Upregulated DEGs", "Downregulated DEGs"],
                                x_title="-logP", shared_xaxes=True)
            fig.add_trace(go.Bar(x=toplot_up['-logP'], y=toplot_up.index,
                                 orientation='h', marker_color="#EF553B"),
                          row=1, col=1)
            fig.add_trace(go.Bar(x=toplot_down['-logP'], y=toplot_down.index,
                                 orientation='h', marker_color="#636EFA"),
                          row=2, col=1)
            fig.update_yaxes(tickmode='linear',tick0=0, dtick=0)


        elif data_up_trunc is not None and data_down_trunc is None:
            toplot_up = data_up_sig.sort_values("-logP", ascending=True).tail(10)
            with enrichr_results_exp:
                st.write("Enrichr with upregulated DEGs")
                st.dataframe(data_up_trunc)
                enrichr_download = [data_up_trunc]
                st.download_button(label="Download Enrichr dataframe", data=to_excel(enrichr_download),
                                   file_name="enrichr_upDEGs_analysis.xlsx")
            fig = go.Figure(go.Bar(x=toplot_up['-logP'], y=toplot_up.index, orientation='h', marker_color="#EF553B"))
            fig.update_xaxes(title="-logP")
            fig.update_yaxes(title="Term")

        else:
            toplot_down = data_down_sig.sort_values("-logP", ascending=True).tail(10)
            with enrichr_results_exp:
                st.write("Enrichr with downregulated DEGs")
                st.dataframe(data_down_trunc)
                enrichr_download = [data_down_trunc]
                st.download_button(label="Download Enrichr dataframe", data=to_excel(enrichr_download),
                                   file_name="enrichr_downDEGs_analysis.xlsx")
            fig = go.Figure(go.Bar(x=toplot_down['-logP'], y=toplot_down.index, orientation='h', marker_color="#636EFA"))
            fig.update_xaxes(title="-logP")
            fig.update_yaxes(title="Term", dtick=0)

        fig.update_layout(title="Enriched Pathways (Top 10), adjpvalue < 0.05", title_x=0.5, showlegend=False, yaxis={'tickmode':'linear'})
        enrichr_results_exp.info("If nothing was plotted in the bar chart, the pathways did not meet the cutoff of adjusted p-value < 0.05")
        st.plotly_chart(fig, use_container_width=True)

##################################### Correlation Matrix (using original df) #########################################
def corr_matrix(dfx):
    st.subheader("Correlation matrix")
    fc_regex = "log2Fold[-_\s]?Changes?|log2FC"

    updated_df_list = []
    for k, v in dfx.items():
        filtered = v.filter(regex=fc_regex)
        df_new = filtered.add_prefix(f"{k}_")
        df_new = df_new.loc[~df_new.index.duplicated(keep='first')]
        updated_df_list.append(df_new)

    concat_fc = pd.concat(updated_df_list, axis=1)
    concat_fc = concat_fc.dropna()
    concat_corr = concat_fc.corr()

    # plot
    mask = np.triu(np.ones_like(concat_corr, dtype=bool))
    df_mask = concat_corr.mask(mask)
    z_raw = df_mask.to_numpy()
    z_text = np.around(z_raw, decimals=3)  # Only show rounded value (full value on hover)
    corr_matrix = ff.create_annotated_heatmap(z=z_raw,
                                              annotation_text=z_text,
                                              x=df_mask.columns.tolist(),
                                              y=df_mask.columns.tolist(),
                                              colorscale=px.colors.diverging.RdBu_r,
                                              zmid=0,
                                              hoverinfo='z',  # Shows hoverinfo for null values
                                              showscale=True, ygap=1, xgap=1
                                              )

    corr_matrix.update_xaxes(side="bottom")
    if is_tp == 1:
        titles = "Time series correlation matrix (based on log2FC values)"
    else:
        titles = "Comparison corrlation matrix (based on log2FC values)"
    corr_matrix.update_layout(
        title_text=titles,
        title_x=0.5,
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
    for i in range(len(corr_matrix.layout.annotations)):
        corr_matrix.layout.annotations[i].font.size = 14
        if corr_matrix.layout.annotations[i].text == 'nan':
            corr_matrix.layout.annotations[i].text = ""

    cmatrix = st.success("Plot complete!")
    time.sleep(0.25)
    cmatrix.empty()
    st.info(
        "If the plot is too small, please hover over the plot and click the expand button on the top right corner of the plot.")
    st.plotly_chart(corr_matrix, use_container_width=True)

##################################### Choose app to build dashboard ##################################################
choose_app = st.sidebar.multiselect("Choose an app to render in the main page 👉",
                                    options=["volcano plot", "DEGs", "clustergram", "enrichr", "correlation matrix"])

for c in choose_app:
    with st.spinner("🔨Building your dashboard 🔨"):
        time.sleep(0.25)
        if c == "volcano plot":
            list_of_days = timepoints(df_dict)
            colorlist = n_colors(list_of_days)
            dfx = check_log(df_dict)
            volcano(dfx, list_of_days, colorlist)

        elif c == "DEGs":
            list_of_days = timepoints(df_dict)
            colorlist = n_colors(list_of_days)
            dfx = check_log(df_dict)
            degs(dfx, list_of_days, colorlist)
            postdeg = st.sidebar.expander("Expand to plot clustergram based on DEGs", expanded=False)
            with postdeg:
                plot_deg_clust = st.checkbox("Plot DEGs in Gene Clustergram", value=False)
            if plot_deg_clust:
                deg_cluster(proportions, dfx)

        elif c == "clustergram":
            dfx = check_log(df_dict)
            clust_expand = st.sidebar.expander("Expand for user-input gene clustergram", expanded=False)
            clustergram(dfx)

        elif c == 'enrichr':
            enrichr_exp = st.sidebar.expander("Expand for Enrichr pathway analysis", expanded=False)
            if "DEGs" in choose_app:
                select_dataset = select_enrichr_dataset()
                genelist = genes_used(premade_dict=proportions)
                run_enrichr = enrichr_exp.button("Run Enrichr")
                if run_enrichr:
                    if degs_but_manual == 0:
                        execute_enrichr(genelist=genelist, select_dataset=select_dataset, use_degs=True)

                    elif degs_but_manual == 1:
                        execute_enrichr(genelist=genelist, select_dataset=select_dataset, use_degs=False)
                else:
                    st.stop()
            else:
                select_dataset = select_enrichr_dataset()
                genelist = genes_used(premade_dict=None)
                if len(genelist) != 0:
                    execute_enrichr(genelist=genelist, select_dataset=select_dataset, use_degs=False)
                else:
                    st.stop()

        elif c == "correlation matrix":
            dfx = check_log(df_dict)
            corr_matrix(dfx)
