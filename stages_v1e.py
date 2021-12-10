#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import time
import math
import re
import base64
from io import BytesIO
import inflect
import dateparser
from datetime import datetime

import gseapy as gp
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import dendrogram, linkage
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


# Update: added internal date-gene conversion
# Update: Improved caching for enrichr and prerank to prevent slow-down of app when no changes are made to enrichr/prerank results
# Update: Fixed bug where volcano plot was unable to be freely manipulated (negative log did not change the graph at (0,0))
# Update v1e: Changes to clustergram to set fold change and have the legend include log2FC

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


def get_table_download_link(df, purpose):  # keeping just in case download button fails
    """Generates a link allowing the data in a given panda dataframe to be downloaded
    in:  dataframe
    out: href string
    """
    val = to_excel(df)
    b64 = base64.b64encode(val)
    return f'<a href="data:application/octet-stream;base64,{b64.decode()}" download="{purpose}.xlsx">' \
           f'📥 Download {purpose} as Excel file 📥</a>'  # decode b'abc' => abc


st.title("STAGEs Dashboard \U0001F4CA")

################################################# Documentation ########################################################
documentation = st.sidebar.checkbox("Read the Docs", value=False, key='documentation')
################################################# File Uploader ########################################################
df_query = st.sidebar.file_uploader(
    'Upload your .csv/.xlsx file. A demo dataset will be uploaded if no files are uploaded',
    accept_multiple_files=True)

df_dict = {} # initial, but should use cleaned_dict after uploading and qc checks
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
            selected_sheet = st.sidebar.multiselect(label="* Select which sheet to read in", options=x.keys())
            for i in selected_sheet:
                data = x[i]
                df_dict[i] = data
                df_names.append(i)

else:
    testdata = st.experimental_memo(pd.read_csv)("demo_dataframe_corrected.csv", index_col=0)
    testname = "Demo"
    df_dict[testname] = testdata
    df_names.append(testname)

for df in df_dict.values():
    df.index = df.index.astype(str, copy=False) # expand to format actual dates from excel sheets as text


####### Important ########
deg_dict = {}  ###########
proportions = {}  ########
##########################

# gene_symbols = pd.read_csv("/Users/clara/Desktop/Actual Work/Correcting Date Genes/gene_date.csv") # local
# gene_symbols = pd.read_csv("gene_date.csv")
# old_symbols = gene_symbols.iloc[:, 0].tolist()
# new_symbols = gene_symbols.iloc[:, 3].tolist()

@st.experimental_memo
def clean_ref():
    for_ref = pd.read_csv("/Users/clara/Dropbox/Streamlit_app/Date Gene Converter/hgnc-symbol-check.csv") # local
    # for_ref = pd.read_csv("hgnc-symbol-check2.csv") # github
    for_ref.reset_index(drop=True,inplace=True)
    for_ref.columns = for_ref.iloc[0,:]
    for_ref.drop(index=0, inplace=True)
    for_ref.drop(columns="Match type", inplace=True)
    for_ref.rename(columns={"Input":"Previous Symbol"}, inplace=True)
    return for_ref
reference_symbols = clean_ref()

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

    Demo examples are provided. You can try out the demo examples to familiarise yourself with the apps before uploading your dataset.

    ## Data safety and security
    The data you upload is safe and is never stored anywhere.

    ## Contributors
    These apps are jointly made by myself (Kuan Rong Chan), Clara Koh, Justin Ooi and Gabrielle Lee from Duke-NUS, Department of Emerging Infectious Diseases. I am also thankful for Eugenia Ong and Ayesa Syenina from VIREMICS for their constructive feedback. These apps are now free for everyone to use, but for a limited period of time as we are constantly upgrading the apps. For more details on what we do, feel free to visit us at [omicsdiary.com](https://omicsdiary.com/).
        ''')


if documentation:
    read_docs()

################################################# Date-Gene Conversion ################################################
#######################################################################################################################
cleaned_dict = {}
def qc_df(dfs):
    with st.spinner("Checking uploaded dataframe for dates in gene column"):
        time.sleep(0.5)
    p = inflect.engine()
    corrected = {"Dec-01_1st": "DELEC1", "01-Dec_1st": "DELEC1", "Mar-03_1st": "MARCHF3", "03-Mar_1st": "MARCHF3",
                 "Mar-04_1st": "MARCHF4", "04-Mar_1st": "MARCHF4", "Mar-05_1st": "MARCHF5", "05-Mar_1st": "MARCHF5",
                 "Mar-06_1st": "MARCHF6", "06-Mar_1st": "MARCHF6", "Mar-07_1st": "MARCHF7", "07-Mar_1st": "MARCHF7",
                 "Mar-08_1st": "MARCHF8", "08-Mar_1st": "MARCHF8", "Mar-09_1st": "MARCHF9", "09-Mar_1st": "MARCHF9",
                 "Mar-10_1st": "MARCHF10", "10-Mar_1st": "MARCHF10", "Mar-11_1st": "MARCHF11", "11-Mar_1st": "MARCHF11",
                 "Sep-15_1st": "SELENOF", "15_Sep_1st": "SELENOF", "Sep-01_1st": "SEPTIN1", "01-Sep_1st": "SEPTIN1",
                 "Sep-02_1st": "SEPTIN2", "02-Sep_1st": "SEPTIN2", "Sep-03_1st": "SEPTIN3", "03-Sep_1st": "SEPTIN3",
                 "Sep-04_1st": "SEPTIN4", "04-Sep_1st": "SEPTIN4", "Sep-05_1st": "SEPTIN5", "05-Sep_1st": "SEPTIN5",
                 "Sep-06_1st": "SEPTIN6", "06-Sep_1st": "SEPTIN6", "Sep-07_1st": "SEPTIN7", "07-Sep_1st": "SEPTIN7",
                 "Sep-08_1st": "SEPTIN8", "08-Sep_1st": "SEPTIN8", "Sep-09_1st": "SEPTIN9", "09-Sep_1st": "SEPTIN9",
                 "Sep-10_1st": "SEPTIN10", "10-Sep_1st": "SEPTIN10", "Sep-11_1st": "SEPTIN11", "11-Sep_1st": "SEPTIN11",
                 "Sep-12_1st": "SEPTIN12", "12-Sep_1st": "SEPTIN12", "Sep-14_1st": "SEPTIN14", "14-Sep_1st": "SEPTIN14"
                 }
    ################ Contains dates and March-01/March-02 and have to be resolved ####################
    def march_resolver(df):
        find = [g for g in df.index.tolist() if re.search("Mar|Apr|Sept?|Oct|Dec", g)]
        formatted = {}
        for d in find:
            zero_pad = re.search("[0-9]{2}", d)
            num = re.findall("[0-9]*", d)
            og_num = [x for x in num if x != ""]
            month = re.findall("[A-Za-z]*", d)
            og_month = [x for x in month if x != ""]
            if not zero_pad:
                a = f"{og_month[0]}-0{og_num[0]}" # still can't use dateparser bc time fmt requires zero-padded number
                formatted[d] = a
            else:
                a = f"{og_month[0]}-{og_num[0]}"
                formatted[d] = a

        each_df_exp = st.expander(f"Expand to resolve naming issues for {k} dataframe", expanded=False)
        found = df.loc[find]
        found.rename(index=formatted, inplace=True)
        found = found.drop_duplicates() # ensures that there aren't duplicate rows (not just duplicate row names)
        found.reset_index(drop=False, inplace=True)
        index_name = found.columns.tolist()[0]
        found[index_name] += found.groupby(index_name).cumcount().add(1).map(p.ordinal).radd('_')
        found.set_index(index_name, inplace=True)
        mar1 = [f for f in found.index.tolist() if re.search("Mar-0?1_1st|0?1-Mar_1st|Mar-0?1_2nd|Mar-0?1_2nd", f)]
        if len(mar1) !=0:
            mar1_df = found.loc[mar1]

            with each_df_exp:
                st.write(f"**MAR01 Genes: {k} Dataframe**")
                st.info("Genes like MARCH1 and MARC1 have to be differentiated by function as they are both corrected to Mar-01 in Excel."
                                 " Check HGNC symbol reference in the sidebar for reference. 👈")
                st.dataframe(mar1_df.astype(str))

                first_mar01_fx = st.selectbox(f"Select the name and function that {mar1[0]} corresponds to for {k} dataframe",
                                              options=["MTARC1: mitochondrial amidoxime reducing component 1",
                                                      "MARCHF1: membrane associated ring-CH-type finger 1"])
            # function below can still apply to genes with only 1 MAR-01 gene because the dictionary will only match those found in the data
            if first_mar01_fx == "MTARC1: mitochondrial amidoxime reducing component 1":
                first_mar01 = "MTARC1"
                second_mar01 = "MARCHF1"
            elif first_mar01_fx == "MARCHF1: membrane associated ring-CH-type finger 1":
                first_mar01 = "MARCHF1"
                second_mar01 = "MTARC1"

            corrected["Mar-01_1st"] = first_mar01
            corrected["01-Mar_1st"] = first_mar01
            corrected["Mar-01_2nd"] = second_mar01
            corrected["01-Mar_2nd"] = second_mar01

        mar2 = [f for f in found.index.tolist() if re.search("Mar-0?2_1st|0?2-Mar_1st|Mar-0?2_2nd|02-Mar_2nd", f)]
        if len(mar2) !=0:
            mar2_df = found.loc[mar2]

            each_df_exp.write(f"**MAR02 Genes: {k} Dataframe**")
            each_df_exp.info(
                "Genes like MARCH2 and MARC2 have to be differentiated by function as they are both corrected to Mar-01 in Excel."
                " Check old and new HGNC symbols in the sidebar for reference. 👈")
            each_df_exp.dataframe(mar2_df.astype(str))

            first_mar02_fx = each_df_exp.selectbox(f"Select the name and function that {mar2[0]} corresponds to for {k} dataframe",
                                          options=[
                                              "MTARC2: mitochondrial amidoxime reducing component 2",
                                              "MARCHF2: membrane associated ring-CH-type finger 2"])

            if first_mar02_fx == "MTARC2: mitochondrial amidoxime reducing component 2":
                first_mar02 = "MTARC2"
                second_mar02 = "MARCHF2"
            elif first_mar02_fx == "MARCHF2: membrane associated ring-CH-type finger 2":
                first_mar02 = "MARCHF2"
                second_mar02 = "MTARC2"

            corrected["Mar-01_1st"] = first_mar01
            corrected["01-Mar_1st"] = first_mar01
            corrected["Mar-01_2nd"] = second_mar01
            corrected["01-Mar_2nd"] = second_mar01
            corrected["Mar-02_1st"] = first_mar02
            corrected["02-Mar_1st"] = first_mar02
            corrected["Mar-02_2nd"] = second_mar02
            corrected["02-Mar_2nd"] = second_mar02

        else:
            corrected["Mar-01_1st"] = first_mar01
            corrected["01-Mar_1st"] = first_mar01
            corrected["Mar-01_2nd"] = second_mar01
            corrected["01-Mar_2nd"] = second_mar01

        found["Gene"] = pd.Series(corrected)  # in order to rename just change this to found.rename(corrected)
        found.reset_index(drop=True, inplace=True)  # remove the gene index with the dates
        found.rename(columns={"Gene": "gene"}, inplace=True)  # rename the incoming column to be used as index
        found.set_index('gene', inplace=True)  # set the index of these date genes using corrected names
        df = df.drop(index=find)  # drop the date genes from the main df
        df2 = pd.concat([df, found], axis=0)  # join these genes back to the main df
        df2.sort_index(axis=0, ascending=True, inplace=True)  # sort alphabetically
        df2.reset_index(drop=False, inplace=True)
        df2.rename(columns={'index': index_name}, inplace=True)
        df2.set_index(index_name, inplace=True)
        cleaned_dict[k] = df2
        return

    ############ Contains dates but no march-01/march-02 and thus nothing to resolve ##############
    def date_resolver(df, date_search):
        formatted = {}
        for d in date_search:
            zero_pad = re.search("[0-9]{2}", d)
            num = re.findall("[0-9]*", d)
            og_num = [x for x in num if x != ""]
            month = re.findall("[A-Za-z]*", d)
            og_month = [x for x in month if x != ""]
            if not zero_pad:
                a = f"{og_month[0]}-0{og_num[0]}" # still can't use dateparser as python time fmts only read zero-padded no.
                formatted[d] = a
            else:
                a = f"{og_month[0]}-{og_num[0]}"
                formatted[d] = a
        found = df.loc[date_search]
        found.rename(index=formatted, inplace=True)
        found = found.drop_duplicates()  # ensures that there aren't duplicate rows (not just duplicate row names)
        found.reset_index(drop=False, inplace=True)
        index_name = found.columns.tolist()[0]
        found[index_name] += found.groupby(index_name).cumcount().add(1).map(p.ordinal).radd('_')
        found.set_index(index_name, inplace=True)
        found.rename(index=corrected, inplace=True)  # in order to rename just change this to found.rename(corrected)
        df = df.drop(index=date_search)  # drop the date genes from the main df
        df2 = pd.concat([df, found], axis=0)  # join these genes back to the main df
        df2.sort_index(axis=0, ascending=True, inplace=True)  # sort alphabetically
        cleaned_dict[k] = df2
        return

    ############################## Dates are only numbers ##########################################
    def numeric_date(k, df, numdate):
        num_exp = st.expander(f"Expand if {k}'s date format is numerical (eg. yyyy/mm/dd)")
        date_fmt = num_exp.radio(f"Select the format that {k} dataframe is in",
                                 options=["yyyy-dd-mm", "yyyy-mm-dd", "dd-mm-yyyy", "mm-dd-yyyy"])
        info_stored = num_exp.radio(f"Select how {k}'s dates should be read to derive gene names (Hover '?' for help)",
                                    options=['month-year', 'month-day'],
                                    help='For example, 2001-03-09 (yyyy-mm-dd) may either be Mar-01 (MARCHF1) or Mar-09 (MARCHF9).',
                                    index=1)
        num_exp.info("If you're unsure about the above option, check your dataframe and select 'month-year'. "
                     "We recommend you to check the converted dataframe to ensure that the dates are converted correctly. If unsuccessful, <NA> symbols will populate at the bottom of the dataframe.")
        if info_stored == 'month-year':
            strfmt = '%b-%y'
        else:
            strfmt = '%b-%d'
        a = []
        df.reset_index(drop=False, inplace=True)
        index_nm = df.columns.tolist()[0]
        for n in numdate:
            indices = df.index[df[index_nm] == n].tolist()
            if indices not in a:
                a.append(indices)
        flattened = [item for sublist in a for item in sublist]
        found = df.iloc[flattened, :]
        found.set_index(index_nm, inplace=True)
        num_exp.write(f"**{k} dataframe**")
        num_exp.dataframe(found)

        if date_fmt == "yyyy-dd-mm":
            extracted = {n: (dateparser.parse(n, date_formats=["%Y-%d-%m"])).strftime(strfmt) for n in numdate}
        elif date_fmt == "yyyy-mm-dd":
            extracted = {n: (dateparser.parse(n, date_formats=["%Y-%m-%d"])).strftime(strfmt) for n in numdate}
        elif date_fmt == "dd-mm-yyyy":
            extracted = {n: (dateparser.parse(n, date_formats=["%d-%m-%Y"])).strftime(strfmt) for n in numdate}
        elif date_fmt == "mm-dd-yyyy":
            extracted = {n: (dateparser.parse(n, date_formats=["%m-%d-%Y"])).strftime(strfmt) for n in numdate}

        df.set_index(index_nm, inplace=True)
        df.rename(index=extracted, inplace=True)
        return df

    ############################ Just old symbols and no date issues ###############################
    def nodates(df):
        corrected = {}
        for i in range(len(reference_symbols)):
            key = reference_symbols.iloc[i, 0]
            value = reference_symbols.iloc[i, 1]
            corrected[key] = value

        df.rename(index=corrected, inplace=True)
        cleaned_dict[k] = df
        return
################################################# Code Flow ##########################################################

# since it's quite likely that old symbols don't exist together with dates (because once opened in excel all are dates),
# this is an all-or-none approach where if date search picks up sth, old search will be empty
# if both lists are empty, nothing is wrong, then cleaned dict[k] = df_dict[k], where k is df with no error

    ismar, isnums = 0, 0

    for k,df in dfs.items():
        date_search = [g for g in df.index.tolist() if re.search("(Mar|Apr|Sept?|Oct|Dec)", g)] # dates
        old_symbols = list(reference_symbols['Previous Symbol'])
        old_search = list(set(df.index.tolist()).intersection(set(old_symbols))) # easy way to find old symbols in df index
        if len(date_search) != 0:
            march_search = [m for m in date_search if re.search("^Mar-0?1|0?1-Mar|Mar-0?2|0?2-Mar", m)] # only march genes
            if len(march_search) != 0:
                ismar += 1
                if ismar == 1:
                    st.subheader("Resolve Duplicate Gene Symbols")
                march_resolver(df) # will have all the date genes (like it is now)
            else:
                date_resolver(df, date_search) # will exclude march

        elif len(old_search) != 0:
            nodates(df) # converts old to new (eg. DEC1 -> DELEC1)

        elif len(date_search) == 0 and len(old_search) == 0:
            numdate = [g for g in df.index.tolist() if re.search("^\d*[-/]?\W", g)]
            if len(numdate) != 0:
                isnums += 1
                if isnums == 1:
                    st.subheader("Resolve Date Format")
                renamed = numeric_date(k,df,numdate)
                march_search = [m for m in renamed.index.tolist() if re.search("^Mar-0?1|0?1-Mar|Mar-0?2|0?2-Mar", m)]  # only march genes
                generic_date = [g for g in renamed.index.tolist() if re.search("(Mar|Apr|Sept?|Oct|Dec)", g)]
                if len(march_search) != 0:
                    ismar += 1
                    if ismar == 1:
                        st.subheader("Resolve Duplicate Gene Symbols")
                        march_resolver(renamed)  # will have all the date genes (like it is now)
                else:
                    date_resolver(renamed, generic_date)  # will exclude march
            else:
                noerror = st.success(f"No errors detected for {k} dataframe")
                time.sleep(0.25)
                noerror.empty()
                cleaned_dict[k] = df
        else:
            noerror = st.success(f"No errors detected for {k} dataframe")
            time.sleep(0.25)
            noerror.empty()
            cleaned_dict[k] = df
    return

################################################## Grabs Timepoint ####################################################
# Is a variable to allow the following functions to render either timepoint or comparison in plot titles
is_tp = 1


def timepoints(dfs):
    pattern_days = "[-_\s]{1}(min[0-9]*[.]?[0-9]+)|([0-9]*[.]?[0-9]+min)[-_\s]?|(hr[0-9]*[.]?[0-9]+)|([0-9]*[.]?[0-9]+hr?)|(D[ay]*[0-9]*[.]?[0-9]+)|([0-9]*[.]?[0-9]+D[ay]*)"
    tpcomps = []
    for k, df in dfs.items():
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
def check_log(dfs):
    df_dict2 = {}
    for k, df in dfs.items():
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
def volcano(dfs, list_of_days, colorlist):
    st.subheader("Volcano plot analysis")
    vol_expand = st.sidebar.expander("Expand for volcano plot", expanded=False)
    reset = vol_expand.checkbox("Reset to default settings", value=False)
    xaxes = vol_expand.slider("Choose log2 fold-change boundaries for volcano plot",
                              help="The app will plot the values between the user-set range",
                              min_value=-5.0, max_value=5.0, step=0.1, value=(0.0, 0.0))
    if reset:
        xaxes = (0.0, 0.0)
    yaxes = vol_expand.slider("Choose negative log10 p-value boundaries for volcano plot",
                              help="The app will plot the values greater than the user-set value",
                              min_value=0.0, max_value=5.0, step=0.1, value=0.0)
    interactive_volcano = vol_expand.checkbox(label="Show interactive volcano plot", value=False,
                                              help="Facilitates gene name display on hover. This may cause lag")

    if is_tp == 1:
        tp_or_comp = "time-points (vs baseline)"
    else:
        tp_or_comp = "comparisons"

    top10annotation, bottom10annotation = [], []

    if len(dfs) == 1:
        volcano1 = go.Figure()
        fig = plt.figure()
        for k, df in dfs.items():
            for tp, col in zip(list_of_days, colorlist):
                complabels = tp.replace("_", " ").replace("-", " ")
                #### selecting the required FC and pval for plotting
                FC_col_name = list([col for col in df.columns if tp in col if "log2FC" in col])
                fold_changes = df[FC_col_name[0]]
                pval_col_name = list([col for col in df.columns if tp in col if "negative_log_pval" in col])
                # edited to negative instead of neg so we can pass all the dfs through data formatter
                # (to include the log2ratio and -log pval)
                pvals = df[pval_col_name[0]]
                for_annotation = pd.concat([fold_changes, pvals], axis=1)
                if xaxes != (0.0, 0.0) and yaxes != (0.0):
                    user_filter = for_annotation[(for_annotation[pval_col_name[0]] >= yaxes) &
                                                 (for_annotation[FC_col_name[0]].between(xaxes[0], xaxes[1],
                                                                                         inclusive='both')
                                                  )]  # bracket over fc filters
                elif xaxes == (0.0, 0.0) and yaxes != (0.0):
                    user_filter = for_annotation[(for_annotation[pval_col_name[0]] >= yaxes)]

                elif xaxes != (0.0, 0.0) and yaxes == (0.0):
                    user_filter = for_annotation[(for_annotation[FC_col_name[0]].between(xaxes[0], xaxes[1],
                                                                                         inclusive='both')
                                                  )]  # bracket over fc filters
                else:
                    user_filter = for_annotation

                top_10 = user_filter.sort_values(by=FC_col_name[0], ascending=False).head(10)
                bottom_10 = user_filter.sort_values(by=FC_col_name[0], ascending=True).head(10)
                bottom10annotation.append(
                    bottom_10.rename(columns={FC_col_name[0]: "log2FC", pval_col_name[0]: "negative_log_pval"}))
                top10annotation.append(
                    top_10.rename(columns={FC_col_name[0]: "log2FC", pval_col_name[0]: "negative_log_pval"}))

                plt.grid(b=True, which="major", axis="both", alpha=0.3)
                plt.scatter(user_filter[FC_col_name[0]], user_filter[pval_col_name[0]], alpha=0.7, label=complabels)
                plt.title(f"Volcano plot across {tp_or_comp}", loc='center')
                plt.legend(bbox_to_anchor=(1.04, 1), loc='upper left')
                plt.xlabel('log2(Fold-change)')
                plt.ylabel('-log10(p-value)')
                plt.axhline(y=0, color='r', linestyle='dashed')
                plt.axvline(x=0, linestyle='dashed')

                if interactive_volcano:
                    volcano1.add_trace(go.Scatter(x=user_filter[FC_col_name[0]], y=user_filter[pval_col_name[0]],
                                                  mode='markers',
                                                  name=complabels, hovertext=list(user_filter.index),
                                                  line=dict(color=col)
                                                  )
                                       )
            annotationconcat_top = pd.concat(top10annotation, axis=0)
            annotationconcat_top = annotationconcat_top.sort_values(by=["log2FC"], ascending=False).head(10)

            annotationconcat_bottom = pd.concat(bottom10annotation, axis=0)
            annotationconcat_bottom = annotationconcat_bottom.sort_values(by=["log2FC"], ascending=True).head(10)

            topoverall = list(annotationconcat_top.index)
            bottomoverall = list(annotationconcat_bottom.index)

            for i in range(len(annotationconcat_top)):
                plt.annotate(text=topoverall[i], xy=(annotationconcat_top.iloc[i, 0], annotationconcat_top.iloc[i, 1]),
                             xytext=(0, 3), horizontalalignment='center', textcoords='offset points',
                             fontsize=7)  # top 10
            for i in range(len(annotationconcat_bottom)):
                plt.annotate(text=bottomoverall[i],
                             xy=(annotationconcat_bottom.iloc[i, 0], annotationconcat_bottom.iloc[i, 1]),
                             xytext=(0, 3), horizontalalignment='center', textcoords='offset points',
                             fontsize=7)  # bottom 10

        volcano1.update_layout(showlegend=True,
                               title=f"Interactive volcano across {tp_or_comp}",
                               legend_title_text="Timepoint",
                               font=dict(family='Arial', size=14),
                               xaxis_title="log2(Fold-change)", yaxis_title="-log10(p-value)")
    else:
        if len(dfs) % 2 == 0:
            nrows = math.ceil(len(dfs) / 2)
            volcano1 = make_subplots(rows=nrows, cols=2, subplot_titles=(list(dfs.keys())),
                                     x_title="log2(Fold-Change)", y_title="-log10(p-value)")
            v_row, v_col = 1, 1
            j = 1
            fig, axes = plt.subplots(nrows=nrows, ncols=2, sharex=True, sharey=True)
        elif math.ceil(len(dfs) % 3) == 0:
            nrows = math.ceil(len(dfs) / 2)
            volcano1 = make_subplots(rows=nrows, cols=2, subplot_titles=(list(dfs.keys())),
                                     x_title="log2(Fold-Change)", y_title="-log10(p-value)")
            v_row, v_col = 1, 1
            j = 1
            fig, axes = plt.subplots(nrows=nrows, ncols=3, sharex=True, sharey=True)

        for k, df, in dfs.items():
            for tp, clr in zip(list_of_days, colorlist):
                complabels = tp.replace("_", " ").replace("-", " ")
                #### selecting the required FC and pval for plotting
                FC_col_name = list([col for col in df.columns if tp in col if "log2FC" in col])
                fold_changes = df[FC_col_name[0]]
                pval_col_name = list([col for col in df.columns if tp in col if "negative_log_pval" in col])
                # edited to negative instead of neg so we can pass all the dfs through data formatter
                # (to include the log2ratio and -log pval)
                pvals = df[pval_col_name[0]]

                for_annotation = pd.concat([fold_changes, pvals], axis=1)
                if xaxes != (0.0, 0.0) and yaxes != (0.0):
                    user_filter = for_annotation[(for_annotation[pval_col_name[0]] >= yaxes) &
                                                 (for_annotation[FC_col_name[0]].between(xaxes[0], xaxes[1],
                                                                                         inclusive='both')
                                                  )]  # bracket over fc filters
                elif xaxes == (0.0, 0.0) and yaxes != (0.0):
                    user_filter = for_annotation[(for_annotation[pval_col_name[0]] >= yaxes)]

                elif xaxes != (0.0, 0.0) and yaxes == (0.0):
                    user_filter = for_annotation[(for_annotation[FC_col_name[0]].between(xaxes[0], xaxes[1],
                                                                                         inclusive='both')
                                                  )]  # bracket over fc filters
                else:
                    user_filter = for_annotation

                top_10 = user_filter.sort_values(by=FC_col_name[0], ascending=False).head(10)
                bottom_10 = user_filter.sort_values(by=FC_col_name[0], ascending=True).head(10)
                bottom10annotation.append(
                    bottom_10.rename(columns={FC_col_name[0]: "log2FC", pval_col_name[0]: "negative_log_pval"}))
                top10annotation.append(
                    top_10.rename(columns={FC_col_name[0]: "log2FC", pval_col_name[0]: "negative_log_pval"}))

                if len(dfs) % 2 == 0:
                    plt.subplot(nrows, 2, j)
                elif len(dfs) % 3 == 0:
                    plt.subplot(nrows, 3, j)
                plt.grid(b=True, which="major", axis="both", alpha=0.3)
                plt.scatter(user_filter[FC_col_name[0]], user_filter[pval_col_name[0]], alpha=0.7, label=complabels)
                plt.axhline(y=0, color='r', linestyle='dashed')
                plt.axvline(x=0, linestyle='dashed')

                if interactive_volcano:
                    volcano1.add_trace(go.Scatter(x=user_filter[FC_col_name[0]], y=user_filter[pval_col_name[0]],
                                                  mode='markers',
                                                  name=complabels, hovertext=list(df.index),
                                                  line=dict(color=clr)
                                                  ),
                                       row=v_row, col=v_col
                                       )
            annotationconcat_top = pd.concat(top10annotation, axis=0)
            annotationconcat_top = annotationconcat_top.sort_values(by=["log2FC"], ascending=False).head(10)

            annotationconcat_bottom = pd.concat(bottom10annotation, axis=0)
            annotationconcat_bottom = annotationconcat_bottom.sort_values(by=["log2FC"], ascending=True).head(10)

            topoverall = list(annotationconcat_top.index)
            bottomoverall = list(annotationconcat_bottom.index)

            for i in range(len(annotationconcat_top)):
                plt.annotate(text=topoverall[i], xy=(annotationconcat_top.iloc[i, 0], annotationconcat_top.iloc[i, 1]),
                             xytext=(0, 3), horizontalalignment='center', textcoords='offset points',
                             fontsize=7)  # top 10
            for i in range(len(annotationconcat_bottom)):
                plt.annotate(text=bottomoverall[i],
                             xy=(annotationconcat_bottom.iloc[i, 0], annotationconcat_bottom.iloc[i, 1]),
                             xytext=(0, 3), horizontalalignment='center', textcoords='offset points',
                             fontsize=7)  # bottom 10

            top10annotation.clear()
            bottom10annotation.clear()

            j += 1
            v_col += 1
            if (len(dfs) % 2 == 0) and j > 2:
                j = 1
            if (len(dfs) % 3 == 0) and j > 3:
                j = 1

            if (len(dfs) % 2 == 0) and v_col > 2:
                v_col = 1
                v_row += 1
            if (len(dfs) % 3 == 0) and v_col > 3:
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
def degs(dfs, list_of_days, colorlist):
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
        df_list = [dfs[x] for x in df_names]
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

    if len(dfs) == 1:
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
        if len(dfs) % 2 == 0:
            nrows = math.ceil(len(dfs) / 2)
            stacked1 = make_subplots(rows=nrows, cols=2, subplot_titles=(list(dfs.keys())),
                                     x_title=tp_or_comp.title(), y_title='Number of DEGs', shared_yaxes=True)
        elif math.ceil(len(dfs)) % 3 == 0:
            nrows = math.ceil(len(dfs) / 3)
            stacked1 = make_subplots(rows=row_no, cols=3, subplot_titles=(list(dfs.keys())),
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
                proportions[f"DOWN_{name}_{tp}"] = downreg_deg

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
            if len(dfs) % 2 == 0 and stacked_col > 2:
                stacked_col = 1
                stacked_row += 1
            elif math.ceil(len(dfs) % 3) == 0 and stacked_col > 3:
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
        # st.download_button(label="Download DEGs", data=to_excel(deg_to_dl), file_name="DEGs.xlsx")
        st.markdown(get_table_download_link(deg_to_dl, "DEGs"), unsafe_allow_html=True)


############################################### Extract DEGs from deg_dict #############################################
def deg_cluster(proportions, log_dfx):
    st.subheader("Pathway Clustergram from DEGs")
    deglist = []  # first list to add the selected DEGs
    remove_dupes = []  # second list to remove duplicate genes
    temp = []  # third list to add log-filtered datasets to be concatenated
    proportion_keys = list(proportions.keys())
    proportion_keys.remove("upcount")
    proportion_keys.remove("downcount")

    select_deg_dicts = postdeg.multiselect("Select DEGs to plot", options=sorted(proportion_keys, key=str.casefold))
    resetter = postdeg.checkbox("Default settings", help="Do not filter by log2 fold-change cutoff", value=True, key='degbased')
    fc_slider = postdeg.slider("Adjust log2 fold-change here", help="The app will plot the values between the user-set range",
                              min_value=-5.0, max_value=5.0, step=0.1, value=(-1.0,1.0), key='degbased')

    f_width = postdeg.slider("Change clustergram width (in inches)", min_value=5, max_value=20,
                             step=1, value=10)
    f_height = postdeg.slider("Change clustergram height (in inches)", min_value=5, max_value=50,
                              step=1, value=10)
    
    for l in select_deg_dicts:
        if not resetter:
            fc_filter = proportions[l][(proportions[l].iloc[:,1].between(fc_slider[0], fc_slider[1],inclusive='both'))]
        else:
            fc_filter = proportions[l]

        degs = fc_filter.index.tolist()
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
        g = sns.clustermap(specific_cluster, cmap='vlag', method='average', figsize=(f_width, f_height),
                           cbar_pos=(0.01, 0.1, 0.03, 0.1), center=0,
                           col_cluster=True, yticklabels=True,
                           cbar_kws={'label': 'log2FC'})
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=11)
        # with st.expander("Expand for pathway clustergram dataframe from DEGs", expanded=False):
        #     st.write("**User-Input Pathway Clustergram**")
        #     st.dataframe(specific_cluster)
        #     st.markdown(get_table_download_link([specific_cluster], "pathway_clustergram"), unsafe_allow_html=True)
        st.pyplot(g)
        # col_linkage = g.dendrogram_col.linkage
        # st.write("Column Linkage", col_linkage)
        # denfig = plt.figure()
        # dendro_col = dendrogram(col_linkage, orientation='top')
        # st.pyplot(denfig)
        # st.write("Row Linkage", g.dendrogram_row.linkage)
        # denrowfig = plt.figure()
        # dendro_row = dendrogram(g.dendrogram_row.linkage, orientation='left')
        # st.pyplot(denrowfig)

    else:
        st.warning("Please choose more than 1 DEG or increase the log2 fold-change limit.")


####################################################### Clustergram #################################################
def clustergram(dfx):
    st.subheader("Pathway clustergram")
    dfx_keys = list(dfx.keys())
    fc_regex = "log2Fold[-_\s]?Changes?|log2FC"
    temp = []

    # user input
    with clust_expand:
        select_df = st.multiselect("Select dataframes to use", options=dfx_keys)
        all_df = st.checkbox("All dataframes", value=False)
        gene_list = st.text_area(label="Input list of genes here",
                                 help="Please use one of the following delimiters: line breaks, commas, or semicolons")
        resetter = st.checkbox("Default settings", value=True, help="Do not filter by log2 fold-change cutoff", key='userclust')
        fc_slider = st.slider("Adjust log2 fold-change here", help="The app will plot the values between the user-set range",
                                    min_value=-5.0, max_value=5.0, step=0.1, value=(-1.0,1.0), key='userclust')
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
        gene_cluster = filter_on.loc[gene_final]
        if not resetter:
            specific_cluster = gene_cluster[gene_cluster.iloc[:,1].between(fc_slider[0], fc_slider[1],inclusive='both')]
        else:
            specific_cluster = gene_cluster
        # with st.expander("Expand for pathway clustergram dataframe", expanded=False):
        #     st.write("**User-Input Pathway Clustergram**")
        #     st.dataframe(specific_cluster)
        #     st.markdown(get_table_download_link([specific_cluster], "pathway_clustergram"), unsafe_allow_html=True)

        try:
            # clustergram
            g = sns.clustermap(specific_cluster, cmap="vlag",
                            method='average', figsize=(g_width, g_height),
                            cbar_pos=(0.01, 0.1, 0.03, 0.15),
                            center=0, col_cluster=True, yticklabels=True, cbar_kws={'label':'log2FC'})
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=11)
            st.pyplot(g)
            # st.write("Column Linkage", g.dendrogram_col.linkage)
            # denfig = plt.figure()
            # dendro_col = dendrogram(g.dendrogram_col.linkage)
            # st.pyplot(denfig)
            # st.write("Row Linkage", g.dendrogram_row.linkage)
        except ValueError:
            st.warning("Increase the log2FC limit or add more genes.")


################################################# Enrichr ##############################################################
degs_but_manual = 0  # If is 0, means user is using DEGs for Enrichr, if is 1, user chooses to add genes manually

def select_enrichr_dataset():
    geneset_dict = {
        "Blood Transcriptomic Modules (BTM)": "BTM.gmt",
        "Reactome 2021": "Reactome.gmt",
        "Vaccinomics (In-house)": "Vaccinomics.gmt", "GO Cellular Component 2021": "GO_Cellular_Component_2021",
        "GO Biological Process 2021": "GO_Biological_Process_2021",
        "GO Molecular Function 2021": "GO_Molecular_Function_2021",
        "KEGG 2021 Human": "KEGG_2021_Human"
    }

    # Selecting genesets (BTM or reactome) to plot from a list
    geneset = enrichr_exp.radio(label='Select a geneset', options=geneset_dict.keys(), key='enrichr')
    return geneset_dict[geneset]


def genes_used(premade_dict=None):
    deglist, remove_dupes, temp, up_enrichr, down_enrichr = [], [], [], [], []
    if premade_dict is not None:
        choose_genetype = enrichr_exp.selectbox("Select whether to use DEGs or manually add gene list",
                                                options=["DEGs", "Add manually"])
        if choose_genetype == "DEGs":
            proportion_keys = list(premade_dict.keys())
            proportion_keys.remove("upcount")
            proportion_keys.remove("downcount")
            uplist = []
            downlist = []

            select_DEGs = enrichr_exp.multiselect("Select DEGs to use",
                                                  options=sorted(proportion_keys, key=str.casefold))
            upkeys = [s for s in select_DEGs if re.search("UP", s)]
            downkeys = [s for s in select_DEGs if re.search("DOWN", s)]

            for s in upkeys:
                ups = premade_dict[s].index.tolist()
                uplist.append(ups)

            for s in downkeys:
                downs = premade_dict[s].index.tolist()
                downlist.append(downs)

            flattenedup = [val for sublist in uplist for val in sublist]  # all the upDEGs in 1 list
            flatteneddown = [val for sublist in downlist for val in sublist]  # all the downDEGs in 1 list

            # converter = defaultdict(list)  # yes it should remain as list do not change
            #
            # for val in flattenedup:
            #     converter[val] = val
            #
            # for val in flatteneddown:
            #     converter[val] = val
            #
            # for o, n in zip(old_symbols, new_symbols):
            #     converter[n] = o
            #
            # for k in flattenedup:  # why this: if not, the converter will also include all the date genes that may not be a DEG
            #     t = converter[k]
            #     if t not in up_enrichr:
            #         up_enrichr.append(t)
            #
            # for k in flatteneddown:
            #     t = converter[k]
            #     if t not in down_enrichr:
            #         down_enrichr.append(t)
            #
            # gene_final = [up_enrichr, down_enrichr]
            gene_final = [flattenedup, flatteneddown]

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
        select_enrichr = st.experimental_memo(gp.enrichr)(gene_list=genelist,
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
        fig.update_layout(title="Enrichr analysis of query genes", title_x=0.5, yaxis={'tickmode': 'linear'})
        st.plotly_chart(fig, use_container_width=True)


    else:
        ups = genelist[0]
        downs = genelist[1]
        remove_dupes_up, remove_dupes_down = [], []
        data_up_trunc, data_down_trunc = None, None
        for f in ups:
            if f not in remove_dupes_up:
                remove_dupes_up.append(f)
        for f in downs:
            if f not in remove_dupes_down:
                remove_dupes_down.append(f)

        if len(ups) != 0:
            gene_final_up = sorted(remove_dupes_up)
            enrichr_up = st.experimental_memo(gp.enrichr)(gene_list=gene_final_up,
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
            enrichr_down = st.experimental_memo(gp.enrichr)(gene_list=gene_final_down,
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
                st.markdown(get_table_download_link(enrichr_download, "enrichr"), unsafe_allow_html=True)
                # st.download_button(label="Download Enrichr dataframe", data=to_excel(enrichr_download),
                #                    file_name="enrichr_updownDEGs_analysis.xlsx")

            fig = make_subplots(rows=2, cols=1, subplot_titles=["Upregulated DEGs", "Downregulated DEGs"],
                                x_title="-logP", shared_xaxes=True)
            fig.add_trace(go.Bar(x=toplot_up['-logP'], y=toplot_up.index,
                                 orientation='h', marker_color="#EF553B"),
                          row=1, col=1)
            fig.add_trace(go.Bar(x=toplot_down['-logP'], y=toplot_down.index,
                                 orientation='h', marker_color="#636EFA"),
                          row=2, col=1)
            fig.update_yaxes(tickmode='linear', tick0=0, dtick=0)


        elif data_up_trunc is not None and data_down_trunc is None:
            toplot_up = data_up_sig.sort_values("-logP", ascending=True).tail(10)
            with enrichr_results_exp:
                st.write("Enrichr with upregulated DEGs")
                st.dataframe(data_up_trunc)
                enrichr_download = [data_up_trunc]
                st.markdown(get_table_download_link(enrichr_download, "enrichr_upDEGs"), unsafe_allow_html=True)
                # st.download_button(label="Download Enrichr dataframe", data=to_excel(enrichr_download),
                #                    file_name="enrichr_upDEGs_analysis.xlsx")
            fig = go.Figure(go.Bar(x=toplot_up['-logP'], y=toplot_up.index, orientation='h', marker_color="#EF553B"))
            fig.update_xaxes(title="-logP")
            fig.update_yaxes(title="Term")

        else:
            toplot_down = data_down_sig.sort_values("-logP", ascending=True).tail(10)
            with enrichr_results_exp:
                st.write("Enrichr with downregulated DEGs")
                st.dataframe(data_down_trunc)
                enrichr_download = [data_down_trunc]
                st.markdown(get_table_download_link(enrichr_download, "enrichr_downDEGs"), unsafe_allow_html=True)
                # st.download_button(label="Download Enrichr dataframe", data=to_excel(enrichr_download),
                #                    file_name="enrichr_downDEGs_analysis.xlsx")
            fig = go.Figure(
                go.Bar(x=toplot_down['-logP'], y=toplot_down.index, orientation='h', marker_color="#636EFA"))
            fig.update_xaxes(title="-logP")
            fig.update_yaxes(title="Term", dtick=0)

        fig.update_layout(title="Enriched Pathways (Top 10), adjpvalue < 0.05", title_x=0.5, showlegend=False,
                          yaxis={'tickmode': 'linear'})
        enrichr_results_exp.info(
            "If nothing was plotted in the bar chart, the pathways did not meet the cutoff of adjusted p-value < 0.05")
        st.plotly_chart(fig, use_container_width=True)


############################################ Prerank and Visualisation #################################################
########## Choose the prerank geneset to use #############
def select_prerank_dataset():
    geneset_dict = {
        "Blood Transcriptomic Modules (BTM)": "BTM.gmt",
        "Reactome 2021": "Reactome.gmt",
        "Vaccinomics (In-house)": "Vaccinomics.gmt",
        # "GO Cellular Component 2021": "GO_Cellular_Component_2021",
        # "GO Biological Process 2021": "GO_Biological_Process_2021",
        # "GO Molecular Function 2021": "GO_Molecular_Function_2021",
        # "KEGG 2021 Human": "KEGG_2021_Human"
    }

    # Selecting genesets (BTM or reactome) to plot from a list
    geneset = prernk_exp.radio(label='Select a geneset', options=geneset_dict.keys(), key="prerank")
    return geneset_dict[geneset]


######### Obtain columns for ranking ###################
def find_cols(df, timepoints):
    col_storage = {}  # yes this pun was intentional
    fc_regex = "log2Fold[-_\s]?Changes?[-_\s]|log2FC"
    df.reset_index(drop=False, inplace=True)
    for tp in timepoints:
        df_genes = df.iloc[:, 0]
        df_FC = df.filter(regex=re.compile(fc_regex, re.IGNORECASE))
        df_FC_tp = df_FC.filter(regex=re.compile(r"[-_\s]?{}[-_\s]?".format(tp), re.IGNORECASE))
        df_FC_tp = df_FC_tp.iloc[:, 0]  ## remove day10s from day1. giving only genes as index and one timepoint

        prerank = pd.concat([df_genes, df_FC_tp], axis=1)
        prerank = prerank.sort_values(prerank.columns[0], ascending=False)
        col_name1 = (prerank.columns.values)[0]
        col_name2 = (prerank.columns.values)[1]
        prerank = prerank.rename(columns={col_name1: "0", col_name2: "1"})

        col_storage[tp] = prerank
    return col_storage


########### Run Prerank #############################
# @st.cache(suppress_st_warning=True)
def execute_prerank(col_dict, geneset):
    prerank_results_dict = {}
    for key, data in col_dict.items():
        running = st.experimental_memo(gp.prerank)(rnk=data,
                             gene_sets=geneset,
                             processes=6,
                             permutation_num=100,  # reduce number to speed up testing
                             outdir=None,
                             seed=6,
                             no_plot=True)

        prerank_results_dict[key] = running

    prerank_results_exp = st.expander("Expand for prerank dataframe", expanded=False)

    prerank_export = []

    for key, result in prerank_results_dict.items():
        terms = result.res2d.index
        ledge = result.res2d.ledge_genes
        nes = result.res2d.nes
        fdr = result.res2d.fdr

        ranked = pd.concat([ledge, nes, fdr], axis=1)
        indexcol = list(ranked.index)
        removed_chars = [i.replace('"', "") for i in indexcol]
        ranked.index = removed_chars
        prerank_export.append(ranked)

        if applyfdr:
            pos_nes = ranked[(ranked["nes"] > 0) & (ranked["fdr"] < 0.05)]
            neg_nes = ranked[(ranked["nes"] < 0) & (ranked["fdr"] < 0.05)]
            neg_nes["negative nes"] = neg_nes["nes"] * -1
        else:
            pos_nes = ranked[ranked["nes"] > 0]
            neg_nes = ranked[ranked["nes"] < 0]
            neg_nes["negative nes"] = neg_nes["nes"] * -1

        pos_nes_sort = pos_nes.sort_values(by=['nes']).tail(10)
        pos_nes_sort.reset_index(inplace=True)
        pos_nes_sort.rename(columns={"index": "Term"}, inplace=True)
        top_pos = len(pos_nes_sort)  # For formatting plot title

        neg_nes_sort = neg_nes.sort_values(by=['negative nes']).tail(10)
        neg_nes_sort.reset_index(inplace=True)
        neg_nes_sort.rename(columns={"index": "Term"}, inplace=True)
        top_neg = len(neg_nes_sort)  # For formatting plot title

        pos = px.bar(pos_nes_sort, x="nes", y="Term", orientation="h",
                     title=f"Positive enrichment ({key})")
        neg = px.bar(neg_nes_sort, x="negative nes", y="Term", orientation="h",
                     title=f"Negative enrichment ({key})")

        pos.update_traces(marker_color="#EF553B")
        neg.update_traces(marker_color="#636EFA")

        pos.update_layout(title=f"Top {top_pos} positively enriched pathways",
                          xaxis_title='NES',
                          font=dict(
                              family='Arial', size=16)
                          )
        neg.update_layout(title=f"Top {top_neg} negatively enriched pathways",
                          xaxis_title='NES',
                          font=dict(
                              family='Arial', size=16)
                          )
        st.plotly_chart(pos)
        st.plotly_chart(neg)

        prerank_results_exp.write(key)
        prerank_results_exp.dataframe(ranked)

    prerank_results_exp.markdown(get_table_download_link(prerank_export, "GSEAprerank"), unsafe_allow_html=True)
    # prerank_results_exp.download_button(label="Download prerank dataframe",
    #                                     data=to_excel(prerank_export),
    #                                     file_name="prerank_analysis.xlsx")

    return


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
                                    options=["volcano plot", "DEGs", "enrichr", "GSEA prerank", "pathway clustergram",
                                             "correlation matrix"])
qc_df(df_dict)

if st.sidebar.checkbox("Show uploaded/demo dataframe"):
    for k, v in cleaned_dict.items():
        st.write(f"**{k} dataframe**", v)

for c in choose_app:
    with st.spinner("🔨Building your dashboard 🔨"):
        time.sleep(0.25)
        if c == "volcano plot":
            list_of_days = timepoints(cleaned_dict)
            colorlist = n_colors(list_of_days)
            dfx = check_log(cleaned_dict)
            volcano(dfx, list_of_days, colorlist)

        elif c == "DEGs":
            list_of_days = timepoints(cleaned_dict)
            colorlist = n_colors(list_of_days)
            dfx = check_log(cleaned_dict)
            degs(dfx, list_of_days, colorlist)
            postdeg = st.sidebar.expander("Expand to plot clustergram based on DEGs", expanded=False)
            with postdeg:
                plot_deg_clust = st.checkbox("Plot DEGs in pathway clustergram", value=False)
            if plot_deg_clust:
                deg_cluster(proportions, dfx)

        elif c == 'enrichr':
            enrichr_exp = st.sidebar.expander("Expand for Enrichr pathway analysis", expanded=False)
            if "DEGs" in choose_app:
                select_dataset = select_enrichr_dataset()
                genelist = genes_used(premade_dict=proportions)
                run_enrichr = enrichr_exp.checkbox("Run Enrichr", value=False)
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
                run_enrichr = enrichr_exp.checkbox("Run Enrichr", value=False)
                if run_enrichr:
                    execute_enrichr(genelist=genelist, select_dataset=select_dataset, use_degs=False)
                else:
                    st.stop()

        elif c == 'GSEA prerank':
            list_of_days = timepoints(cleaned_dict)
            dfx = check_log(cleaned_dict)

            prernk_exp = st.sidebar.expander("Expand for Prerank Analysis", expanded=False)
            select_df = prernk_exp.selectbox("Select your dataset to use", options=dfx.keys())
            if is_tp == 1:
                fillin = "timepoints"
            else:
                fillin = "comparisons"
            use_tps = prernk_exp.multiselect(f"Select {fillin} for prerank", options=list_of_days)

            prerank_dataset = select_prerank_dataset()  # same datasets as enrichr
            prerank_cols = find_cols(dfx[select_df], use_tps)
            applyfdr = prernk_exp.checkbox("Apply FDR < 0.05 cutoff to the bar plots", value=False)
            run_prerank = prernk_exp.checkbox("Run GSEA prerank", value=False)
            if run_prerank:
                execute_prerank(prerank_cols, prerank_dataset)
            else:
                st.stop()

        elif c == "pathway clustergram":
            dfx = check_log(cleaned_dict)
            clust_expand = st.sidebar.expander("Expand for user-input pathway clustergram", expanded=False)
            clustergram(dfx)

        elif c == "correlation matrix":
            dfx = check_log(cleaned_dict)
            corr_matrix(dfx)