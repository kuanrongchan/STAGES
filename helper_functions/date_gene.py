import pandas as pd
import numpy as np
import re
import time

import base64
from io import BytesIO

import inflect
import dateparser
from datetime import datetime

import streamlit as st
from streamlit_tags import st_tags, st_tags_sidebar

from helper_functions.session_state import ss

def qc_df(df_dict):
    ss.initialise_state({'date_fmt': 'yyyy-dd-mm',
                         'info_stored': 'month-day',
                         'first_mar01_fx': "MTARC1: mitochondrial amidoxime reducing component 1",
                         'first_mar02_fx': "MTARC2: mitochondrial amidoxime reducing component 2"
                         })
    ########################################### HGNC Reference Table ####################################################
    @st.cache_data()
    def clean_ref():
        for_ref = pd.read_csv("accessory_files/hgnc-symbol-check2.csv") # github
        for_ref.reset_index(drop=True,inplace=True)
        for_ref.columns = for_ref.iloc[0,:]
        for_ref.drop(index=0, inplace=True)
        for_ref.drop(columns="Match type", inplace=True)
        for_ref.rename(columns={"Input":"Previous Symbol"}, inplace=True)
        return for_ref
    reference_symbols = clean_ref()
    ######################################### Variables used throughout code #############################################
    p = inflect.engine()
    cleaned_dict = {}

    corrected = {"DEC-01_1st": "DELEC1", "01-DEC_1st":"DELEC1", "MAR-03_1st": "MARCHF3", "03-MAR_1st":"MARCHF3",
                "MAR-04_1st": "MARCHF4", "04-MAR_1st":"MARCHF4", "MAR-05_1st": "MARCHF5", "05-MAR_1st":"MARCHF5",
                "MAR-06_1st": "MARCHF6", "06-MAR_1st":"MARCHF6", "MAR-07_1st": "MARCHF7", "07-MAR_1st":"MARCHF7",
                "MAR-08_1st": "MARCHF8", "08-MAR_1st":"MARCHF8", "MAR-09_1st": "MARCHF9", "09-MAR_1st":"MARCHF9",
                "MAR-10_1st": "MARCHF10", "10-MAR_1st":"MARCHF10", "MAR-11_1st": "MARCHF11", "11-MAR_1st":"MARCHF11",
                "SEP-15_1st": "SELENOF", "15_SEP_1st":"SELENOF", "SEP-01_1st": "SEPTIN1", "01-SEP_1st":"SEPTIN1",
                "SEP-02_1st": "SEPTIN2", "02-SEP_1st":"SEPTIN2", "SEP-03_1st": "SEPTIN3", "03-SEP_1st":"SEPTIN3",
                "SEP-04_1st": "SEPTIN4", "04-SEP_1st":"SEPTIN4", "SEP-05_1st": "SEPTIN5", "05-SEP_1st":"SEPTIN5",
                "SEP-06_1st": "SEPTIN6", "06-SEP_1st":"SEPTIN6", "SEP-07_1st": "SEPTIN7", "07-SEP_1st":"SEPTIN7",
                "SEP-08_1st": "SEPTIN8", "08-SEP_1st":"SEPTIN8", "SEP-09_1st": "SEPTIN9", "09-SEP_1st":"SEPTIN9",
                "SEP-10_1st": "SEPTIN10", "10-SEP_1st":"SEPTIN10", "SEP-11_1st": "SEPTIN11", "11-SEP_1st":"SEPTIN11",
                "SEP-12_1st": "SEPTIN12", "12-SEP_1st":"SEPTIN12", "SEP-13_1st":"SEPTIN7P2", "13-SEP_1st":"SEPTIN7P2",
                "SEP-14_1st": "SEPTIN14", "14-SEP_1st":"SEPTIN14",

                "Dec-01_1st": "DELEC1", "01-Dec_1st":"DELEC1", "Mar-03_1st": "MARCHF3", "03-Mar_1st":"MARCHF3",
                "Mar-04_1st": "MARCHF4", "04-Mar_1st":"MARCHF4", "Mar-05_1st": "MARCHF5", "05-Mar_1st":"MARCHF5",
                "Mar-06_1st": "MARCHF6", "06-Mar_1st":"MARCHF6", "Mar-07_1st": "MARCHF7", "07-Mar_1st":"MARCHF7",
                "Mar-08_1st": "MARCHF8", "08-Mar_1st":"MARCHF8", "Mar-09_1st": "MARCHF9", "09-Mar_1st":"MARCHF9",
                "Mar-10_1st": "MARCHF10", "10-Mar_1st":"MARCHF10", "Mar-11_1st": "MARCHF11", "11-Mar_1st":"MARCHF11",
                "Sep-15_1st": "SELENOF", "15_Sep_1st":"SELENOF", "Sep-01_1st": "SEPTIN1", "01-Sep_1st":"SEPTIN1",
                "Sep-02_1st": "SEPTIN2", "02-Sep_1st":"SEPTIN2", "Sep-03_1st": "SEPTIN3", "03-Sep_1st":"SEPTIN3",
                "Sep-04_1st": "SEPTIN4", "04-Sep_1st":"SEPTIN4", "Sep-05_1st": "SEPTIN5", "05-Sep_1st":"SEPTIN5",
                "Sep-06_1st": "SEPTIN6", "06-Sep_1st":"SEPTIN6", "Sep-07_1st": "SEPTIN7", "07-Sep_1st":"SEPTIN7",
                "Sep-08_1st": "SEPTIN8", "08-Sep_1st":"SEPTIN8", "Sep-09_1st": "SEPTIN9", "09-Sep_1st":"SEPTIN9",
                "Sep-10_1st": "SEPTIN10", "10-Sep_1st":"SEPTIN10", "Sep-11_1st": "SEPTIN11", "11-Sep_1st":"SEPTIN11",
                "Sep-12_1st": "SEPTIN12", "12-Sep_1st":"SEPTIN12", "Sep-13_1st":"SEPTIN7P2", "13-Sep_1st":"SEPTIN7P2",
                "Sep-14_1st": "SEPTIN14", "14-Sep_1st":"SEPTIN14"
                    }

    misidentified = None
    ############################### Path after initial regex ############################################################

    ################ Contains dates and March-01/March-02 and have to be resolved ####################
    def march_resolver(dfs):
        find = [g for g in dfs.index.tolist() if re.search("^MAR-|^APR-|^SEPT?-|^OCT-|^DEC-", g, flags=re.I)]
        global misidentified
        misidentified = ";".join(find)
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
        mar1 = [f for f in found.index.tolist() if re.search("Mar-0?1_1st|0?1-Mar_1st|Mar-0?1_2nd|Mar-0?1_2nd", f, flags=re.I)]
        if len(mar1) !=0:
            mar1_df = found.loc[mar1]

            with each_df_exp:
                st.write(f"**MAR01 Genes: {k} Dataframe**")
                st.info("Genes like MARCH1 and MARC1 have to be differentiated by another identifier (e.g. Gene description) as they are both corrected to Mar-01 in Excel."
                                " Check HGNC symbol reference in the sidebar for reference. 👈")
                st.dataframe(mar1_df.astype(str))

                first_mar01_fx_opts = ["MTARC1: mitochondrial amidoxime reducing component 1", "MARCHF1: membrane associated ring-CH-type finger 1"]
                first_mar01_fx = st.selectbox(f"Select the name and function that {mar1[0]} corresponds to for {k} dataframe",
                                            options=first_mar01_fx_opts,
                                            index = first_mar01_fx_opts.index(st.session_state['first_mar01_fx']))

            # function below can still apply to genes with only 1 MAR-01 gene because the dictionary will only match those found in the data
            if first_mar01_fx == "MTARC1: mitochondrial amidoxime reducing component 1":
                first_mar01 = "MTARC1"
                second_mar01 = "MARCHF1"
            elif first_mar01_fx == "MARCHF1: membrane associated ring-CH-type finger 1":
                first_mar01 = "MARCHF1"
                second_mar01 = "MTARC1"

            corrected["MAR-01_1st"], corrected['Mar-01_1st'], corrected["01-MAR_1st"], corrected["01-Mar_1st"] = first_mar01, first_mar01, first_mar01, first_mar01
            corrected["MAR-01_2nd"], corrected["Mar-01_2nd"], corrected["01-MAR_2nd"], corrected["01-Mar_2nd"] = second_mar01, second_mar01, second_mar01, second_mar01

        mar2 = [f for f in found.index.tolist() if re.search("Mar-0?2_1st|0?2-Mar_1st|Mar-0?2_2nd|02-Mar_2nd", f, flags=re.I)]
        if len(mar2) !=0:
            mar2_df = found.loc[mar2]

            each_df_exp.write(f"**MAR02 Genes: {k} Dataframe**")
            each_df_exp.info(
                "Genes like MARCH2 and MARC2 have to be differentiated by another identifier (e.g. Gene description) as they are both corrected to Mar-01 in Excel."
                " Check old and new HGNC symbols in the sidebar for reference. 👈")
            each_df_exp.dataframe(mar2_df.astype(str))

            first_mar02_fx_opts = ["MTARC2: mitochondrial amidoxime reducing component 2", "MARCHF2: membrane associated ring-CH-type finger 2"]
            first_mar02_fx = each_df_exp.selectbox(f"Select the name and function that {mar2[0]} corresponds to for {k} dataframe",
                                        options=first_mar02_fx_opts,
                                        index=first_mar02_fx_opts.index(st.session_state['first_mar02_fx']))

            if first_mar02_fx == "MTARC2: mitochondrial amidoxime reducing component 2":
                first_mar02 = "MTARC2"
                second_mar02 = "MARCHF2"
            elif first_mar02_fx == "MARCHF2: membrane associated ring-CH-type finger 2":
                first_mar02 = "MARCHF2"
                second_mar02 = "MTARC2"

            corrected["MAR-02_1st"], corrected["Mar-02_1st"], corrected["02-MAR_1st"], corrected["02-Mar_1st"]  = first_mar02, first_mar02, first_mar02, first_mar02
            corrected["MAR-02_2nd"], corrected["Mar-02_2nd"], corrected["02-MAR_2nd"], corrected["02-Mar_2nd"] = second_mar02, second_mar02, second_mar02, second_mar02

        # else:
        #     corrected["MAR-01_1st"], corrected['Mar-01_1st'], corrected["01-MAR_1st"], corrected["01-Mar_1st"] = first_mar01, first_mar01, first_mar01, first_mar01
        #     corrected["MAR-01_2nd"], corrected["Mar-01_2nd"], corrected["01-MAR_2nd"], corrected["01-Mar_2nd"] = second_mar01, second_mar01, second_mar01, second_mar01
        #     corrected["MAR-02_1st"], corrected["Mar-02_1st"], corrected["02-MAR_1st"], corrected["02-Mar_1st"]  = first_mar02, first_mar02, first_mar02, first_mar02
        #     corrected["MAR-02_2nd"], corrected["Mar-02_2nd"], corrected["02-MAR_2nd"], corrected["02-Mar_2nd"] = second_mar02, second_mar02, second_mar02, second_mar02

        found["Gene"] = pd.Series(corrected)  # in order to rename just change this to found.rename(corrected)
        found.reset_index(drop=True, inplace=True)  # remove the gene index with the dates
        found.rename(columns={"Gene": "gene"}, inplace=True)  # rename the incoming column to be used as index
        found.set_index('gene', inplace=True)  # set the index of these date genes using corrected names
        dfs = dfs.drop(index=find)  # drop the date genes from the main df
        df2 = pd.concat([dfs, found], axis=0)  # join these genes back to the main df
        df2.sort_index(axis=0, ascending=True, inplace=True)  # sort alphabetically
        df2.reset_index(drop=False, inplace=True)
        df2.rename(columns={'index': index_name}, inplace=True)
        df2.set_index(index_name, inplace=True)
        cleaned_dict[k] = df2
        return

    ############ Contains dates but no march-01/march-02 and thus nothing to resolve ##############
    def date_resolver(df, date_search):
        global misidentified
        misidentified = ";".join(date_search)
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
    def numeric_date(k,df,numdate):
        num_exp = st.expander(f"Expand if {k}'s date format is numerical (eg. yyyy/mm/dd)")

        date_fmt_opts = ["yyyy-dd-mm", "yyyy-mm-dd", "dd-mm-yyyy", "mm-dd-yyyy"]
        info_stored_opts = ['month-year', 'month-day']

        date_fmt = num_exp.radio(f"Select the format that {k} dataframe is in",
                    options=date_fmt_opts, index = date_fmt_opts.index(st.session_state['date_fmt']))
        ss.save_state({'date_fmt':date_fmt})

        info_stored = num_exp.radio(f"Select how {k}'s dates should be read to derive gene names (Hover '?' for help)",
                                    options= info_stored_opts,
                                    help='For example, 2001-03-09 (yyyy-mm-dd) may either be Mar-01 (MARCHF1) or Mar-09 (MARCHF9).',
                                    index=info_stored_opts.index(st.session_state['info_stored']))
        ss.save_state({'info_stored':info_stored})

        num_exp.info("If you're unsure about the above option, check the converted dataframe and select 'month-year.' "
                    "We recommend you to check the converted dataframe to ensure that the dates are converted correctly. If unsuccessful, <NA> symbols will populate at the bottom of the converted dataframe.")
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
        found = df.iloc[flattened,:]
        found.set_index(index_nm, inplace=True)
        num_exp.write(f"**{k} dataframe**")
        num_exp.dataframe(found)

        if date_fmt == "yyyy-dd-mm":
            extracted = {n:(dateparser.parse(n, date_formats=["%Y-%d-%m"])).strftime(strfmt) for n in numdate}
        elif date_fmt == "yyyy-mm-dd":
            extracted = {n:(dateparser.parse(n, date_formats=["%Y-%m-%d"])).strftime(strfmt) for n in numdate}
        elif date_fmt == "dd-mm-yyyy":
            extracted = {n:(dateparser.parse(n, date_formats=["%d-%m-%Y"])).strftime(strfmt) for n in numdate}
        elif date_fmt == "mm-dd-yyyy":
            extracted = {n:(dateparser.parse(n, date_formats=["%m-%d-%Y"])).strftime(strfmt) for n in numdate}

        global misidentified
        misidentified = ";".join(extracted.keys())
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

        global misidentified
        misidentified = ";".join(corrected.keys())
        df.rename(index=corrected, inplace=True)
        cleaned_dict[k] = df
        return

    ################################################# Code Flow ##########################################################

    # since it's quite likely that old symbols don't exist together with dates (because once opened in excel all are dates),
    # this is an all-or-none approach where if date search picks up sth, old search will be empty
    # if both lists are empty, nothing is wrong, then cleaned dict[k] = df_dict[k], where k is df with no error

    ismar, isnums = 0, 0

    for k,df in df_dict.items():
        date_search = [g for g in df.index.tolist() if re.search("^Mar-|^Apr-|^Sept?-|^Oct-|^Dec-", g, flags=re.I)] # dates
        old_symbols = list(reference_symbols['Previous Symbol'])
        old_search = list(set(df.index.tolist()).intersection(set(old_symbols))) # easy way to find old symbols in df index
        if len(date_search) != 0:
            march_search = [m for m in date_search if re.search("^Mar-0?1|^0?1-Mar|^Mar-0?2|^0?2-Mar", m, flags=re.I)] # only march genes
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
            numdate = [g for g in df.index.tolist() if re.search("^\d*[-/]?\W", g, flags=re.I)]
            if len(numdate) != 0:
                isnums += 1
                if isnums == 1:
                    st.subheader("Resolve Date Format")
                renamed = numeric_date(k,df,numdate)
                march_search = [m for m in renamed.index.tolist() if re.search("^Mar-0?1|^0?1-Mar|^Mar-0?2|^0?2-Mar", m, flags=re.I)]  # only march genes
                generic_date = [g for g in renamed.index.tolist() if re.search("^Mar-|^Apr-|^Sept?-|^Oct-|^Dec-", g, flags=re.I)]
                if len(march_search) != 0:
                    ismar += 1
                    if ismar == 1:
                        st.subheader("Resolve Duplicate Gene Symbols")
                        march_resolver(renamed)  # will have all the date genes (like it is now)
                else:
                    date_resolver(renamed, generic_date)  # will exclude march
            else:
                noerr = st.success(f"No errors detected for {k} dataframe")
                cleaned_dict[k] = df
                time.sleep(1)
                noerr.empty()
        else:
            noerr = st.success(f"No errors detected for {k} dataframe")
            cleaned_dict[k] = df
            time.sleep(1)
            noerr.empty()
    return cleaned_dict