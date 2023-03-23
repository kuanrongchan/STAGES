import streamlit as st
import pandas as pd
from helper_functions.date_gene import qc_df
from helper_functions.session_state import ss

from io import StringIO
from collections import defaultdict

class FileUploads():
    def read_xfile(self, df_query, ss_excel):
        '''
        Parameter
        ---------
        df_query: from the st.file_uploader output
        ss_excel: potential session state key for any possible excel files

        Returns
        -------

        dict| keys as names of file, values as pd DataFrame
        '''

        df_dict = {}
        for d in df_query:
            head, sep, tail = str(d.name).partition(".")
            if tail == 'csv':
                data = st.cache_data(pd.read_csv)(d, index_col=0)
                df_dict[head] = data

            elif tail == 'txt':
                data = st.cache_data(pd.read_csv)(d, sep='\t', index_col=0)
                df_dict[head] = data

            elif tail == 'xlsx':
                x = st.cache_data(pd.read_excel)(d, index_col=0, sheet_name=None, engine='openpyxl')
                if ss_excel == "df_excel":
                    selected_sheet = st.multiselect(label="Select which sheet to read in", options=list(x.keys()), default = st.session_state[ss_excel])
                    if len(selected_sheet) != 0:
                        for i in selected_sheet:
                            df_dict[f"{head}_{i}"] = x[i]
                        ss.save_state({ss_excel: selected_sheet})
                    else:
                        ss.save_state({ss_excel: st.session_state[ss_excel]})
                else:
                    selected_meta = st.multiselect(label="Select which sheet to read in for metadata", options=list(x.keys()), default = st.session_state[ss_excel])
                    if len(selected_meta) != 0:
                        for i in selected_meta:
                            df_dict[f"{head}_{i}"] = x[i]
                        ss.save_state({ss_excel: selected_meta})
                    else:
                        ss.save_state({ss_excel: st.session_state[ss_excel]})
        return df_dict
    
    def capslock_genes(self, df_dict):
        for df in df_dict.values():
            df.index = df.index.astype(str, copy=False) # expand to format actual dates from excel sheets as text
            df.index = df.index.str.upper()

        cleandict = qc_df(df_dict)

        return cleandict
    
    def gmt_to_dict(self, add_geneset_in):
        geneset_dicts = defaultdict(dict)
        for geneset in add_geneset_in:
            upload_to_str = StringIO(st.session_state['add_geneset_in'][0].getvalue().decode("utf-8"))
            fname = geneset.name
            for line in upload_to_str.readlines():
                break_mod = line.split("\t")
                modname = break_mod[0]
                genes = break_mod[2:None]
                geneset_dicts[fname.replace(".gmt","")][modname] = genes
        return dict(geneset_dicts)
    
fileuploads = FileUploads()