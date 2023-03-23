import streamlit as st
import pandas as pd
from io import BytesIO

from helper_functions.downloads import file_downloads
from helper_functions.preprocessing import tested
from helper_functions.session_state import ss
from helper_functions.uploads import fileuploads

# ###################################################### SESSION STATES ##################################################
st.session_state.update(st.session_state)
ss.initialise_state(state_dict = {'file_type':'Ratios and P-values',
                                  'df_in':None,
                                  'demo':False,
                                  'view_df':False,
                                  'anova_dict':None,
                                  'expr_dict':None,
                                  'meta_dict':None,
                                  'meta_in':None,
                                  'meta_excel':None,
                                  'df_excel':None
                                  })
# ########################################################################################################################

st.header("File Uploader")
# ################################################# File Uploader ########################################################
file_opts = st.sidebar.expander("File Upload Options", expanded = True)
use_demo = file_opts.checkbox("Use demo dataset", value=st.session_state['demo'], on_change=ss.binaryswitch, args = ('demo',))

ftype_list = ['RNAseq Counts', 'Log2-Normalised Data', 'Ratios and P-values']
file_type = file_opts.radio(label="Select data type for upload", options = ftype_list,
                             index = ftype_list.index(st.session_state['file_type']))
ss.save_state(dict(file_type = file_type))

df_query = file_opts.file_uploader('Upload your file', type = ['csv', 'txt', 'xlsx'], accept_multiple_files=True, help='Note that excel files take longer to upload')

clear = file_opts.button("Clear cache", on_click=ss.clear_output)
if clear:
    st.experimental_rerun()
# # Save the df_query so that the data remains stored as the conditions will be met via proxy (df_in session state)
# ## Note that the df_in is not actually able to be read in as data directly, but rather is just a number for my conditions to be met
# ## In this case, st.cache_data around the read_xfile method is still required
if len(df_query) != 0:
    ss.save_state(dict(df_in = df_query))
elif st.session_state['demo']:
    ss.save_state(dict(df_in = None))
else:
    st.session_state['df_in'] = st.session_state['df_in']


# # Complicated conditions
# ## 1. have files uploaded and not using demo
# ## 2. no files uploaded and using demo
# ## 3. have files uploaded and using demo
# ## 4. no files uploaded and no demo
if st.session_state['df_in'] is not None and st.session_state['demo'] is False:
    with file_opts:
        cleandict = fileuploads.read_xfile(st.session_state['df_in'], ss_excel = 'df_excel')
    cleandict = fileuploads.capslock_genes(cleandict)

    if file_type == "Ratios and P-values":
        ss.save_state({'anova_dict':cleandict})

    else:
        with file_opts:
            metadata = st.file_uploader(label="Upload gene expression's metadata here", type = ['csv', 'txt', 'xlsx'], accept_multiple_files=True)
            if len(metadata) != 0:
                ss.save_state({'expr_dict':cleandict,'meta_in':metadata})
            else:
                ss.save_state({'expr_dict':cleandict, 'meta_in':st.session_state['meta_in']})

            if st.session_state['meta_in'] is not None:
                meta_dict = fileuploads.read_xfile(st.session_state['meta_in'], ss_excel = 'meta_excel')
                ss.save_state({'meta_dict':meta_dict})

elif st.session_state['df_in'] is None and st.session_state['demo'] is True:
    testdata = pd.read_csv("accessory_files/demo_dataframe_corrected.csv", index_col=0)
    cleandict = {"Demo":testdata}
    ss.save_state(dict(anova_dict = cleandict))

elif st.session_state['df_in'] is not None and st.session_state['demo'] is True:
    st.warning("Deselect demo dataset to use uploaded file or clear cache to remove your file!")

else:
    ss.save_state({'anova_dict':None, 'expr_dict':None, 'meta_dict':None})


view_df = file_opts.checkbox("View demo/uploaded gene expression dataset", value = st.session_state['view_df'], on_change=ss.binaryswitch, args=('view_df', ))

if view_df is True:
    main_expr, meta_expr, anova_expr = st.tabs(['Gene Expression Data', 'Metadata', 'Ratios and P-value Data'])
    exprdict, metadatadict, anovadict = st.session_state['expr_dict'], st.session_state['meta_dict'], st.session_state['anova_dict']

    if exprdict is not None:
        for k,v in exprdict.items():
            main_expr.subheader(k)
            main_expr.dataframe(v)
    else:
        main_expr.info("No gene expression counts uploaded")

    if metadatadict is not None:
        for k,v in metadatadict.items():
            meta_expr.subheader(k)
            meta_expr.dataframe(v)
    else:
        meta_expr.info("No metadata uploaded")

    if anovadict is not None:
        for k,v in anovadict.items():
            anova_expr.subheader(k)
            anova_expr.dataframe(v)
    else:
        anova_expr.info("No ratios and p-values uploaded")

    # _ = {main_expr.subheader(k):main_expr.dataframe(v) for k,v in exprdict.items()} if exprdict is not None else main_expr.info("No gene expression counts uploaded")
    # _ = {meta_expr.subheader(k):meta_expr.dataframe(v) for k,v in metadatadict.items()} if metadatadict is not None else meta_expr.info("No metadata uploaded")
    # _ = {anova_expr.subheader(k):anova_expr.dataframe(v) for k,v in anovadict.items()} if anovadict is not None else anova_expr.info("No ratios and p-values uploaded")

################################### DEBUGGING STAGE #############################################
# ss.initialise_state(dict(viewdf = False, use_demo = False, cleandict = None))

# demodf = file_opts.checkbox("Use Demo", value = st.session_state['use_demo'], on_change=ss.binaryswitch, args = ('use_demo', ))

# if demodf is True:
#     data = pd.read_csv("accessory_files/demo_dataframe_corrected.csv", index_col=0)
#     ss.save_state(dict(cleandict = {'demo':data}))
# else:
#     if len(df_query) != 0:
#         with file_opts:
#             cleandict = st.cache_data(experimental_allow_widgets = True)(fileuploads.read_xfile)(df_query)
#         cleandict = fileuploads.capslock_genes(cleandict)
#         ss.save_state(dict(cleandict = cleandict))
#     elif st.session_state['cleandict'] is None:
#         ss.save_state(dict(cleandict = None))
#     else:
#         st.session_state['cleandict'] = st.session_state['cleandict']

# if file_opts.checkbox("View dataframe", value=st.session_state['viewdf'], on_change=ss.binaryswitch, args = ('viewdf',)):
#     st.write(st.session_state['cleandict'])


################################# SELF-HELP GUIDE TO PREVENT INSANITY ################################
# General Logic for session state:
## For our use case, all widgets and outputs that should be saved, must be initialised
## For widgets, use an on_change function (eg. binaryswitch) or save state right after the input
## If something is de-selected, demo dataset for instance, have to clear the dataset from session state to prevent it from overwriting the user's input

# Logic for file uploader:
## Save the inputs from file uploader to session_state, then when trying to read in the file, use the session_state ver. Not the df_query as it will be None/0 when rerun
