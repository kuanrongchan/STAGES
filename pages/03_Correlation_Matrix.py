import streamlit as st
from helper_functions.session_state import ss
from helper_functions.corr_matrix import cmatrix
from helper_functions.preprocessing import tested

st.session_state.update(st.session_state)

ss.initialise_state({'corr_mtd': "pearson",
                     'corr_matrix_plot':None})

st.header("Correlation Matrix")
corr_exp = st.sidebar.expander("Expand for correlation matrix", expanded=False)
corr_opts = ['pearson', 'kendall', 'spearman', 'phik']

if "log_dict_ready" in st.session_state:
    if st.session_state['log_dict_ready'] is not None:
        log_dict = st.session_state['log_dict_ready']
        
        corr_mtd = corr_exp.selectbox("Choose the correlation coefficient to use", options=corr_opts, format_func = lambda x: x.title(), index = corr_opts.index(st.session_state['corr_mtd']))
        ss.save_state({'corr_mtd':corr_mtd})

        mtx = cmatrix.corr_matrix(log_dict, method = st.session_state['corr_mtd'])
        ss.save_state({'corr_matrix_plot':mtx})
        st.info("If the plot is too small, please hover over the plot and click the expand button on the top right corner of the plot.")
        st.plotly_chart(st.session_state['corr_matrix_plot'], theme = None, use_container_width=True)

        st.subheader("Pre-processed Data")
        data_exp = st.expander("Expand to show pre-processed data", expanded=False)
        for k,v in st.session_state['ready'].items():
            data_exp.write(f"**{k}**")
            data_exp.dataframe(v)
    
    else:
        st.error("Perhaps you forgot to run through the pre-processing step?")
else:
    st.error("Perhaps you forgot to run through the pre-processing step?")
