# Housekeeping
from helper_functions.session_state import ss

# Helper functions
from helper_functions.degs import preDE, DE
from helper_functions.downloads import file_downloads

# Streamlit
import streamlit as st

st.session_state.update(st.session_state)
ss.initialise_state({'reset_volcano':False,
                    'xaxes_volcano':(0.0,0.0),
                     'yaxes_volcano':0.0,
                     'interactive_volcano':False,
                     'volcano_plots_static':None,
                     'volcano_plots_interactive':None,
                     'cdf_pthresh':0.05,
                     'cdf_linemode':'lines',
                     'cdf_plot':None,
                     'bar_pval':0.05,
                     'bar_fc':1.30,
                     'bar_width':800,
                     'bar_height':600,
                     'barplot':None,
                     'degs':None})

st.header("Differential Expression Analysis")

try:
    bar_t, cdf_t,  volcano_t, data_t = st.tabs(["Differential Expression Bar Plots", "Cumulative Distribution Function", "Volcano Plot", "DEG identity"])
    use_corrected_pval_fmt = "adjusted p-value" if st.session_state['use_corrected_pval'] else "p-value"
    ######### BAR PLOT #################
    deg_opts = st.sidebar.expander("Differential expression bar plot options", expanded=True)
    bar_pval = deg_opts.number_input(f"Choose {use_corrected_pval_fmt} threshold for differentially expressed genes", min_value = 0.00, max_value = 1.00, step = 0.01, value = st.session_state['bar_pval'])
    bar_fc = deg_opts.number_input(label="Adjust fold-change cutoff here ", value=st.session_state['bar_fc'], min_value=0.0, max_value=20.0, step=0.1)
    bar_width = deg_opts.number_input(label="Adjust bar plot width (in px)", min_value=300, max_value=1200, value=st.session_state['bar_width'], step=50)
    bar_height = deg_opts.number_input(label="Adjust bar plot height (in px)", min_value=300, max_value=1200, value=st.session_state['bar_height'], step=50)

    ss.save_state({'bar_pval':bar_pval,
                   'bar_fc':bar_fc,
                   'bar_width':bar_width,
                   'bar_height':bar_height})

    stacked1, proportions = DE.degs(st.session_state['log_dict_ready'],
                                    st.session_state['comparisons'],
                                    pval_cutoff=st.session_state['bar_pval'],
                                    fc_cutoff=st.session_state['bar_fc'],
                                    u_width=st.session_state['bar_width'],
                                    u_height=st.session_state['bar_height'],
                                    use_corrected_pval=st.session_state['use_corrected_pval'])

    ss.save_state({'degs':proportions,
                   'barplot':stacked1})

    with bar_t:
        st.plotly_chart(stacked1, theme=None, use_container_width=False)
        file_downloads.create_pdf(stacked1, "stacked_DEG_plot", graph_module='plotly')

    ########## CDF PLOT ###############
    line_options = ["lines", "markers", "lines+markers"]
    cdf_exp = st.sidebar.expander("Cumulative distribution function options", expanded=True)
    cdf_pthresh = cdf_exp.number_input(f"Choose {use_corrected_pval_fmt} threshold for cumulative distribution function plot",
                                    min_value = 0.00,
                                    max_value = 1.00,
                                    step = 0.01,
                                    value = st.session_state['cdf_pthresh'])
    cdf_linemode = cdf_exp.selectbox("Choose line mode", options=line_options,
                                    format_func=lambda x: x.title().replace("+", " & "),
                                    index = line_options.index(st.session_state['cdf_linemode']))
    ss.save_state({'cdf_pthresh': round(cdf_pthresh,2), 'cdf_linemode':cdf_linemode})
    cdf_plot = preDE.deg_cdf(st.session_state['ready'],
                            st.session_state['comparisons'],
                            pval=st.session_state['cdf_pthresh'],
                            markermode=st.session_state['cdf_linemode'],
                            use_corrected_pval=st.session_state['use_corrected_pval'])
    ss.save_state({'cdf_plot':cdf_plot})

    with cdf_t:
        st.plotly_chart(st.session_state['cdf_plot'], theme=None, use_container_width=True)
        file_downloads.create_pdf(st.session_state['cdf_plot'], "cumulative_density_DEGs", graph_module='plotly')

######### VOLCANO PLOT ##################
    vol_opts = st.sidebar.expander("Volcano plot options", expanded=True)

    reset = vol_opts.checkbox("Reset to default settings", value=st.session_state['reset_volcano'], on_change=ss.binaryswitch, args=('reset_volcano', ))
    if reset:
        ss.save_state({'xaxes_volcano':(0.0, 0.0), 'yaxes_volcano': 0.0})

    xaxes = vol_opts.slider("Choose log2 fold-change boundaries for volcano plot",
                                help="The app will plot the values between the user-set range",
                                min_value=-10.0, max_value=10.0, step=0.1, value= st.session_state['xaxes_volcano'])
    ss.save_state({'xaxes_volcano':xaxes})

    yaxes = vol_opts.number_input(f"Choose -log10 {use_corrected_pval_fmt} boundaries for volcano plot",
                                help="The app will plot the values less than or equal to the user-set value",
                                min_value=0.0, max_value=50.0, step=0.1, value=st.session_state['yaxes_volcano'])
    ss.save_state({'yaxes_volcano':yaxes})

    interactive_volcano = vol_opts.checkbox(label="Show interactive volcano plot", value=st.session_state['interactive_volcano'],
                                            help="Facilitates gene name display on hover. This may cause lag", on_change=ss.binaryswitch, args=('interactive_volcano', ))

    vol_plot, iplot = preDE.volcano(
        user_log=st.session_state['log_dict_ready'],
        comparison_dict=st.session_state['comparisons'],
        xaxes=st.session_state['xaxes_volcano'],
        yaxes=st.session_state['yaxes_volcano'],
        interactive_volcano= st.session_state['interactive_volcano'],
        use_corrected_pval=st.session_state['use_corrected_pval']
        )
    ss.save_state({'volcano_plots_static':vol_plot,
                   'volcano_plots_interactive':iplot})

    with volcano_t:
        if st.session_state['interactive_volcano']:
            st.pyplot(st.session_state['volcano_plots_static'])
            file_downloads.create_pdf(st.session_state['volcano_plots_static'], "volcano_plot", graph_module='pyplot')

            st.plotly_chart(st.session_state['volcano_plots_interactive'], theme = None, use_container_width=True)
            file_downloads.create_pdf(st.session_state['volcano_plots_interactive'], "interactive_volcano_plot", graph_module="plotly")

        else:
            st.pyplot(st.session_state['volcano_plots_static'])
            file_downloads.create_pdf(st.session_state['volcano_plots_static'], "volcano_plot", graph_module='pyplot')

    ####### DATA OUTPUT ###############
    with data_t:
        st.info("Users may select the gene names within this dataframe and copy them for downstream analysis.")
        for k, v in st.session_state['degs'].items():
            st.write(f"**{k}**", v)
        st.download_button(label="Download DEGs", data=file_downloads.to_excel(st.session_state['degs'].values()), file_name="DEGs.xlsx")

except AttributeError:
    st.error("Perhaps you forgot to run through the pre-processing page?")