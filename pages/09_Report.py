import streamlit as st
import streamlit.components.v1 as components

import jinja2
import pdfkit
import base64
from io import BytesIO

import plotly
import matplotlib
from PIL import Image
from datetime import datetime
import pytz

from helper_functions.downloads import file_downloads


st.header("Report Generation")
# Get time
today = datetime.now(tz=pytz.timezone("Asia/Singapore"))
dt_string = today.strftime("%d %B %Y %I:%M:%S %p (GMT%Z)")

st.info("If you skipped any of the steps in the process, the report will still display the header but not the plots as they were not run.")

# Converting all the plots to output to buffer with BytesIO

all_plots_bytes = {}

for key in ['corr_matrix_plot','cdf_plot', 'barplot', 'enrichr_plots', 'prerank_plots']:
    if key in st.session_state:
        if st.session_state[key] is not None:
            to_bytes = file_downloads.plot_to_bytes(st.session_state[key], graph_module='plotly', format='png')
            all_plots_bytes[key] = to_bytes
        else:
            all_plots_bytes[key] = None
    else:
        all_plots_bytes[key] = None

for key in ['volcano_plots_static', 'clustergram_plot']:
    if key in st.session_state:
        if st.session_state[key] is not None:
            to_bytes = file_downloads.plot_to_bytes(st.session_state[key], graph_module='pyplot', format='png')
            all_plots_bytes[key] = to_bytes
        else:
            all_plots_bytes[key] = None
    else:
        all_plots_bytes[key] = None

# cmatrix = file_downloads.plot_to_bytes(st.session_state['corr_matrix_plot'], graph_module="plotly", format="png")
# volplot = file_downloads.plot_to_bytes(st.session_state['volcano_plots'][0], graph_module="pyplot", format="png") # As a static plot, use the matplotlib one only
# cdf = file_downloads.plot_to_bytes(st.session_state['cdf_plot'], graph_module="plotly", format="png")
# barplot = file_downloads.plot_to_bytes(st.session_state['barplot'], graph_module="plotly", format="png")
# clustergram = file_downloads.plot_to_bytes(st.session_state['clustergram_plot'], graph_module="pyplot", format="png")
# enrichr = file_downloads.plot_to_bytes(st.session_state['enrichr_plots'], graph_module="plotly", format="png")
# prerank = file_downloads.plot_to_bytes(st.session_state['prerank_plots'], graph_module="plotly", format="png")

if 'string_plots' in st.session_state:
    string = {}
    for k,v in st.session_state['string_plots'].items():
        in_memory_file = BytesIO(v)
        im = Image.open(in_memory_file)
        im = im.resize((round(im.size[0]*0.3), round(im.size[1]*0.3)))
        b = BytesIO()
        im.save(b, 'png')
        data = base64.b64encode(b.getbuffer()).decode("ascii")
        string[k] = data
else:
    string = {"None":None}

if 'use_corrected_pval' in st.session_state:
    pval_fmt = "adjusted p-value" if st.session_state['use_corrected_pval'] else "p-value"
else:
    pval_fmt = ""

session_data = {}
if 'enr_showall' in st.session_state:
    if st.session_state['enr_showall']:
        session_data['enr_showX'] = "All"
    else:
        session_data['enr_showX'] = f"Top {st.session_state['enr_showX']}"
else:
    session_data['enr_showX'] = "No"

for key in ['bar_pval', 'bar_fc', 'geneset_enr', 'enr_genedict', 'geneset_prerank', 'prerank_choose_col', 'prerank_showX']:
    if key in st.session_state:
        session_data[key] = st.session_state[key]
    else:
        if key == 'enr_genedict':
            session_data[key] = {}
        else:
            session_data[key] = None

templateLoader = jinja2.FileSystemLoader(searchpath="accessory_files/")
templateEnv = jinja2.Environment(loader=templateLoader)
TEMPLATE_FILE = "output_report_template.html"
template = templateEnv.get_template(TEMPLATE_FILE)


outputText = template.render(date = dt_string,
                             cmatrix = all_plots_bytes['corr_matrix_plot'],
                             volplot = all_plots_bytes['volcano_plots_static'],
                             cdf= all_plots_bytes['cdf_plot'],
                             barplot=all_plots_bytes['barplot'],
                             pval_fmt=pval_fmt,
                             bar_pval=session_data['bar_pval'],
                             bar_fc=session_data['bar_fc'],
                             clustergram=all_plots_bytes['clustergram_plot'],
                             geneset_enr=session_data['geneset_enr'],
                             enr_genedict=session_data['enr_genedict'],
                             enr_showX=session_data['enr_showX'],
                             enrichr=all_plots_bytes['enrichr_plots'],
                             geneset_prerank=session_data['geneset_prerank'],
                             prerank_choose_col=session_data['prerank_choose_col'],
                             prerank_showX=session_data['prerank_showX'],
                             prerank=all_plots_bytes['prerank_plots'],
                             string_dict=string)
html_file = open("STAGES_report.html", 'w')
html_file.write(outputText)
html_file.close()


HTML_inapp = open("STAGES_report.html", 'r', encoding='utf-8')
source_code = HTML_inapp.read()
components.html(source_code, height = 900, scrolling=True)

pdf_out = pdfkit.from_string(outputText, False)
st.download_button("Download STAGES report as PDF here", data=pdf_out, file_name="STAGES_report.pdf", mime='application/octet-stream')