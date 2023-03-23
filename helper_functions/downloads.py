import pandas as pd
import numpy as np
import time
import math
import re
import base64
from PIL import Image
from io import BytesIO
import zipfile

import streamlit as st
from streamlit_tags import st_tags, st_tags_sidebar


class DLs():
    def convert_df(self, df):
        return df.to_csv().encode('utf-8')

    def to_excel(self, df, sheetnames = None):
        output = BytesIO()
        writer = pd.ExcelWriter(output, engine='xlsxwriter')
        if sheetnames is not None:
            for d, i in zip(df, sheetnames):
                d.to_excel(writer, sheet_name = i[:30])
        else:
            for d, i in zip(df, range(len(df))):
                d.to_excel(writer, sheet_name=f'Sheet {i + 1}')
        writer.close()
        processed_data = output.getvalue()
        return processed_data


    def get_table_download_link(self, df, purpose, sheetnames = None):  # downloads without needing to reset the whole scripts
        """Generates a link allowing the data in a given panda dataframe to be downloaded
        in:  dataframe
        out: href string
        """
        val = self.to_excel(df, sheetnames=sheetnames)
        b64 = base64.b64encode(val)
        return f'<a style = "border:1px solid #31333f33; border-radius:0.25rem; background-color:#f9f9fb; text-decoration:none; color:black; padding:0.50rem 0.75rem" href="data:application/octet-stream;base64,{b64.decode()}" download="{purpose}.xlsx">' \
            f'📥 Download {purpose} as Excel file 📥</a>'  # decode b'abc' => abc

    # def zip_file(dfs, keys):
    #     with zipfile.ZipFile("cleanedfiles.zip", 'w') as compress:
    #         for df, k in zip(dfs, keys):
    #             file = df.to_csv().encode('utf-8')
    #             compress.writestr(f"cleaned_{k}.csv", file)

    #     with open("cleanedfiles.zip", "rb") as fp:
    #           btn = st.download_button(
    #               label="Download ZIP",
    #               data=fp,
    #               file_name="cleanedfiles.zip",
    #               mime="application/octet-stream"
    #               )
    def zip_imgs(self, imgsave_dict, zipfilename = "STRING_network.zip"):
        '''
        Parameters
        ----------
        imgsave_dict: dict | keys containing filename (without extension), values containing response.content from STRING
        zipfilename: str | name for the zip file
        '''
        with zipfile.ZipFile(zipfilename, 'w') as compress: # First create a zip file
            for k,v in imgsave_dict.items():
                buf = BytesIO(v) # Get the bytes from this PNG file
                compress.writestr(f"{k}.png", buf.getvalue()) # Add the bytes object.getvalue() to the zip
        
        with open(zipfilename, "rb") as fp: # Open to allow download
              btn = st.download_button(
                  label="Download STRING Network as ZIP",
                  data=fp,
                  file_name=zipfilename,
                  mime="application/zip"
                  )


    def create_pdf(self, fig, fn, graph_module = "pyplot"):
        buf = BytesIO()
        if graph_module == 'pyplot':
            fig.savefig(buf, format = 'pdf', bbox_inches = 'tight')
        elif graph_module == 'plotly':
            fig.write_image(file = buf, format = 'pdf', engine='kaleido')
        st.download_button(
            label = f"Download {fn.replace('_', ' ')} as pdf",
            data = buf,
            file_name = f"{fn}.pdf",
            mime='application/pdf'
        )

    def plot_to_bytes(self, fig, graph_module="pyplot", format="png"):
        buf = BytesIO()
        if graph_module == "pyplot":
            fig.savefig(buf, format = format, bbox_inches="tight", dpi=300)
        elif graph_module == 'plotly':
            fig.write_image(file = buf, format = format, scale=3)
        
        data = base64.b64encode(buf.getbuffer()).decode("ascii")
        return data

file_downloads = DLs()