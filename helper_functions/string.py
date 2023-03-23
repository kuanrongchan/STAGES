import requests
from PIL import Image
from io import BytesIO

import time
import streamlit as st

class StringDb():
    '''
    This class only helps to query stringdb based on user's genes or DEGs.

    Might expand its capacities to label etc in the future.
    '''
    def string_query(_self, gene_dict):
        string_api_url = "https://version-11-5.string-db.org/api"
        output_format = "highres_image"
        method = "network"
        request_url = "/".join([string_api_url, output_format, method])

        ims = {}
        content_tozip = {}
        for k,v in gene_dict.items():       
            gene_ready  = "%0d".join(v)

            params = {
                "identifiers" : gene_ready,
                "species" : 9606, # species NCBI identifier 
                "network_flavor": "confidence", # show confidence links
                "caller_identity" : "stages", # your app name
                "block_structure_pics_in_bubbles":1
                }
            if len(gene_ready) != 0:
                response = st.cache_resource(requests.post)(request_url, data=params)
                in_memory_file = BytesIO(response.content)
                im = Image.open(in_memory_file)
                content_tozip[k] = response.content
                ims[k] = im
            else:
                ims[k] = None
            time.sleep(1)
        return ims, content_tozip

stages_str = StringDb()