import streamlit as st

class States():
    def initialise_state(self, state_dict):
        for k,v in state_dict.items():
            if k not in st.session_state:
                st.session_state[k] = v
            else:
                pass

    def save_state(self, state_dict):
        for k,v in state_dict.items():
            st.session_state[k] = v
    
    def binaryswitch(self, key):
        if st.session_state[key] is True:
            st.session_state[key] = False
        else:
            st.session_state[key] = True

    def clear_output(self):
        for k in st.session_state:
            del st.session_state[k]
ss = States()