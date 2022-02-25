import streamlit as st
import numpy as np
import pandas as pd

# Randomly fill a dataframe and cache it
@st.cache(allow_output_mutation=True)
def get_dataframe():
    return pd.DataFrame(
        np.random.randn(50, 20),
        columns=('col %d' % i for i in range(20)))


df = get_dataframe()

# Create row, column, and value inputs
row = st.number_input('row', max_value=df.shape[0])
col = st.number_input('column', max_value=df.shape[1])
value = st.number_input('value')

# Change the entry at (row, col) to the given value
df.values[row][col] = value

# And display the result!
st.dataframe(df)