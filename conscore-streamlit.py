import streamlit as st
from Bio import Blast
from Bio import Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from io import StringIO
from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline
import pandas as pd
import tempfile
import os
import numpy as np
from pandas import DataFrame

st.title('ConScore Calculation')

st.header('Submit MSA')

#can insert code here later that will allow file upload & do the alignment internally; will probably not do this bc i think the psa is 'bad' (doesn't match clustal omega)
    
##below is core of code for clustal and consurf uploads & alignment in df

msa_file = st.file_uploader("",type='fa')
if msa_file is not None:
    st.success("MSA file uploaded")
else:
    st.info("please upload your .fasta MSA file")

try:
    temp = msa_file.getvalue().decode("utf-8") ##decodes characters correctly but still has too long file name issue
except AttributeError:
    pass
#st.text(temp)

#declaring variables outside of button if statement so i can access them after the button step
df = ''
df1 = ''
df2 = ''
df_exploded = ''


alignment = AlignIO.read(StringIO(temp), "fasta")
#st.write(alignment)
targets = []
for record in alignment:
    targets.append(record.id)
#st.write(targets)

if st.button('make msa df & freq df'):
    msa_df = pd.DataFrame(alignment, index = targets)
    st.write(msa_df)

    #make frequency df
    aa = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','-']
    freq_df = pd.DataFrame(index = aa, columns = msa_df.columns)
    for i in range(len(msa_df.columns)):
        for a in range(len(aa)):
            freqaa = (msa_df.iloc[:,i] == aa[a]).sum()
            #st.write(freqaa)
            freq_df.iloc[a,i] = freqaa
    st.write(freq_df)

    #make table with consensus scores
    score_df = pd.DataFrame(index = ['consensus aa','score'], columns = msa_df.columns)
    st.write(score_df)
    for i in range(len(freq_df.columns)):
        st.write(str(msa_df.iloc[:,i].mode())[5])
        score_df.iloc[0,i] = str(msa_df.iloc[:,i].mode())[5]
        total = freq_df.iloc[:,i].sum()
        top_aa = freq_df.iloc[:,i].max()
        percent = top_aa/total
        score_df.iloc[1,i] = percent
    st.write(score_df)
        
    #call target and determine it's score
    @st.fragment()
    def frag():
        target = st.text_input('enter your target ID', 'ex: P22259')
        scores = []
        for i in range(len(score_df)):
            if msa_df.iloc[msa_df.index.get_loc(target),i] == score_df.iloc[0,i]:
                scores.append(score_df.iloc[1,i])
        conscore = scores.sum()
        st.write(conscore)
    frag()
        

@st.fragment()
def frag():
    if st.button('create consurf dataframe & align with clustal dataframe'):
        consurf_df = pd.read_csv(consurf_file)
        #st.write(consurf_df)
        consurf_df = consurf_df[['SEQ','COLOR']]
        consurf_df = consurf_df.iloc[1:].reset_index(drop=True)
        #st.write(consurf_df)
        
        #combine dataframes; can concat OR just create the new COLOR one based on presence/absence of letter in each row
        
        #df_combined = pd.concat([df_exploded, consurf_df], axis=1)
        #st.write(df_combined)
                    

        for idx, aa in enumerate(df_exploded['Project Standard Seq']):
            #st.write(aa)
            if aa == '-':
                #st.write(idx)
                gap = idx
                #st.write(gap)
                #st.write(gap - 0.5)
                #consurf_df.loc[gap] = ''
                line = DataFrame({"SEQ": '', "COLOR": 0}, index=[gap -0.5])
                consurf_df = pd.concat([consurf_df, line])
                consurf_df = consurf_df.sort_index().reset_index(drop=True)
        #st.write(consurf_df)
        df_combined = pd.concat([df_exploded, consurf_df], axis=1)
        df_combined = df_combined.iloc[:-1]
        #st.write(df_combined)

        #create new column (evoscore) and fill cells
        df_combined['EvoScore'] = ''
        #st.write(df_combined)
        #st.write(df_combined.dtypes)
        for idx,i in enumerate(df_combined['COLOR']):
            if i < 4:
                df_combined.iloc[idx,5] = 0
            if i >= 4:
                df_combined.iloc[idx,5] = i
            if df_combined.iloc[idx,0] == df_combined.iloc[idx,1]:
                df_combined.iloc[idx,5] = 0
        st.write(df_combined)
        st.write(df_combined.dtypes)
        df_combined['EvoScore'] = df_combined['EvoScore'].astype(float)
        #st.write(df_combined.dtypes)
        evoscore = df_combined['EvoScore'].sum()
        st.write('EvoScore = ' + str(evoscore))
        
        df_combined['Weighted EvoScore'] = ''
        for idx, i in enumerate(df_combined['COLOR']):
            if i < 4:
                df_combined.iloc[idx,6] = 0
            if i >= 4:
                if df_combined.iloc[idx,2] == '*':
                    df_combined.iloc[idx,6] = 0
                if df_combined.iloc[idx,2] == ':':
                    df_combined.iloc[idx,6] = i*0.5
                if df_combined.iloc[idx,2] == '.':
                    df_combined.iloc[idx,6] = i*0.75
                if df_combined.iloc[idx,2] == ' ':
                    df_combined.iloc[idx,6] = i
        st.write(df_combined)
        st.write(df_combined[['Project Standard Seq', 'Target Seq', 'EvoScore', 'Weighted EvoScore']])
        weighted_evoscore = df_combined['Weighted EvoScore'].sum()
        st.write('Weighted EvoScore = ' + str(weighted_evoscore))

            
        
frag()


#@st.fragment()
#def PSA_download():
 #   with open('clustalPSA.aln') as f:
  #      st.download_button('download PSA', f)
#PSA_download()
