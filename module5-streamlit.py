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

st.title('EvoScore Calculation')

st.header('Submit Sequences for PSA')

#can insert code here later that will allow file upload & do the alignment internally; will probably not do this bc i think the psa is 'bad' (doesn't match clustal omega)
    
##below is core of code for clustal and consurf uploads & alignment in df

psa_file = st.file_uploader("",type='clustal_num')
if psa_file is not None:
    st.success("PSA file uploaded")
else:
    st.info("please upload your clustal .clustal file")


consurf_file = st.file_uploader('',type='csv')
if consurf_file is not None:
    st.success('consurf file uploaded')
else:
    st.info('please upload the consurf excel file')


#temp = psa_file.read() ##adds 'b in front of file & other character issues (adds /n etc)
bytes = psa_file.getvalue() ##adds 'b in front of file & other character issues (adds /n etc)
#st.write(bytes)
temp = psa_file.getvalue().decode("utf-8") ##decodes characters correctly but still has too long file name issue
#st.text(temp)

#declaring variables outside of button if statement so i can access them after the button step
df = ''
df1 = ''
df2 = ''
df_exploded = ''

if st.button('read in clustal alignment file'):
    ps_line = ""
    target_line = ""
    temp_split = temp.splitlines()
    # Process the captured Clustal output text line by line
    for line in temp_split[1:]:
        # The conservation line is identifiable by its spacing
        if 'P22259' in line:
            seq = line[26:86]
            ps_line += seq
        if line[:1].isalpha():
            if 'P22259' not in line:
                seq = line[26:86]
                target_line += seq
    ps_line = "".join(char for char in ps_line if not char.isdigit())
    target_line = "".join(char for char in target_line if not char.isdigit())
    #ps_line = ps_line.strip()
    #target_line = target_line.strip()
    st.text(ps_line)
    st.text(target_line)
            
    #alignment = AlignIO.read(temp, 'clustal')
    #alignment = AlignIO.read(StringIO(temp), "clustal")
    #alignment = AlignIO.read('ctei_clustal.aln', 'clustal')
    #st.write(alignment)
    
    #convert clustal alignment to individual sequence strings
    
    #seq1 = str(alignment[0].seq)
    #seq2 = str(alignment[1].seq)
    
    #st.write(seq1)
    #st.write(seq2)
    
    #convert strings to pandas dataframe
    
    data = {'Target Seq': [target_line],
                'Project Standard Seq': [ps_line]}
    df = pd.DataFrame(data)
    #st.write(df)
    df1 = df['Target Seq'].str.split('').explode().reset_index(drop=True)
    #st.write(df1)
    df2 = df['Project Standard Seq'].str.split('').explode().reset_index(drop=True)
    #st.write(df2)
    df_exploded = pd.concat([df1, df2], axis=1)
    #st.write(df_exploded)
    #df_exploded['color'] = 0
    #df_exploded = df_exploded.iloc[1:].reset_index(drop=True) #moving this to after the conservation symbols are added
    #st.write(df_exploded)
       
    #declaring more variables outside button if statement
    consurf_df = ''

    conservation_line = ""
    # Process the captured Clustal output text line by line
    for line in temp_split:
        # The conservation line is identifiable by its spacing
        if any(char in "*:. " for char in line):
            # The line contains symbols, but not sequence data
            if line.startswith(" "):
                # Extract only the symbols, stripping whitespace
                symbols = line[26:] #.strip() #not using .strip() bc i need spaces; dropping the first 26 bc that seems to be the number of excess at the front in the file (should be consistent btwn alignments)
                conservation_line += symbols
    #st.write(conservation_line)
    co_data = {'conservation': [conservation_line]}
    df_symbols = pd.DataFrame(co_data)
    df_symbols = df_symbols['conservation'].str.split('').explode().reset_index(drop=True)
    df_exploded = pd.concat([df_exploded, df_symbols], axis=1)
    df_exploded = df_exploded.iloc[1:-1].reset_index(drop=True)
    #st.write(df_exploded)
    

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
