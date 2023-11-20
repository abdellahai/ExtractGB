import streamlit as st
from Bio import SeqIO
import pandas as pd
import io
from Bio.SeqUtils import GC
from Bio.Seq import Seq
@st.cache_data
def convert_df(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return df.to_csv().encode('utf-8')


def GetGens (records):
    concatenated_sequence = Seq("")
    data = [
        {'Gene group':'NADH dehydrogenase', 'Genes':[], 'Number of genes':0, 'Size':0},
        {'Gene group':'Succinate dehydrogenase', 'Genes':[], 'Number of genes':0, 'Size':0},
        {'Gene group':'Ubichinol cytochrom c reductase', 'Genes':[], 'Number of genes':0, 'Size':0},
        {'Gene group':'Cytochrom c oxidase', 'Genes':[], 'Number of genes':0, 'Size':0},
        {'Gene group':'ATP synthase', 'Genes':[], 'Number of genes':0, 'Size':0},
        {'Gene group':'Cytochrom biogenesis', 'Genes':[], 'Number of genes':0, 'Size':0},
        {'Gene group':'Ribosomal proteins SSU', 'Genes':[], 'Number of genes':0, 'Size':0},
        {'Gene group':'NRibosomal proteins LSU', 'Genes':[], 'Number of genes':0, 'Size':0},
        {'Gene group':'Maturase', 'Genes':[], 'Number of genes':0, 'Size':0},
        {'Gene group':'ORFs', 'Genes':[], 'Number of genes':0, 'Size':0},
        {'Gene group':'photosystem', 'Genes':[],'Number of genes':0,'Size':0},
        {'Gene group':'RNA polymerase', 'Genes':[],'Number of genes':0,'Size':0},
        {'Gene group':'RuBisCO', 'Genes':[],'Number of genes':0,'Size':0},
        {'Gene group':'Transfer RNA', 'Genes':[], 'Number of genes':0, 'Size':0},
        {'Gene group':'Ribosomal RNA', 'Genes':[], 'Number of genes':0, 'Size':0},
        {'Gene group':'Other genes', 'Genes':[], 'Number of genes':0, 'Size':0},
    ]
    for record in records:
        dna_sequence = record.seq
        for feature in record.features:
            if 'gene' in feature.qualifiers:
                gene = feature.qualifiers['gene'][0]
                if 'nad' in gene or 'ndh' in gene:
                    data[0]['Genes'].append(gene)
                    data[0]['Size'] += len(feature.extract(record.seq))
                elif 'sdh' in gene:
                    data[1]['Genes'].append(gene)
                    data[1]['Size'] += len(feature.extract(record.seq))
                elif 'cob' in gene or 'pet' in gene:
                    data[2]['Genes'].append(gene)
                    data[2]['Size'] += len(feature.extract(record.seq))
                elif 'cox' in gene:
                    data[3]['Genes'].append(gene)
                    data[3]['Size'] += len(feature.extract(record.seq))
                elif 'atp' in gene:
                    data[4]['Genes'].append(gene)
                    data[4]['Size'] += len(feature.extract(record.seq))
                elif 'ccm' in gene or 'ccs' in gene:
                    data[5]['Genes'].append(gene)
                    data[5]['Size'] += len(feature.extract(record.seq))
                elif 'rps' in gene:
                    data[6]['Genes'].append(gene)
                    data[6]['Size'] += len(feature.extract(record.seq))
                elif 'rpl' in gene:
                    data[7]['Genes'].append(gene)
                    data[7]['Size'] += len(feature.extract(record.seq))
                elif 'mat' in gene:
                    data[8]['Genes'].append(gene)
                    data[8]['Size'] += len(feature.extract(record.seq))
                elif 'orf' in gene:
                    data[9]['Genes'].append(gene)
                    data[9]['Size'] += len(feature.extract(record.seq))
                elif 'psa' in gene or 'psb' in gene:
                    data[10]['Genes'].append(gene)
                    data[10]['Size']+= len(feature.extract(record.seq))
                elif 'rpo' in gene :
                    data[11]['Genes'].append(gene)
                    data[11]['Size']+= len(feature.extract(record.seq))
                elif 'rbc' in gene :
                    data[12]['Genes'].append(gene)
                    data[12]['Size']+= len(feature.extract(record.seq))
                elif 'rrn' in gene:
                    data[14]['Genes'].append(gene)
                    data[14]['Size'] += len(feature.extract(record.seq))
                else:
                    data[15]['Genes'].append(gene)
                    data[15]['Size'] += len(feature.extract(record.seq))
            
            elif 'product' in feature.qualifiers:
                product  = feature.qualifiers['product'][0]
                if "TRN" in product.upper():
                    data[13]['Genes'].append(product)
                    data[13]['Size'] += len(feature.extract(record.seq))
                elif "ribosomal RNA" in product:
                    data[14]['Genes'].append(product)
                    data[14]['Size'] += len(feature.extract(record.seq))
        concatenated_sequence += dna_sequence
    for d in data:
        d['Genes'] = list(set(d['Genes']))                
    for d in data:
        d['Number of genes'] = len(d['Genes'])
    features = pd.DataFrame(data)    
    trnas = features[features['Gene group'] == 'Transfer RNA']['Number of genes'].sum()
    rrnas = features[features['Gene group'] == 'Ribosomal RNA']['Number of genes'].sum()
    pcg = features['Number of genes'].sum() - (trnas + rrnas)
    description = pd.DataFrame(
        [
            {'Parameter':'GC%', 'Value': GC(concatenated_sequence)},
            {'Parameter':'Coding region ratio%', 'Value': (features["Size"].sum()/len(concatenated_sequence)) * 100},
            {'Parameter':'Number of Protein coding genes', 'Value': pcg},
            {'Parameter':'Number of Transfer RNAs', 'Value': trnas},
            {'Parameter':'Number of Ribosomal RNAs', 'Value': rrnas}
            
        ]
    )
    return features,description
def main():
    st.set_page_config(layout='wide')
    st.title("GenBank File Classifier")
    hide_streamlit_style = """
            <style>
            #MainMenu {visibility: hidden;}
            footer {visibility: hidden;}
            </style>
            """
    st.markdown(hide_streamlit_style, unsafe_allow_html=True) 

    # Upload GenBank file
    uploaded_file = st.file_uploader("Upload a GenBank file", type=["gbk", "gb"])

    if uploaded_file is not None:
        # Read the content of the file in chunks
        chunk_size = 1024  # Adjust the chunk size as needed
        content = ""

        while True:
            chunk = uploaded_file.read(chunk_size)
            if not chunk:
                break
            content += chunk.decode('utf-8')

        records = SeqIO.parse(io.StringIO(content), "genbank")
        features, description = GetGens(records)
        with st.container():
            col1, _ , col2 = st.columns([3,0.01,1.5])
            with col1:
                st.subheader('Genes table')
                st.dataframe(features, hide_index=True)
                csv_features = convert_df(features)
                st.download_button(
                    label="Download genes table",
                    data=csv_features,
                    file_name='features.csv',
                    mime='text/csv',
                )
            with col2:
                st.subheader('Mitochondrion description')
                st.dataframe(description, hide_index=True)
                csv_description = convert_df(description)
                st.download_button(
                    label="Download description table",
                    data=csv_description,
                    file_name='description.csv',
                    mime='text/csv',
                )
        
        

# Run the Streamlit app
if __name__ == "__main__":
    main()
