import subprocess
try:
    subprocess.check_call(["pip", "install", 'biopython'])
    print(f"Successfully installed biopython")
except subprocess.CalledProcessError:
    print(f"Failed to install biopython")
import streamlit as st
from Bio import Entrez, SeqIO
import pandas as pd
import zipfile
import io
import base64

# Function to fetch protein coding sequences from GenBank
def fetch_sequences(accession_numbers):
    seqdic = {'cds tag':[], 'organism':[], 'cds seq': []}
    for accession in accession_numbers:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        # Extract protein coding sequences
        for feature in record.features:
            if feature.type == "CDS":
                cds_name = feature.qualifiers.get("gene", ["UnknownGene"])[0]
                species_name = record.annotations["organism"]
                seqdic['cds tag'].append(cds_name)
                seqdic['organism'].append(species_name)
                seqdic['cds seq'].append(str(feature.extract(record.seq)))
    seqdf = pd.DataFrame(seqdic)
    return seqdf

def find_shared(seqdf: pd.DataFrame):
    sequences = {}
    seqdf['cds tag'] = seqdf['cds tag'].astype('category')
    for cds_tag in seqdf['cds tag'].cat.categories.tolist():
        sub_df = seqdf[seqdf['cds tag'] == cds_tag].reset_index()
        if sub_df.shape[0] == 16:
            sequences[cds_tag] = []
            for _, r in sub_df.iterrows():
                sequences[cds_tag].append(SeqIO.SeqRecord(
                    seq=r['cds seq'],
                    id=f'{cds_tag}_{r["organism"].replace(" ","_")}'
                ))
    return sequences
def create_zip_file(sequences: dict):
    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, 'a') as zipf:
        for cds, seq_records in sequences.items():
            if len(seq_records) == 16:
                # Merge all seq_records into a single multi-FASTA string
                multi_fasta_str = "\n".join(f">{seq_record.id}\n{seq_record.seq}" for seq_record in seq_records)
                # Write the multi-FASTA string to the zip file
                zipf.writestr(f"{cds}_merged.fasta", multi_fasta_str)
    return zip_buffer

# Streamlit app
def main():
    st.title("GenBank Sequence Downloader")

    # Get user input
    accession_numbers = st.text_area("Enter GenBank Accession Numbers (comma-separated):").split(',')

    # Button to trigger the download
    if st.button("Download Sequences"):
        if accession_numbers:
            df = fetch_sequences(accession_numbers)
            dfs = find_shared(df)
            zip_buffer = create_zip_file(dfs)

            # Encode the zip buffer as base64
            zip_base64 = base64.b64encode(zip_buffer.getvalue()).decode('utf-8')

            # Provide download link for the zip file
            st.success("Files created and zipped successfully!")
            st.markdown(f"**[Download Zip File](data:application/zip;base64,{zip_base64})**")
        else:
            st.warning("Please provide valid accession numbers.")

if __name__ == "__main__":
    main()
