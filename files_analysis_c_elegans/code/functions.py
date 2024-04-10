import ribopy
from ribopy import Ribo
from Fasta import FastaFile
from ribopy.core.get_gadgets import get_region_boundaries, get_reference_names
import pandas as pd

def find_sequence(ribo_object, reference_file_path):
    transcript_np = ribo_object.transcript_names
    fasta = FastaFile(reference_file_path)
    
    fasta_dict = {e.header: e.sequence for e in fasta}
    sequence_dict = {
        transcript: fasta_dict[transcript] for transcript in transcript_np
    }
    return sequence_dict

def get_psite_offset(ribo_object, exp, mmin, mmax) :    
    df = (ribo_object.get_metagene("start", experiments = exp,\
                                   range_lower= mmin, range_upper= mmax,\
                                   sum_lengths = False,\
                                   sum_references = True))

    p_site = {}

    for index, row in df.iterrows():
        max_value_index = row.iloc[35:41].idxmax()
        offset = -1 * max_value_index + 1

        p_site[index[1]] = offset 
    return p_site

def find_cds_boundaries(bed_file):
    bed_df = pd.read_csv(bed_file, sep='\t', header=None, names=['Transcript', 'Start', 'End', 'Type', 'Score', 'Strand'])
    filtered_df = bed_df[bed_df['Type'] == 'CDS'].reset_index()
    return filtered_df