import ribopy
from ribopy import Ribo
from functions import get_psite_offset, find_sequence, find_cds
import numpy as np
import pandas as pd
import multiprocessing
from collections import Counter, defaultdict
import time

# Initialize variables
ribo_path   = '/Users/reikotachibana/Documents/files_analysis_c_elegans/ribo_files/Rep_5.ribo'
reference_path = '/Users/reikotachibana/Documents/files_analysis_c_elegans/celegans_reference/appris_celegans_v1_selected_new.fa'
exp= "WT_1cell_D"
min_len = 26
max_len = 40
ribo_object = Ribo(ribo_path)
sequence = find_sequence(ribo_object, reference_path)
alias = False
transcript_list = ribo_object.transcript_names
# transcript_list = ['Y110A7A.10.1|cdna|chromosome:WBcel235:I:5107833:5110183:1|gene:WBGene00000001.1|gene_biotype:protein_coding|transcript_biotype:protein_coding|gene_symbol:aap-1',
# 'F27C8.1.2|cdna|chromosome:WBcel235:IV:9598986:9601695:-1|gene:WBGene00000002.1|gene_biotype:protein_coding|transcript_biotype:protein_coding|gene_symbol:aat-1']
offset = get_psite_offset(ribo_object, exp, min_len, max_len)

# Dictionary of CDS ranges {transcript: (start, stop)}
cds_range = {} 
for transcript in transcript_list:
    cds_range[transcript] = find_cds(sequence[transcript])
    start, stop = cds_range[transcript]

# Finds coverage using p-site offset for a transcript
def find_coverage(ribo_object, transcript, cds_range, offset, exp, min_len, max_len, alias):
    start, stop = cds_range[transcript]
    coverages = [
        ribo_object.get_coverage(experiment=exp, range_lower=i, range_upper=i, alias=alias)
        [transcript][start - offset[i] : stop - offset[i]]
        for i in range(min_len, max_len + 1)
    ]

    if any(coverage.size == 0 for coverage in coverages):
        # If any coverage is empty, return an array of zeros
        return np.zeros_like(coverages[0])
    else:
        # Otherwise, sum the coverages
        coverage = sum(coverages, np.zeros_like(coverages[0]))
        return coverage


# Finds codon occupancy of each footprint for a transcript
def find_codon_occupancy(coverage, cds_range, sequence, transcript):
    start, stop = cds_range[transcript]
    
    # Initialize a dictionary to store the number of reads for each codon
    codon_reads = defaultdict(int)
    
    # Iterate over the range of the CDS in steps of 3 (since codons are of length 3)
    for i in range(start, stop + 1, 3):
        codon = sequence[transcript][i:i+3] # Get the codon from the sequence
        total_reads = sum(coverage[i-start:i-start+3])  # Sum up the coverage values for the codon
        codon_reads[codon] += total_reads
    
    return codon_reads

def process_transcript(transcript):
    coverage = find_coverage(ribo_object, transcript, cds_range, offset, exp, min_len, max_len, alias)
    
    if sum(coverage) > 0:
        codon_occ = find_codon_occupancy(coverage, cds_range, sequence, transcript)
        return codon_occ
    else:
        return Counter()


if __name__ == '__main__':
    start_time = time.time()

    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    results = []
    for transcript in transcript_list:
        start, stop = cds_range[transcript]
        if start is None or stop is None:
            continue  # Skip this transcript if start or stop is None

        result = pool.apply_async(process_transcript, args=(transcript,))
        results.append(result)

    pool.close()
    pool.join()

    codon_occupancy = Counter()
    for result in results:
        codon_occupancy.update(result.get())

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Time taken: {elapsed_time} seconds")
    
        # Convert the dictionary to a DataFrame
    df = pd.DataFrame(sorted(codon_occupancy.items()), columns=['Codon', 'Reads'])

    # Save the DataFrame as a CSV file
    df.to_csv('codon_reads.csv', index=False)

    # Print the DataFrame
    print(df)

