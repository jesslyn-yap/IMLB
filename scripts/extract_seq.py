import pandas as pd
from pyfaidx import Fasta
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor

_genome = None # Global variable for genome

def init_worker():
    """Initialize worker with global genome object"""
    global _genome
    _genome = Fasta("genome/hg19.fa")

def seq_single_parallel(args):
    """Function to obtain sequence for one gene"""

    chr_, bin_start, bin_end = args
    strand = _genome[chr_][bin_start:bin_end].seq
    strand = strand.upper().replace('R', 'N').replace('Y', 'N')

    return strand

def process_chunk(chunk, num_workers):
    """Function to obtain sequences in a single chunk in parallel"""

    args_iter = list(zip(chunk["Chromosome"], chunk["Start"], chunk["End"]))
    seqs = []
    
    with ProcessPoolExecutor(max_workers=num_workers, initializer=init_worker) as executor:
        for strand in executor.map(seq_single_parallel, args_iter):
            seqs.append(strand)

    res = chunk.copy()
    res['seq'] = seqs
    return res

def sequential_process_chunk(chunk):
    """Function to obtain sequences in a single chunk sequentially"""

    genome = Fasta("genome/hg19.fa")
    seqs = []
    
    for _, gene_row in chunk.iterrows():
        chr_ = gene_row['Chromosome']
        bin_start = gene_row['Start']
        bin_end = gene_row['End']

        strand = genome[chr_][bin_start:bin_end].seq
        strand = strand.upper().replace('R', 'N').replace('Y', 'N')
        seqs.append(strand)
    
    res = chunk.copy()
    res['seq'] = seqs
    return res

def parallel_seq_chunked(df: pd.DataFrame, num_workers=8, chunk_size=20000, use_parallel=True):
    """
    Function to obtain sequences in a df in chunks
    
    Args:
        df: Dataframe which contains 'Chromosome', 'Start', 'End' columns
        num_workers: Number of parallel processes
        chunk_size: Number of rows in df processed at once
        use_parallel: Whether to parallelize the process of obtaining the sequences
    
    Returns:
        df 
    """

    chunks = [df[i:i + chunk_size] for i in range(0, len(df), chunk_size)]
    processed_chunks = []
    
    for chunk in tqdm(chunks, desc="Processing chunks"):
        if use_parallel:
            try:
                processed_chunk = process_chunk(chunk, num_workers)
            except:
                processed_chunk = sequential_process_chunk(chunk)
        else:
            processed_chunk = sequential_process_chunk(chunk)
        
        processed_chunks.append(processed_chunk)
    
    return pd.concat(processed_chunks, ignore_index=True)