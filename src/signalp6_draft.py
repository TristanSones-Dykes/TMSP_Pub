

import pandas as pd
import numpy as np
import re
from Bio import SeqIO # read and write fasta files, for the sequences 

import biolib
signalp_6 = biolib.load('DTU/SignalP-6')
print(signalp_6.cli(args='--help')) 


from pathlib import Path
from typing import Union, List
from Bio import SeqIO


def split_fasta(filename: Union[str, Path], nperfile: int) -> List[Path]:
    """
    Split a FASTA file into multiple FASTA files, each containing up to nperfile sequences.

    Parameters
    ----------
    filename : str or pathlib.Path
        Path to the input FASTA file.
    nperfile : int
        Number of sequences per output file (must be >= 1).

    Returns
    -------
    list[pathlib.Path]
        Paths to the created output FASTA files.

    Notes
    -----
    - Output files are named by appending a numeric suffix to the input stem, e.g.:
      input.fasta -> input_001.fasta, input_002.fasta, ...
    - Each input record is written exactly once to exactly one output file.
    
    Example Usage
    -----   
    outs = split_fasta("input.fasta", nperfile=100)
    print(f"Wrote {len(outs)} files:")
    for p in outs:
        print(" -", p)
    """
    if nperfile < 1:
        raise ValueError("nperfile must be >= 1")

    in_path = Path(filename)
    if not in_path.exists():
        raise FileNotFoundError(in_path)

    out_paths: List[Path] = []
    records_buffer = []
    file_index = 1

    def write_chunk(chunk, idx: int):
        out_path = in_path.with_name(f"{in_path.stem}_{idx:03d}{in_path.suffix or '.fasta'}")
        with out_path.open("w") as handle:
            SeqIO.write(chunk, handle, "fasta")
        out_paths.append(out_path)

    with in_path.open("r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            records_buffer.append(record)
            if len(records_buffer) == nperfile:
                write_chunk(records_buffer, file_index)
                records_buffer = []
                file_index += 1

    # Write any remaining sequences
    if records_buffer:
        write_chunk(records_buffer, file_index)
    
    return out_paths

SC_outs = split_fasta("data/Proteins/full/S_Cerevisiae.fasta", nperfile=500)
print(f"Wrote {len(SC_outs)} files:")
for p in SC_outs:
    print(" -", p)

## the next lines of code did not work
job_test = signalp_6.run(args = "--fastafile" + str(SC_outs[0]) + "--output_dir results/signalp6/S_Cerevisiae_001 --format txt --organism eukarya --mode fast")
job_test2 = signalp_6.cli(args = "--fastafile" + str(SC_outs[1]) + "--output_dir results/signalp6/S_Cerevisiae_001 --format txt --organism eukarya --mode fast")


# run_signalp_batch =