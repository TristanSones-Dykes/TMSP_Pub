from typing import Optional

import pandas as pd
import numpy as np
import re


def extract_deepTMHMM(origin: str, max_regions: Optional[int] = None) -> pd.DataFrame:
    with open(origin, "r") as f:
        fileString = f.read()

    # split into list of proteins, remove first line of file
    proteins = fileString.split("//")

    # split each protein into list of lines, remove first and last lines of each protein
    proteins = [
        [re.sub("\s+", " ", line).strip() for line in protein.split("\n")[1:-1]]  # type: ignore
        for protein in proteins
    ]

    # extract regions for each protein and filter out proteins with more than max_regions regions
    # and proteins with no regions
    protein_dat = [protein[2:] for protein in proteins]
    col_2 = [[line.split(" ")[1] for line in protein] for protein in protein_dat]
    region_count = [
        (protein.count("TMhelix"), protein.count("signal")) for protein in col_2
    ]

    match max_regions:
        case None:
            max_regions = max([max(x) for x in region_count])
        case int():
            pass
        case _:
            raise TypeError("max_regions must be an integer or None")

    tm = [i for i, x in enumerate(region_count) if 0 < x[0] <= max_regions]
    sig = [i for i, x in enumerate(region_count) if 0 < x[1] <= max_regions]
    tm_sig = list(np.intersect1d(tm, sig))

    # create dataframe for predicted transmembrane window sizes
    window_df = pd.DataFrame(
        columns=["seqid", "start", "end", "window_type", "protein_type"]
    )
    window_df = window_df.astype(
        {
            "seqid": str,
            "start": int,
            "end": int,
            "window_type": str,
            "protein_type": str,
        }
    )  # setting datatypes

    # grab protein data with tm, sig, and both regions
    tm_proteins, sig_proteins, tm_sig_proteins = (
        [protein_dat[i] for i in tm if i not in tm_sig],
        [protein_dat[i] for i in sig if i not in tm_sig],
        [protein_dat[i] for i in tm_sig],
    )

    for protein in sig_proteins:
        sig_row = protein[0].split(" ")

        seqid = sig_row[0]
        start = int(sig_row[2])
        end = int(sig_row[3])
        window_type = "SP"
        protein_type = "SP"

        # append to dataframe
        new_row = pd.DataFrame(
            [[seqid, start, end, window_type, protein_type]],
            columns=["seqid", "start", "end", "window_type", "protein_type"],
        )
        window_df = pd.concat([window_df, new_row], ignore_index=True)

    for protein in tm_proteins:
        col_2 = [line.split(" ")[1] for line in protein]
        tm_rows = [i for i, x in enumerate(col_2) if x == "TMhelix"]

        for tm_row in tm_rows:
            tm_row = protein[tm_row].split(" ")

            seqid = tm_row[0]
            start = int(tm_row[2])
            end = int(tm_row[3])
            window_type = "TM"
            protein_type = "TM"

            # append to dataframe
            new_row = pd.DataFrame(
                [[seqid, start, end, window_type, protein_type]],
                columns=["seqid", "start", "end", "window_type", "protein_type"],
            )
            window_df = pd.concat([window_df, new_row], ignore_index=True)

    for protein in tm_sig_proteins:
        col_2 = [line.split(" ")[1] for line in protein]
        tm_rows = [i for i, x in enumerate(col_2) if x == "TMhelix"]
        sig_row = protein[0].split(" ")

        for tm_row in tm_rows:
            tm_row = protein[tm_row].split(" ")

            seqid = tm_row[0]
            start = int(tm_row[2])
            end = int(tm_row[3])
            window_type = "TM"
            protein_type = "SP+TM"

            # append to dataframe
            new_row = pd.DataFrame(
                [[seqid, start, end, window_type, protein_type]],
                columns=["seqid", "start", "end", "window_type", "protein_type"],
            )
            window_df = pd.concat([window_df, new_row], ignore_index=True)

        seqid = sig_row[0]
        start = int(sig_row[2])
        end = int(sig_row[3])
        window_type = "SP"
        protein_type = "SP+TM"

        # append to dataframe
        new_row = pd.DataFrame(
            [[seqid, start, end, window_type, protein_type]],
            columns=["seqid", "start", "end", "window_type", "protein_type"],
        )
        window_df = pd.concat([window_df, new_row], ignore_index=True)

    return window_df
