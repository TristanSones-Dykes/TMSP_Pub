from bs4 import BeautifulSoup as bs
import mechanicalsoup
import pandas as pd
import os


# ------ Description ------
# This script contains functions to run phobius transmembrane
# prediction on a list of proteins in FASTA format by web scraping
# the form of the online tool

# url for phobius
url = "https://phobius.sbc.su.se/index.html"


# function to run phobius on a FASTA formatted list of proteins
# reads from file by default
def phobius(
    origin: str, isString: bool = False, fullSignal: bool = False
) -> pd.DataFrame:
    if isString:
        fileString = origin
    else:
        with open(origin, "r") as f:
            fileString = f.read()

    # instantiate browser, open url, select form
    browser = mechanicalsoup.StatefulBrowser(soup_config={"features": "html.parser"})
    browser.open(url)
    browser.select_form("form")

    # fill in form and check
    browser["format"] = "nog"
    if isString:
        browser["protseq"] = fileString
    else:
        file = open(origin, "rb")
        browser.form.set("protfile", file)

    response = browser.submit_selected()

    # parse response to get all data in <pre> tag
    soup = bs(response.text, "html.parser")
    out = soup.find("pre")

    if out is None:
        raise ValueError(
            "Could not find <pre> tag in response. The website structure might have changed or the request failed."
        )

    # right -> left:
    # splits output text into list of protein outputs
    # then splits each protein output into list of non-empty lines and collapses whitespace
    out_split = [
        [y.split() for y in x.split("\n") if y != ""] for x in out.text.split("//")
    ][:-1]

    assert len(out_split) == len(fileString.split(">")) - 1, (
        "number of proteins in output does not match number of proteins in input"
    )

    # split into predicted and unpredicted, and predicted into transmembrane and signal
    predicted = []
    unpredicted = []
    for protein in out_split:
        if "TRANSMEM" in [x[1] for x in protein]:
            predicted.append((protein, "TRANSMEM"))
        elif "SIGNAL" in [x[1] for x in protein]:
            predicted.append((protein, "SIGNAL"))
        else:
            unpredicted.append(protein)

    # populate dataframe with predicted transmembrane window sizes
    window_df = pd.DataFrame(columns=["ID", "start", "end", "type"])
    # setting datatypes
    window_df = window_df.astype({"ID": str, "start": int, "end": int, "type": str})

    for protein, window_type in predicted:
        # extract category column and protein ID
        col_2 = [x[1] for x in protein]
        proteinID = col_2[0]

        if window_type == "TRANSMEM":
            # grab transmembrane row
            window_row = col_2.index("TRANSMEM")
        elif window_type == "SIGNAL":
            # search region rows for "H-REGION."
            if not fullSignal:
                region_rows = [i for i, x in enumerate(col_2) if x == "REGION"]
                window_row = region_rows[1]
            else:
                # if fullSignal is true, grab signal row
                window_row = col_2.index("SIGNAL")

        start = int(protein[window_row][2])
        end = int(protein[window_row][3])

        # append to dataframe
        new_row = pd.DataFrame(
            [[proteinID, start, end, window_type]],
            columns=["ID", "start", "end", "type"],
        )
        window_df = pd.concat([window_df, new_row], ignore_index=True)

    for protein in unpredicted:
        # extract category column and protein ID
        col_2 = [x[1] for x in protein]
        proteinID = col_2[0]

        # append to dataframe
        new_row = pd.DataFrame(
            [[proteinID, 0, 0, "NONE"]], columns=["ID", "start", "end", "type"]
        )
        window_df = pd.concat([window_df, new_row], ignore_index=True)

    # delete temp file
    if isString:
        os.remove("temp.fasta")

    return window_df
