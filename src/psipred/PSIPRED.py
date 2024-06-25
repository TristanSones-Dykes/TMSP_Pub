import os
import requests
import json
import time
import pandas as pd
from Bio import SeqIO

# define urls for submission and download
submission_uri = 'http://bioinf.cs.ucl.ac.uk/psipred/api/submission'
download_uri = 'http://bioinf.cs.ucl.ac.uk/psipred/api'

# split input fasta into individual sequence files of length 60
if not os.path.exists('temp'):
    os.makedirs('temp')

if not os.path.exists('results'):
    os.makedirs('results')

for seq_record in SeqIO.parse("S_Cerevisiae.fa", "fasta"):
    seq_id = seq_record.id
    if len(str(seq_record.seq)) < 61:
        continue
    seq = str(seq_record.seq)[:61]
    with open('temp/'+seq_id+'.fa', 'w') as f:
        f.write('>'+seq_id+'\n'+seq)

# remove all files from temp that have been processed
for seq_file in os.listdir('temp'):
    if os.path.exists('results/'+seq_file+'.horiz'):
        with open('results/'+seq_file+'.horiz', 'r') as f:
            lines = f.read().splitlines()
            first = lines[0]
        if first[0] != ">":
            print('Removing: '+seq_file)
            os.remove('temp/'+seq_file)
        else:
            print('Keeping: '+seq_file+" didn't complete")

# submit each sequence to the server
for seq_file in os.listdir('temp'):
    # setup the request
    payload = {'input_data': (seq_file, open('temp/'+seq_file, 'rb'))}
    data = {'job': 'psipred',
            'submission_name': seq_file,
            'email': "tristansonesdykes@outlook.com"
            }
    
    # send the request
    print("\nSending Request")
    r = requests.post(submission_uri+".json", data=data, files=payload)
    response_data = json.loads(r.text)
    print(response_data)
    time.sleep(10)

    # Please do no poll the server more than once every 30 seconds. If you DOS the server we will block your IP range
    while True:
        print("Polling result for:"+response_data['UUID'])
        result_uri = submission_uri+"/"+response_data['UUID']
        try:
            r = requests.get(result_uri, headers={"Accept":"application/json"}, timeout = 30)
        except requests.exceptions.Timeout:
            print("Timeout occurred")
            continue
        result_data = json.loads(r.text)
        if "Complete" in result_data["state"]:
            break
        else:
            time.sleep(30)

    horiz_path = result_data['submissions'][0]['results'][-1]['data_path']
    # http://bioinf.cs.ucl.ac.uk/psipred/api/submissions/[FILE_NAME]
    # Download the file
    r = requests.get(download_uri + horiz_path)
    with open('results/'+seq_file+'.horiz', 'w') as f:
        f.write(r.text)
    time.sleep(1)

# go through results file and create sequence CSV
def psipred_df(conf_lower: int, length_lower: int, count = False) -> pd.DataFrame:
    if count:
        from collections import Counter

    row_list = []
    res_files = os.listdir('results')
    for res_file in res_files:
        # read file
        with open('results/'+res_file, 'r') as f:
            lines = f.read().splitlines()

        # take predictions and subset using confidence
        conf_1, conf_2 = lines[2][6:], lines[11][6:]
        pred_1, pred_2 = lines[3][6:], lines[12][6:]
        conf = conf_1 + conf_2
        pred = pred_1 + pred_2

        conf_pred = []
        for (confidence, prediction) in zip(conf, pred):
            if int(confidence) >= conf_lower:
                conf_pred.append(prediction)

        length = 0
        if count:
            # return count of helix region predicted AAs
            pred_counts = Counter(conf_pred)
            length = pred_counts.get("H", 0)
        else:
            # add first helix of length >= length_lower
            i = 0

            # search entire sequence
            while i < len(conf_pred) - 1:
                # find next helix
                while conf_pred[i] != "H" and i < len(conf_pred) - 1:
                    i += 1
                # find length of helix
                while conf_pred[i] == "H" and i < len(conf_pred) - 1:
                    length += 1
                    i += 1
                # check if length is long enough
                if length < length_lower:
                    length = 0
                else:
                    break
        
        if length > 0:
            row_list.append({'seqid': res_file[:-9], 'window_length': length})

    seq_df = pd.DataFrame(row_list)
    return seq_df