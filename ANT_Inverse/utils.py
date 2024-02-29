# Lower level

import subprocess
from io import StringIO
import csv
import pandas as pd
import datetime
import os
import logging

log = logging.getLogger(__name__)

def read_runinfo (table):

    with open(table, newline='') as f:
        reader = csv.reader(f)
        data = list(reader)
    fields = [row for row in data[0]]
    data = pd.DataFrame(data[1:],columns=fields)
    data = data.astype({"pgid": int, "ExpID": int, "Submitted": int, "Running": int, "Error": int, "Finished": int,
                        "Restart": int, "Iterations": int, "IterationsDone": int, "gsC": int, "gsA": int,
                        "gaC": int, "gaA": int, "Velocity": int, "Geometry": int, "m": int, "Cstart": int, "n": int, "Aglenstart": int})
    data = data.astype({"InvertFor": str, "GradientCalc": str, "Measurements": str, "SlidingLaw": str,
                        "Mesh": str, "Comments": str})
    data['SubmissionTime'] = pd.to_datetime(data['SubmissionTime'],errors = 'coerce')
    data['ErrorTime'] = pd.to_datetime(data['ErrorTime'],errors='coerce')
    data['FinishedTime'] = pd.to_datetime(data['FinishedTime'],errors='coerce')

    log.info('   ...Successfully read '+table)

    return data

def save_runinfo (df,file_name):

    df.to_csv(file_name, encoding='utf-8', index=False)

    log.info('   ...Successfully saved '+file_name)

def submit_batch (options,command):

    # Call the command and capture the output
    p = subprocess.Popen(command, shell=True, preexec_fn=os.setsid)
    pgid = os.getpgid(p.pid)
    return_code = p.wait()

    if not return_code == 0:
        raise RuntimeError(
            'The process call "{}" returned with code {}.'
            .format(command, return_code))
    else:
        log.info('   ...Submitted new job with pgid '+str(pgid))
        return pgid


def active_jobs (command):

    # Call the command and capture the output
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    stdout, stderr = p.communicate()

    processes = stdout.decode('ascii').split()
    
    return(processes)





