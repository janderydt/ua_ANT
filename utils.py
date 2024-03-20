# Lower level

import subprocess
from io import StringIO
import csv
import pandas as pd
import datetime
import os
import logging

log = logging.getLogger(__name__)

def read_runinfo (table,runtype):

    # Read details from the Run Table
    with open(table, newline='') as f:
        reader = csv.reader(f)
        data = list(reader)
    fields = [row for row in data[0]]
    data = pd.DataFrame(data[1:],columns=fields)

    if "Inverse" in runtype:
        data = data.astype({"pgid": int, "ExpID": int, "Submitted": int, "Running": int, "Error": int, "Finished": int,
                            "Restart": int, "InverseIterationsDone": int, "SpinupYearsDone": float, "gsC": int, "gsA": int,
                            "gaC": int, "gaA": int, "Velocity": int, "startGeometry": int, "m": int, "startC": int, "n": int, "startAglen": int})
        data = data.astype({"InverseIterations": str, "SpinupYears": str, "InvertFor": str, "GradientCalc": str, 
                            "Measurements": str, "SlidingLaw": str,
                            "startMesh": str, "Comments": str})
    elif "Diagnostic" in runtype:
        data = data.astype({"pgid": int, "ExpID": int, "Submitted": int, "Running": int, "Error": int, "Finished": int,
                            "Restart": int, "InverseA": int, "InverseAFill": int, "InverseCycleA":int , 
                            "InverseC": int, "InverseCFill": int, "InverseCycleC": int,
                            "Calv": int, "ISthick": int, "InverseCycleIS": int, "GIthick": int, "InverseCycleGI": int})
        data = data.astype({"BaseMesh": str, "Comments": str})
    else:
        log.info('   ...Unknown run type')

    
    data['SubmissionTime'] = pd.to_datetime(data['SubmissionTime'],errors = 'coerce')
    data['ErrorTime'] = pd.to_datetime(data['ErrorTime'],errors='coerce')
    data['FinishedTime'] = pd.to_datetime(data['FinishedTime'],errors='coerce')

    log.info('   ...Successfully read '+table)

    return data

def save_runinfo (df,file_name):

    # Save Run Table

    df.to_csv(file_name, encoding='utf-8', index=False)

    log.info('   ...Successfully saved '+file_name)

def submit_batch (options,command):

    # Execute new batch command

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

    # Check for active jobs

    # Call the command and capture the output
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    stdout, stderr = p.communicate()

    processes = stdout.decode('ascii').split()
    
    return(processes)





