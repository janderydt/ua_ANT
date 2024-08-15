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
                            "Restart": int, "InverseIterationsDone": int, "Velocity": int, "startGeometry": int, "startC": int, "startAglen": int})
        data = data.astype({"gsC": float, "gsA": float, "gaC": float, "gaA": float, "SpinupYearsDone": float, 
                            "m": float, "muk": float, "priorC": float, "n": float, "priorAGlen": float})
        data = data.astype({"Domain": str, "InverseIterations": str, "SpinupYears": str, "InvertFor": str, "GradientCalc": str, 
                            "Measurements": str, "SlidingLaw": str,
                            "startMesh": str, "Comments": str})
    elif "Diagnostic" in runtype:
        data = data.astype({"pgid": int, "ExpID": int, "Submitted": int, "Running": int, "Error": int, "Finished": int,
                            "Restart": int, "InverseA": int, "InverseCycleA":int, "InverseC": int, "InverseCycleC": int, 
                            "Calv": int, "ISthick": int, "InverseCycleIS": int, "GIthick": int, "InverseCycleGI": int})
        data = data.astype({"Domain": str, "BaseMesh": str, "Comments": str})
    elif "Transient" in runtype:
        data = data.astype({"pgid": int, "ExpID": int, "Submitted": int, "Running": int, "Error": int, "Finished": int,
                            "Restart": int, "AdaptMesh": int, "Geometry": int, "InverseA": int, "InverseCycleA":int, 
                            "InverseC": int, "InverseCycleC": int, "OceForcing": int})
        data = data.astype({"PICO_C1": float, "PICO_gam": float, "LQ_gam": float})
        data = data.astype({"Domain": str, "SMB": str, "BasalMelt": str})
    else:
        log.info('   ...Unknown run type')

    
    data['SubmissionTime'] = pd.to_datetime(data['SubmissionTime'],errors = 'coerce')
    data['ErrorTime'] = pd.to_datetime(data['ErrorTime'],errors='coerce')
    data['FinishedTime'] = pd.to_datetime(data['FinishedTime'],errors='coerce')

    log.info('   ...Successfully read '+table)

    return data

# Convert a list of strings to a single string with elements separated by the given separator character.
def list_with_separator (A, sep):
    s = ''
    for elm in A:
        s += elm + sep
    # Remove the last character
    return s[:-1]  

def save_runinfo (df,file_name):

    # Save Run Table

    df.to_csv(file_name, encoding='utf-8', index=False)

    log.info('   ...Successfully saved '+file_name)

# Submit batch script on local machine
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

# Submit the given SBATCH script on ARCHER2 and return the PBS job ID.
# Optional keyword arguments:
# options: Options object
# input_var: a list of variable definitions to pass with --export=, eg 'UA_CONFIG=<path>'
# afterok: a list of SBATCH job IDs of previously submitted jobs. If it is defined, this job will stay on hold until the given jobs successfully complete.
def submit_job (budget_code, Nnodes, sbatch_script, input_var=None, afterok=None):

    # Construct sbatch call line by line.
    command = 'sbatch'
    # Specify budget
    command += ' -A ' + budget_code
    command += ' -N ' + Nnodes
    if input_var is not None:
        # Add variable definitions
        command += ' --export=ALL,'
        command += list_with_separator(input_var,',')
    if afterok is not None:
        command += ' -W depend=afterok:'
        command += list_with_separator(afterok,':')
    # Specify script
    command += ' ' + sbatch_script
    
    # Call the command and capture the output
    sbatch_id = subprocess.check_output(command, shell=True, text=True)

def active_jobs (command):

    # Check for active jobs

    # Call the command and capture the output
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    stdout, stderr = p.communicate()

    processes = stdout.decode('ascii').split()
    
    return(processes)





