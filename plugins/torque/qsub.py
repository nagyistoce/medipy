import copy
import os
import shutil
import subprocess

import arguments

def submit(job_script, script_arguments, working_directory):
    """ Submit a single job to a Torque master node.
    """
    
    # Copy the job script to the working directory
    remote_job_script = os.path.join(working_directory, os.path.basename(job_script))
    shutil.copy(job_script, remote_job_script)
    
    # Run qsub from the working directory
    environment = "medipy_torque={0}".format(arguments.encode(script_arguments))
    subprocess.call(["qsub", remote_job_script, "-v", environment], cwd=working_directory)

def submit_array(job_script, scripts_arguments_list, working_directory):
    """ Submit an array of jobs to a Torque master node.
    
        scripts_arguments_list must be a list of arguments, each item will be
        passed to one of the individual jobs.
    """
    
    # Copy the job script to the working directory
    remote_job_script = os.path.join(working_directory, os.path.basename(job_script))
    shutil.copy(job_script, remote_job_script)
    
    # Run qsub from the working directory
    environment = "medipy_torque={0}".format(arguments.encode(scripts_arguments_list))
    subprocess.call([
         "qsub", 
         "-t", "0-{0}".format(len(scripts_arguments_list)-1),
         "-v", environment, 
         remote_job_script
    ], cwd=working_directory)
