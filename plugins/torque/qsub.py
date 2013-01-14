import copy
import os
import shutil
import subprocess

import arguments

def submit(job_script, script_arguments, working_directory):
    """ Submit a single job to a Torque master node. The job identifier is 
        returned, unless an error occured. In that case, an exception is raised.
    """
    
    # Copy the job script to the working directory
    remote_job_script = os.path.join(working_directory, os.path.basename(job_script))
    shutil.copy(job_script, remote_job_script)
    
    # Run qsub from the working directory
    environment = "medipy_torque={0}".format(arguments.encode(script_arguments))
    command = ["qsub", remote_job_script, "-v", environment]
    process = subprocess.Popen(command, cwd=working_directory, 
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if process.returncode == 0 :
        job_id = stdout.strip()
        return job_id
    else :
        raise medipy.base.Exception(stderr.strip())

def submit_array(job_script, scripts_arguments_list, working_directory):
    """ Submit an array of jobs to a Torque master node. The job identifier is 
        returned, unless an error occured. In that case, an exception is raised.
    
        scripts_arguments_list must be a list of arguments, each item will be
        passed to one of the individual jobs.
    """
    
    # Copy the job script to the working directory
    remote_job_script = os.path.join(working_directory, os.path.basename(job_script))
    shutil.copy(job_script, remote_job_script)
    
    # Run qsub from the working directory
    environment = "medipy_torque={0}".format(arguments.encode(scripts_arguments_list))
    command = ["qsub", "-t", "0-{0}".format(len(scripts_arguments_list)-1),
               "-v", environment, remote_job_script]
    process = subprocess.Popen(command, cwd=working_directory, 
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if process.returncode == 0 :
        job_id = stdout.strip()
        return job_id
    else :
        raise medipy.base.Exception(stderr.strip())
