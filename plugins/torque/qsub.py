import copy
import os
import shutil
import subprocess

import arguments

def submit(job_script, script_arguments, working_directory, 
           mail=None, mail_on_begin=False, mail_on_abort=True, mail_on_terminate=False):
    """ Submit a single job to a Torque master node. The job identifier is 
        returned, unless an error occured. In that case, an exception is raised.
        
        ``mail`` can be either an e-mail address or a list of e-mail address.
        The execution server will send mail to these addresses about the job.
        
        The three ``mail_on_...`` options defines the set of conditions under 
        which the execution server will send a mail message about the job:
        
            * ``mail_on_begin``: mail is sent when the job begins execution.
            * ``mail_on_abort``: mail is sent when the job is aborted by the batch system.
            * ``mail_on_terminate``: mail is sent when the job terminates.
    """
    
    # Copy the job script to the working directory
    remote_job_script = os.path.join(working_directory, os.path.basename(job_script))
    shutil.copy(job_script, remote_job_script)
    
    # Run qsub from the working directory
    environment = "medipy_torque={0}".format(arguments.encode(script_arguments))
    
    command = _get_command_line(remote_job_script, environment, None,
                                mail, mail_on_begin, mail_on_abort, mail_on_terminate)
    
    process = subprocess.Popen(command, cwd=working_directory, 
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if process.returncode == 0 :
        job_id = stdout.strip()
        return job_id
    else :
        raise medipy.base.Exception(stderr.strip())

def submit_array(job_script, scripts_arguments_list, working_directory, 
                 mail=None, mail_on_begin=False, mail_on_abort=True, mail_on_terminate=False):
    """ Submit an array of jobs to a Torque master node. The job identifier is 
        returned, unless an error occured. In that case, an exception is raised.
    
        scripts_arguments_list must be a list of arguments, each item will be
        passed to one of the individual jobs.
        
        See :func:`submit` for the description of mail options
    """
    
    # Copy the job script to the working directory
    remote_job_script = os.path.join(working_directory, os.path.basename(job_script))
    shutil.copy(job_script, remote_job_script)
    
    # Run qsub from the working directory
    environment = "medipy_torque={0}".format(arguments.encode(scripts_arguments_list))
    
    command = _get_command_line(remote_job_script, environment, 
                                (0, len(scripts_arguments_list)-1),
                                mail, mail_on_begin, mail_on_abort, mail_on_terminate)
    
    process = subprocess.Popen(command, cwd=working_directory, 
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if process.returncode == 0 :
        job_id = stdout.strip()
        return job_id
    else :
        raise medipy.base.Exception(stderr.strip())

def _get_command_line(script, environment, array,
                      mail, mail_on_begin, mail_on_abort, mail_on_terminate):
    """ Return the qsub command line for the given options.
    """
    
    command = ["qsub", script, "-v", environment]
    if array :
        command.extend(["-t", "{0}-{1}".format(*array)])
    command.extend(_parse_mail_options(mail, mail_on_begin, mail_on_abort, mail_on_terminate))
    
    return command

def _parse_mail_options(mail, mail_on_begin, mail_on_abort, mail_on_terminate):
    """ Generate the qsub options for sending email.
    """
    
    options = []
    if mail :
        if isinstance(mail, (list, tuple)) :
            mail = ",".join(mail)
        options.extend(["-M", mail])
        if True in [mail_on_begin, mail_on_abort, mail_on_terminate] :
            value = ""
            if mail_on_begin : 
                value += "b"
            if mail_on_abort : 
                value += "a"
            if mail_on_terminate : 
                value += "e"
            options.extend(["-m", value])
    
    return options
