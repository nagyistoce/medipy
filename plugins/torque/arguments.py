import base64
import cPickle
import os

import medipy.base

def get():
    """ Recover the decoded arguments from the process environment
    """
    
    if "medipy_torque" not in os.environ :
        raise medipy.base.Exception("Variable \"medipy_torque\" does not exist in environment")
    
    arguments = decode(os.environ["medipy_torque"])
    if "PBS_ARRAYID" in os.environ :
        # An array job was submitted, get the correct arguments based on PBS_ARRAY_ID
        arguments = arguments[int(os.environ["PBS_ARRAYID"])]
    
    return arguments
    
def decode(arguments):
    """ Decode the arguments from their ASCII form.
    """
    
    try :
        result = cPickle.loads(base64.b64decode(arguments))
    except :
        raise medipy.base.Exception("Arguments are not properly formatted")
    
    return result

def encode(arguments):
    """ Encode the arguments to ASCII form.
    """
    
    return base64.b64encode(cPickle.dumps(arguments))
