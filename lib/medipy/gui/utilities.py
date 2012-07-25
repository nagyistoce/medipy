##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import re

def underscore_to_camel_case(string):
    elements = string.split("_")
    
    result = ""
    for element in elements :
        result += element[0].upper() + element[1:]
    return result

def remove_access_letter_from_menu(string):
    """ Remove any ampersand that is not followed by another ampersand
    """
    
    result = re.split(re.compile("&(?!&)"), string)
    
    return reduce(lambda x,y : x+y, result, "") 