##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import logging
import sys

from tag import Tag

execfile(sys.argv[1])

name_dictionary = dict([(v[4], t) for t, v in data_dictionary.items() if isinstance(t, int)])

print "name_dictionary = {"
for key in sorted(name_dictionary.keys()) :
    value = name_dictionary[key]
    try :
        str(key)
    except UnicodeEncodeError :
        logging.warning("Skipping %s %s as it contains Unicode characters"%(key, value))
        continue
    print "\"%s\" : %#010x,"%(key, value)
print "}"
