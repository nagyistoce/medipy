##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import sys
import urllib
import xml.dom.minidom

from generate_public_dictionary import underscored_name

def main():

    fd = urllib.urlopen(sys.argv[1])
    input = xml.dom.minidom.parseString(fd.read())
    fd.close()
    
    data_dictionary = {}
    
    for child in input.getElementsByTagName("dict")[0].childNodes :
        if child.nodeType != xml.dom.Node.ELEMENT_NODE or child.nodeName != "entry":
            continue
        
        owner = child.getAttribute("owner")
        
        if owner not in data_dictionary :
            data_dictionary[owner] = {}
        
        owners_dictionary = data_dictionary[owner]
        
        tag = child.getAttribute("group") + child.getAttribute("element")
        
        if "x" not in tag :
            tag = int(tag, 16)
        
        vr = "/".join((child.getAttribute("vr") or "UN").split("_"))
        vm = child.getAttribute("vm") or "1"
        name = child.getAttribute("name")
        retired = (child.getAttribute("retired") == "true")
        
        owners_dictionary[tag] = (vr, vm, name, retired, underscored_name(name))
    
    print "private_dictionaries = {"
    
    for owner in sorted(data_dictionary.keys()) :
        print "    {0} : {{".format(repr(owner))
        owners_dictionary = data_dictionary[owner]
        for key in sorted(owners_dictionary.keys()) :
            value = owners_dictionary[key]
            if isinstance(key, int) :
                print "        %#010x : %s,"%(key, value)
            else :
                print "        \"%8s\" : %s,"%(key, value)
        print "    },"
    print "}"
    
if __name__ == "__main__" :
    main()