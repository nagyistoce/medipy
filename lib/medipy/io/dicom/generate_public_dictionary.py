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

def underscored_name(name):
    """ Return the name in the underscore_separated_name formalism
    """
    
    # Some names finish by a space, which seems incorrect
    name = "".join([x for x in name if x not in r"""!@#$%^&*(),;:.?\|{}[]+-="'/"""])
    name = name.strip()
    name = name.replace(" ", "_")
    name = name.lower()
    
    return name

def main():

    fd = urllib.urlopen(sys.argv[1])
    input = xml.dom.minidom.parseString(fd.read())
    fd.close()
    
    data_dictionary = {}
    
    for child in input.getElementsByTagName("dict")[0].childNodes :
        if child.nodeType != xml.dom.Node.ELEMENT_NODE or child.nodeName != "entry":
            continue
        
        tag = child.getAttribute("group") + child.getAttribute("element")
        
        if "x" not in tag :
            tag = int(tag, 16)
        
        vr = "/".join((child.getAttribute("vr") or "UN").split("_"))
        vm = child.getAttribute("vm") or "1"
        name = child.getAttribute("name")
        retired = (child.getAttribute("retired") == "true")
        
        data_dictionary[tag] = (vr, vm, name, retired, underscored_name(name))
    
    print "data_dictionary = {"
    
    for key in sorted(data_dictionary.keys()) :
        value = data_dictionary[key]
        if isinstance(key, int) :
            print "%#010x : %s,"%(key, value)
        else :
            print "\"%8s\" : %s,"%(key, value)
    
    print "}"
    
if __name__ == "__main__" :
    main()
