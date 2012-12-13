##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
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
    
    for dictionary in input.getElementsByTagName("dict") :
        for entry in dictionary.childNodes :
            if entry.nodeName != "entry" :
                continue
        
            tag = entry.getAttribute("group") + entry.getAttribute("element")
            
            if "x" not in tag :
                tag = int(tag, 16)
            
            vr = "/".join((entry.getAttribute("vr") or "UN").split("_"))
            vm = entry.getAttribute("vm") or "1"
            name = entry.getAttribute("name")
            retired = (entry.getAttribute("retired") == "true")
            
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
