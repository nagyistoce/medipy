##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import re
import sys
import urllib
import xml.dom.minidom

def underscored_name(name):
    """ Return the name in the underscore_separated_name formalism
    """

    elements = re.split(r"[^\w]", name)
    result = "_".join([x for x in elements if x]).lower() 
    
    return result

def main():

    fd = urllib.urlopen(sys.argv[1])
    input = xml.dom.minidom.parseString(fd.read())
    fd.close()
    
    uid_dictionary = {}
    
    for table_element in input.getElementsByTagName("dict")[0].childNodes :
        if table_element.nodeType != xml.dom.Node.ELEMENT_NODE or table_element.nodeName != "table":
            continue
        
        if table_element.getAttribute("ref") != "Table A-1" :
            continue
        
        for uid_element in table_element.childNodes :
            if uid_element.nodeType != xml.dom.Node.ELEMENT_NODE or uid_element.nodeName != "uid":
                continue
            uid = uid_element.getAttribute("value")
            # Typo in the table ...
            uid = re.sub(r"[^0-9.]", "", uid)
            name = uid_element.getAttribute("name")
            type_ = uid_element.getAttribute("type")
            part = uid_element.getAttribute("part")
            retired = (uid_element.getAttribute("retired") == "true")
            
            uid_dictionary[uid] = (name, type_, part, retired, underscored_name(name)) 
    
    print "uid_dictionary = {"
    for key in sorted(uid_dictionary.keys()) :
        value = uid_dictionary[key]
        print "\"{0}\" : {1},".format(key, value)
    print "}"
    
    uid_name_dictionary = {}
    for key, value in uid_dictionary.items() :
        if re.match(r"[0-9a-zA-Z]", value[4]) :
            uid_name_dictionary[value[4]] = key
    
    print "uid_name_dictionary = {"
    for key in sorted(uid_name_dictionary.keys()) :
        value = uid_name_dictionary[key]
        print "\"{0}\" : \"{1}\",".format(key, value)
    print "}"
    
if __name__ == "__main__" :
    main()
