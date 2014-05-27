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
    
def normalize_VR(inputVR):
    """ Rules :
        1- Type1 or Type2 => Type1/Type2
        2- See Note => UN
        3- Can't be null => UN
    """
    
    outputVR = inputVR.replace(" or ", "/")
    outputVR = outputVR.replace("See Note", "UN")
    
    if outputVR == "":
        return "UN"
    
    return outputVR
    
def normalize_VM(inputVM):
    """ Rules
        1- Can't be null => 1
    """
    
    outputVM = inputVM
    if outputVM == "":
        return "1"
        
    return outputVM
    
def normalize_name(inputName):
    """ Rules
        1- Delete non alpha char
    """
    
    outputName = re.split(ur"[\u200b]+", inputName, flags=re.UNICODE)
    
    outputName = "".join(outputName)
    return outputName
    
def normalize_Keyword(inputKeyword):
    """ Rules
        1- Delete non alpha char
        2- Add '_' between words
        3- Character => lower
    """
    
    # Delete space
    outputKeyword = re.split(r"\s", inputKeyword)
    
    temp = []
    # Delete non alpha char
    for word in outputKeyword:
        tempstr = re.split(r"\W", word)
        temp.append("".join(tempstr))
    
    # Add '_' to replace space
    outputKeyword = "_".join(temp)
    return outputKeyword.lower()
    
def get_td_value(tdNode):
    """ Return value for a given column
    
        tdNode doit etre un element xml au format
        <td>
            <para>
                <emphasis>Value</emphasis>
            </para>
        </td>
        Ou un element xml au format
        <td>
            <para>Value</para>
        </td>
    """
    
    # Read tag PARA and EMPHASIS
    para = [node for node in tdNode.getElementsByTagName("para")]
    emphasis = [node for node in para[0].getElementsByTagName("emphasis")]
    
    # if tag EMPHASIS is missing
    if len(emphasis) == 0 or not emphasis[0].firstChild:
        # if tag PARA is missing
        if para[0].childNodes:
            return para[0].firstChild.data
            
    # Else if tag EMPHASIS is not empty
    elif emphasis[0].childNodes:
        return emphasis[0].firstChild.data
        
    # Default value
    return ""
    
def get_tr_values(trNode):
    """ Return all column values for a given row
    """
    
    returnValues = {}
    
    tdList = [node for node in trNode.getElementsByTagName("td")]
    
    # Read all column values
    tag = get_td_value(tdList[0])
    name = get_td_value(tdList[1])
    keyword = get_td_value(tdList[1]) # ! Use name !
    vr = get_td_value(tdList[3])
    vm = get_td_value(tdList[4])
    retired = "RET" in tdList[5].toxml()
    
    # Formatting tag => OxAAAABBBB
    tag = tag[1:5] + tag[6:10]
    if "x" not in tag:
        tag = int(tag, 16)
        
    # Formatting VR
    vr = normalize_VR(vr)
    
    # Formatting VM
    vm = normalize_VM(vm)
    
    # Formatting Keyword
    keyword = normalize_Keyword(keyword)
    
    # Formatting Name
    name = normalize_name(name)
    
    # Return values
    returnValues[tag] = (vr, vm, name, retired, keyword)
        
    return returnValues

def process_table(table):
    dictionary = {}
    
    tbody = table.getElementsByTagName("tbody")[0]
    # Read each row (TR)
    for child in tbody.getElementsByTagName("tr"):
        # Get information from columns (TD)
        dictionary.update(get_tr_values(child))
    
    return dictionary

def main():

    fd = urllib.urlopen(sys.argv[1])
    input = xml.dom.minidom.parseString(fd.read())
    fd.close()
    
    data_dictionary = {}
    
    tables = input.getElementsByTagName("table")
    
    # Registry of DICOM Data Elements
    table = [node for node in tables 
             if node.getAttribute("xml:id") == "table_6-1"][0]
    data_dictionary.update(process_table(table))
    
    # Registry of DICOM File Meta Elements
    table = [node for node in tables 
             if node.getAttribute("xml:id") == "table_7-1"][0]
    data_dictionary.update(process_table(table))
    
    tbody = table.getElementsByTagName("tbody")[0]
    # Read each row (TR)
    for child in tbody.getElementsByTagName("tr"):
        # Get information from columns (TD)
        data_dictionary.update(get_tr_values(child))
    
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
