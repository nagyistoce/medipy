# encoding: utf-8

##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import argparse
import ConfigParser
import csv
import re
import StringIO
import sys
import urllib
import xml.dom.minidom


def normalize_VR(inputVR, tag):
    """ Rules :
        1- Type1 or Type2 => Type1/Type2
        2- See Note => UN
        3- Can't be null => UN
    """
    
    all_vr = ["AE", "AS", "AT", "CS", "DA", "DS", "DT", "FL", "FD", "IS", "LO",
              "LT", "OB", "OD", "OF", "OW", "PN", "SH", "SL", "SQ", "SS", "ST",
              "TM", "UI", "UL", "UN", "US", "UT",
              "OB_or_OW", "US_or_OW", "US_or_SS", "US_or_SS_or_OW",
              None]
    
    if tag in [0xfffee000, 0xfffee00d, 0xfffee0dd]:
        # Item, Item Delimitation Item, Sequence Delimitation Item
        outputVR = None
    elif tag in [0x00189445, 0x00280020]:
        # See PS 3.6, note 6.3
        # For some Data Elements, no Name or Keyword or VR or VM is specified; 
        # these are "placeholders" that are not assigned but will not be reused.
        outputVR = None
    else:
        outputVR = inputVR.replace(" ", "_")
    
    if outputVR not in all_vr:
        raise Exception("Unknown VR for tag {0}: {1!r}".format(
            "{0:08x}".format(tag) if isinstance(tag, int) else tag, inputVR))
    
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
        2- Replace µ by u
        2- Add '_' between words
        3- Character => lower
    """
    
    outputKeyword = normalize_name(inputKeyword)
    outputKeyword = outputKeyword.replace(u"µ", "u")
    
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
    vr = normalize_VR(vr, tag)
    
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

def read_configuration_file(ifilepath):
    """ Parse Configuration file
    """
    
    Config = ConfigParser.ConfigParser()
    Config.read(ifilepath)
    
    return Config
     
def parse_gdcm(iurl, idictType, imainID, iName=None):
    """ Parse private dictionary from GDCM website
    """
    
    # open url and parse text
    fd = urllib.urlopen(iurl)
    input = xml.dom.minidom.parseString(fd.read())
    fd.close()
    
    # return array
    data_dictionary = {}
    
    # process "entry" nodes
    for child in input.getElementsByTagName("dict")[0].childNodes :
        if child.nodeType != xml.dom.Node.ELEMENT_NODE or child.nodeName != "entry":
            continue
        
        # Warning: getAttribute return an unicode string
        owner = child.getAttribute("owner")
        owner = owner.encode("ascii", "ignore")
        
        # add data to dictionary
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
        
        owners_dictionary[tag] = (vr, vm, name, retired, normalize_Keyword(name))
    
    return data_dictionary
    
def parse_dcmtk(iurl, idictType, imainID, iName=None):
    """ Parse private dictionary from DCMTK website
    """
    
    # open url
    fd = urllib.urlopen(iurl)
    input = StringIO.StringIO(fd.read())
    fd.close()
    
    # return array
    data_dictionary = {}
    
    # parse CSV text
    reader = csv.reader(input, delimiter='\t')
    
    # process row
    for row in reader:
        # ignore empty row or comment
        if len(row) == 0 or row[0] == "" or row[0][0] == '#':
            continue
        
        # parse first element ('Group', '"Private dictionary name"', 'Element')
        firstInfo = row[0].replace('(', '')
        firstInfo = firstInfo.replace(')', '')
        firstInfo = firstInfo.replace('"', '')
        firstInfo = firstInfo.replace('"', '')
        firstInfo = firstInfo.split(',')
        
        # get 'Private dictionary name'
        owner = firstInfo[1]
        
        # add data to dictionary
        if owner not in data_dictionary :
            data_dictionary[owner] = {}
            
        owners_dictionary = data_dictionary[owner]
        
        # get 'GroupxxElement': 
        # example: 0023xx01 with Group=0023 and Element=01
        tag = firstInfo[0] + "xx" + firstInfo[2]
        
        vr = row[1]
        vm = row[3]
        name = row[2]
        
        owners_dictionary[tag] = (vr, vm, name, False, name)
    
    return data_dictionary
    
def parse_dicom(iurl, idictType, imainID, iName=None):
    """ Parse private dictionary from FLI-IAM project
    """
    
    # return array
    data_dictionary = {}
    
    # input cannot be null
    if not imainID or imainID == "":
        return data_dictionary
    
    # open url and parse text
    fd = urllib.urlopen(iurl)
    input = xml.dom.minidom.parseString(fd.read())
    fd.close()
    
    tables = input.getElementsByTagName("table")
    
    # Registry of DICOM Data Elements
    table = [node for node in tables 
             if node.getAttribute("xml:id") == imainID][0]
    data_dictionary.update(process_table(table))
    
    # Format dictionary
    if idictType == "private" and iName:
        data_dictionary = { iName: data_dictionary }
        
    return data_dictionary
    
def print_private_dictionary(idictionary):
    """ Display Dictionary
        
        private_dictionary = {
            Private dictionary name 1 : {
                GroupxxElement1 : (VR, VM, Name, Retired, KeyWord),
                GroupxxElement2 : (VR, VM, Name, Retired, KeyWord)
            },
            Private dictionary name 2 : {
                GroupxxElement1 : (VR, VM, Name, Retired, KeyWord),
                GroupxxElement2 : (VR, VM, Name, Retired, KeyWord)
            }
        }
    """
    
    print "private_dictionaries = {"
    
    for owner in sorted(idictionary.keys()) :
        print "    {0} : {{".format(repr(owner))
        owners_dictionary = idictionary[owner]
        for key in sorted(owners_dictionary.keys()) :
            value = owners_dictionary[key]
            if isinstance(key, int) :
                print "        %#010x : %s,"%(key, value)
            else :
                print "        \"%8s\" : %s,"%(key, value)
        print "    },"
    print "}"
    
def print_public_dictionary(idictionary):
    """ Display Dictionary
    
        data_dictionary = {
            0xGroupElement : (VR, VM, Name, Retired, KeyWord),
            0xGroupElement : (VR, VM, Name, Retired, KeyWord),
            0xGroupElement : (VR, VM, Name, Retired, KeyWord)
        }
    """
    
    print "data_dictionary = {"
    
    for key in sorted(idictionary.keys()) :
        value = idictionary[key]
        if isinstance(key, int) :
            print "%#010x : %s,"%(key, value)
        else :
            print "\"%8s\" : %s,"%(key, value)
    
    print "}"
    
def main():
    """ Main function
        
        Read configuration file and
        Create private or public Dictionary
    """
    
    # Parse Argument
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--type', choices=['private', 'public'], help='Type of dictionary', required=True)
    parser.add_argument('configFile', help='Configuration file')
    args = parser.parse_args()
    
    # Available config type
    actype = { "GDCM"  : parse_gdcm,
               "DCMTK" : parse_dcmtk,
               "DICOM" : parse_dicom}
    
    # Read configuration
    conffile = read_configuration_file(args.configFile)
    
    data_dictionary = {}
    
    # Foreach item in configuration file
    for lcurrentSection in conffile.sections():
        # Get element type and url
        if conffile.has_option(lcurrentSection, "type") and conffile.has_option(lcurrentSection, "url"):
            urltype = conffile.get(lcurrentSection, "type")
            
            if urltype in actype.keys():
                dictionarytype = None
                if conffile.has_option(lcurrentSection, "dictionary_type"):
                    dictionarytype = conffile.get(lcurrentSection, "dictionary_type")
                    
                if dictionarytype != args.type:
                    continue
                
                privateName = None
                if conffile.has_option(lcurrentSection, "private_dico_name"):
                    privateName = conffile.get(lcurrentSection, "private_dico_name")
                    
                mainID = None
                if conffile.has_option(lcurrentSection, "mainid"):
                    mainID = conffile.get(lcurrentSection, "mainid")
                
                if not conffile.get(lcurrentSection, "url") == "":
                    data_dictionary.update(actype[urltype](conffile.get(lcurrentSection, "url"), args.type, mainID, privateName))
    
    # Available dictionary type
    adtype = { "private" : print_private_dictionary,
               "public"  : print_public_dictionary }
               
    # Display sorted dictionary
    adtype[args.type](data_dictionary)
    
    

if __name__ == "__main__":
    main()
    
