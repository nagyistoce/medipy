##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import struct

def parse_csa(csa):
    """ Return a dictionary of (tag, items) of the CSA.
    """
    
    format = ("<"  # Little-endian
              "4s" # SV10
              "4s" # \x04\x03\x02\x01
              "I"  # Number of items
              "I"  # Unknown
             )
    size = struct.calcsize(format)
    
    version, unknown_1, number_of_elements, unknown_2 = struct.unpack(format, csa[:size])
    
    start = size
    content = {}
    for _ in range(number_of_elements) :
        (name, items), size = parse_element(csa, start)
        content[name] = items
        start += size
    
    return content

def parse_element(csa, start):
    """ Return a pair (name, items), total_size
    """
    format = ("<"   # Little endian
              "64s" # Name
              "I"   # VM
              "2s"  # VR
              "2s"  # Unknown (end of VR ?)
              "I"   # Syngo datatype
              "I"   # Number of items
              "I"   # Unknown
            )
    size = struct.calcsize(format)
    
    name, vm, vr, unknown_1, syngo_datatype, number_of_items, unknown_2 = struct.unpack(
        format, csa[start:start+size])
    name = name.split("\x00")[0]
    
    total_size = size
    start += size
    items = []
    for i in range(number_of_items) :
        item, size = parse_item(csa, start)
        if i < vm : 
            if vr in ["DS", "FL", "FD"] :
                item = float(item[:-1])
            elif vr in ["IS", "SS", "US", "SL", "UL"] :
                item = int(item[:-1])
            items.append(item)
        start += size
        total_size += size
    
    return (name, items), total_size

def parse_item(csa, start):
    """ Return a pair content, size
    """
    format = ("<"  # Little endian
              "4I" # Length
             )
    header_size = struct.calcsize(format)
    
    length = struct.unpack(format, csa[start:start+header_size])
    
    format = ("<"    # Little endian
              "{0}s" # Content
              "{1}s" # Padding (?)
             ).format(length[1], (4-length[1]%4)%4)
    content_size = struct.calcsize(format)
    content, padding = struct.unpack(format, csa[start+header_size:start+header_size+content_size])
    
    return content, header_size+content_size