##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import re
import os
import subprocess

import medipy.base
from dataset import DataSet
from vr import *

def encapsulate(filename, mime_type=None):
    """ Create a DataSet encapsulating the content of a file. If not provided,
        the MIME type will be guessed from the content of the file.
        
        The Data Set only contains the following elements :
          * Document Title: contains the leaf of the filename.
          * Document Class Code Sequence: this sequence has only one element,
            identifying whether the length of the file is odd or even. This is
            necessary since the OB VR must be padded to be of even length, and
            the expose function must be able to determine the length of the
            original file.
          * MIME Type of Encapsulated Document.
          * Encapsulated Document.
    """
    
    dataset = DataSet()
    
    # Use the original filename as the document title
    dataset.document_title = ST(os.path.basename(filename))
    
    # Save the length parity information, since OB values must be padded by 0
    # if they are of odd length. Use the Document Class Code Sequence for this.
    length_parity = os.stat(filename).st_size%2
    length_parity = "LENGTH_PARITY:{0}".format({0:"EVEN", 1:"ODD"}[length_parity])
    dataset.document_class_code_sequence = SQ([DataSet()])
    dataset.document_class_code_sequence[0].code_value = SH(length_parity)
    dataset.document_class_code_sequence[0].coding_scheme_designator = SH("MEDIPY")
    dataset.document_class_code_sequence[0].code_meaning = LO(length_parity)
    
    # Try to find the MIME type if not provided
    if mime_type is None :
        process = subprocess.Popen(["file", "--mime-type", "--brief", filename],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if process.returncode != 0 :
            raise medipy.base.Exception("Cannot determine MIME type: {0!r}".format(stderr))
        mime_type = stdout
    dataset.mime_type_of_encapsulated_document = LO(mime_type)
    
    # Add the file data to the dataset
    fd = open(filename, "rb")
    dataset.encapsulated_document = OB(fd.read())
    fd.close()
    
    return dataset

def expose(dataset, filename):
    """ Save the Encapsulated Document in the dataset to a file. 
    """
    
    parity = None
    for item in dataset.get("document_class_code_sequence", SQ([])) :
        if item.get("coding_scheme_designator", "") == "MEDIPY" :
            match = re.match(r"LENGTH_PARITY:(EVEN|ODD)", item.code_value.value)
            if match :
                parity = match.group(1)
    
    fd = open(filename, "wb")
    if parity == "ODD" and len(dataset.encapsulated_document)%2==0 :
        fd.write(dataset.encapsulated_document.value[:-1])
    else :
        fd.write(dataset.encapsulated_document.value)
    fd.close()
