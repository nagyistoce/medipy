##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import os       
import re       # regular expressions handling module
import numpy
import scipy

import medipy.base

#   RBNMR     Reads processed Bruker NMR-data.
#
#   SYNTAX    A = rbnmr(Path)    
                                
def rbnmr(path, report_progress = None):
    
    if report_progress is not None :
            report_progress(0.)
    
    # Get back the current directory
    working_dir = os.path.dirname(path)  
    
    # A is the data structure returned by rbnmr()
    A = {}
        
    # Look for the file entitled "proc2s" and put it into A
    if os.path.isfile(os.path.join(working_dir, "proc2s")):
        A["Proc2s"] = readnmrpar(os.path.join(working_dir, "proc2s"))
    else :
        raise medipy.base.Exception("RBNMR : Impossible to find : " + os.path.join(working_dir, "proc2s"))
            
    # Other configuration files like "acqu2s", for instance, are stored higher in the filesystem
    split_list = working_dir.split(os.path.sep)
    dir_path = os.path.sep.join(split_list[:-2])
    
    # Look for the files entitled "procs", "acqus" and "acqu2s" and put them into A  
    if os.path.isfile(os.path.join(dir_path, "acqu2s")):
        A["Acqu2s"] = readnmrpar(os.path.join(dir_path, "acqu2s")) 
    
    if os.path.isfile(os.path.join(dir_path, "acqus")):
        A["Acqus"] = readnmrpar(os.path.join(dir_path, "acqus"))
    else :
        raise medipy.base.Exception("RBNMR : Impossible to find : " + os.path.join(dir_path,"acqus"))               

    if os.path.isfile(os.path.join(working_dir, "procs")):
        A["Procs"] = readnmrpar(os.path.join(working_dir, "procs"))
    else :
        raise medipy.base.Exception("RBNMR : Impossible to find : " + os.path.join(working_dir, "procs"))
    
    # Correct data for NC_proc-parameter
    coeff = numpy.power(2.,-float(A["Procs"]["NC_proc"]))
    # Read imaginary data if the file entitled "1i" exists
    if os.path.isfile(os.path.join(working_dir, "1i")):
        A["IData"] = numpy.fromfile(os.path.join(working_dir, "1i"), numpy.int32) 
        A["IData"] = A["IData"]/coeff
    
    # Set up the X axis
    A["XAxis"] = scipy.linspace(float(A["Procs"]["OFFSET"]),float(A["Procs"]["OFFSET"])-float(A["Procs"]["SW_p"])/float(A["Procs"]["SF"]),int(A["Procs"]["SI"]))

    # The Y axis is set up only if the file entitled "proc2s" exists
    if "Proc2s" in A:
            A["YAxis"] = scipy.linspace(float(A["Proc2s"]["OFFSET"]),float(A["Proc2s"]["OFFSET"])-float(A["Proc2s"]["SW_p"])/float(A["Proc2s"]["SF"]),int(A["Proc2s"]["SI"]))
    
    if report_progress is not None :
            report_progress(0.25)
            
    # Reorder submatrixes (se XWinNMR-manual, chapter 17.5 (95.3))
    SI1 = int(A["Procs"]["SI"])
    SI2 = int(A["Proc2s"]["SI"])
    
    XDIM1 = int(A["Procs"]["XDIM"])
    XDIM2 = int(A["Proc2s"]["XDIM"])

    NoSM = SI1*SI2/(XDIM1*XDIM2)    # Total number of submatrices
    NoSM1 = SI1/XDIM1               # Number of SM along first axis
    NoSM2 = SI2/XDIM2               # Number of SM along second axis
    
    # Get back data from the main file
    data = numpy.fromfile(path, numpy.int32)
    
    # Show data as an array of submatrices
    data = data.reshape(NoSM, XDIM2, XDIM1)
    
    if report_progress is not None :
            report_progress(0.5)
    
    # Re-order the submatrices in the final array
    A["Data"] = numpy.ndarray((SI2, SI1), dtype=data.dtype)
    for j in range(NoSM2) :
        for i in range(NoSM1) :
            submatrix = data[j*NoSM1+i]
            A["Data"][j*XDIM2:(j+1)*XDIM2,i*XDIM1:(i+1)*XDIM1] = submatrix
    # Scale the data
    A["Data"] /= coeff 
    
    if report_progress is not None :
            report_progress(0.75)
    
    # Read the "clevels" file if it exists    
    if os.path.isfile(os.path.join(working_dir, "clevels")):
        L = readnmrpar(os.path.join(working_dir, "clevels"))
        A["Levels"] = []
        if int(L["LEVSIGN"]) == 1:
            # Keep only positive data for A["Levels"]
            for l in L["LEVELS"].split()[1:]:
                if float(l)>0:
                    A["Levels"].append(float(l))
            A["Levels"] = numpy.asarray(A["Levels"])
            A["Levels"] = numpy.reshape(A["Levels"],[len(A["Levels"]), 1])
        
        elif int(L["LEVSIGN"]) == 2:
            # Keep only negative data for A["Levels"]
            for l in L["LEVELS"].split()[1:]:
                if float(l)<0:
                    A["Levels"].append(float(l))
            A["Levels"] = numpy.asarray(A["Levels"])
            A["Levels"] = numpy.reshape(A["Levels"],[len(A["Levels"]), 1])
            
        elif int(L["LEVSIGN"]) == 3:
            # Keep positive and negative data within a fixed range
            for l in L["LEVELS"].split()[1:int(L["MAXLEV"])*2+1]:
                A["Levels"].append(float(l))
            A["Levels"] = numpy.asarray(A["Levels"])
            A["Levels"] = numpy.reshape(A["Levels"],[len(A["Levels"]), 1])
        
    if report_progress is not None :
            report_progress(1.)
        
    # The function returns a data structure
    return A
            
# This method converts a parameter file into a data structure
def readnmrpar(filename):
    
    # P is the data structure returned by readnmrpar(filename)
    P = {}

    # List of all the regular expressions that can be found in files 
    regular_expressions = []
    regular_expressions.append(("^##\\$*(.+)=\\ \\(\\d\\.\\.\\d+\\)(.+)","ParVecVal"))
    regular_expressions.append(("^##\\$*(.+)=\\ \\(\\d\\.\\.\\d+\\)$","ParVec"))
    regular_expressions.append(("^##\\$*(.+)=\\ (.+)","ParVal"))
    regular_expressions.append(("^([^\\$#].*)","Val"))
    regular_expressions.append(("^\\$\\$(.+)","Stamp"))
    regular_expressions.append(("^##\\$*(.+)=","Empty"))

    # Read the input file
    f = open(filename,'r')
    lines = f.readlines()
    f.close()

    TypeOfRow = []
    row_value = []

    # Every single line is compared to all the regular expressions until a match is found
    for i in range(len(lines)):
        stop = 0
        for (regexpr, type_regexpr) in regular_expressions:
             pattern = regexpr
             found = re.match(pattern, lines[i])
             if found is not None and stop == 0:
                stop = 1
                a = re.match(pattern, lines[i])
                TypeOfRow.append((type_regexpr,a.group(1)))
                row_value.append(a.group(0).split()[1:])


    # Set up the structure
    for i in range(len(TypeOfRow)):
        if TypeOfRow[i][0] == "ParVal":
            LastParameterName = TypeOfRow[i][1]
            P[LastParameterName] = " ".join(row_value[i])
        elif TypeOfRow[i][0] == "ParVec":
            LastParameterName = TypeOfRow[i][1]
            P[LastParameterName]=""
        elif TypeOfRow[i][0] == "ParVecVal":
            LastParameterName = TypeOfRow[i][1]
            P[LastParameterName]=" ".join(row_value[i])
        elif TypeOfRow[i][0] == "Stamp":
            if "Stamp" not in P:
                P["Stamp"]=TypeOfRow[i][1]
            else:
                P["Stamp"]=P["Stamp"] + " " + TypeOfRow[i][1]
        elif TypeOfRow[i][0] == "Val":
            if P[LastParameterName] == "":
                P[LastParameterName] = TypeOfRow[i][1]
            else:
                P[LastParameterName] = str(P[LastParameterName]) + " " + TypeOfRow[i][1]
        elif TypeOfRow[i][0] == "Empty":
            pass
    
    return P

def point_to_ppm(point, procs, proc2s):
    """ Return the PPM value as (F1, F2) from an (y,x) point
    """
    
    # It seems that F1 is related to the Y axis, while F2 is related to the X axis
    
    begin = (float(proc2s["OFFSET"]), float(procs["OFFSET"]))
    # End is begin-sw_p/sf, so step is (end-begin)/si, which simplifies to
    # (-sw_p/sf+1)/si
    step = [(-float(p["SW_p"])/float(p["SF"]))/float(p["SI"]) 
            for p in [proc2s, procs] ]
    
    return [begin[i]+step[i]*point[i] for i in (0,1)]

def ppm_to_point(ppm, procs, proc2s):
    """ Return an (y,x) point from the PPM as (F1, F2)
    """
    
    # It seems that F1 is related to the Y axis, while F2 is related to the X axis
    
    begin = (float(proc2s["OFFSET"]), float(procs["OFFSET"]))
    # End is begin-sw_p/sf, so step is (end-begin)/si, which simplifies to
    # (-sw_p/sf+1)/si
    step = [(-float(p["SW_p"])/float(p["SF"]))/float(p["SI"]) 
            for p in [proc2s, procs] ]
    
    return [(ppm[i]-begin[i])/step[i] for i in (0,1)]

if __name__ == "__main__" : 
    
     import medipy.components.io as io
    
     a = io.load("/base_image/nmr/colon-HSQC/51mailio-C2-HSQC/4/pdata/1/2rr")
