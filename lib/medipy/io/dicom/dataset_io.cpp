/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include "dataset_io.h"

#include <fstream>
#include <string>

#include <gdcmBinEntry.h>
#include <gdcmDocEntry.h>
#include <gdcmDocument.h>
#include <gdcmUtil.h>
#include <Python.h>

#include "DataSetBridge.h"

PyObject* read(std::string const & filename)
{
    gdcm::Document document;
    document.SetFileName(filename);
    document.Load();

    if(!document.IsReadable())
    {
        PyErr_SetString(PyExc_Exception, ("Cannot parse "+filename).c_str());
        return NULL;
    }

    // Load /all/ binary entries
    gdcm::DocEntry* entry = document.GetFirstEntry();
    while(entry != NULL)
    {
        gdcm::BinEntry* bin_entry = dynamic_cast<gdcm::BinEntry*>(entry);
        if(bin_entry != NULL)
        {
            document.LoadEntryBinArea(bin_entry);
        }
        entry = document.GetNextEntry();
    }

    DataSetBridge data_set_bridge(&document);
    return data_set_bridge.to_python();
}

unsigned int ComputeGroup0002Length(gdcm::DocEntrySet & document)
{
    unsigned int groupLength = 0;
    bool found0002 = false;

    // for each zero-level Tag in the DCM Header
    gdcm::DocEntry *entry = document.GetFirstEntry();
    while(entry)
    {
        if(entry->GetGroup() == 0x0002)
        {
            found0002 = true;

            if(entry->GetElement() != 0x0000)
            {
                //if ( (vr == "OB")||(vr == "OW")||(vr == "UT")||(vr == "SQ"))
                // (no SQ, OW, UT in group 0x0002;)
                if(entry->GetVR() == "OB")
                {
                    // explicit VR AND (OB, OW, SQ, UT) : 4 more bytes
                    groupLength +=  4;
                }
                groupLength += 2 + 2 + 4 + entry->GetLength();
            }
        }
        else if(found0002)
        {
            break;
        }

        entry = document.GetNextEntry();
    }
    return groupLength;
}

void write(PyObject* dictionary, std::string const & filename)
{
    DataSetBridge data_set_bridge(dictionary);
    gdcm::Document document;
    data_set_bridge.to_gdcm(document);

    // Set the File Meta Information

    // File Meta Information Group Length
    uint32_t groupLength = ComputeGroup0002Length(document);
    document.InsertBinEntry(
        reinterpret_cast<uint8_t*>(&groupLength), sizeof(groupLength),
        0x0002, 0x0000, "UL");

    // File Meta Information Version
    uint8_t file_meta_information_version[2] = {0, 1};
    document.InsertBinEntry(file_meta_information_version, 2, 0x0002, 0x0001, "OB");

    // Media Storage SOP Instance UID
    std::string sop = gdcm::Util::CreateUniqueUID();
    document.InsertValEntry(sop, 0x0002, 0x0003, "UI");

    // Media Storage Implementation Class UID
    document.InsertValEntry(
        gdcm::Util::GetRootUID() + ".147.144.143.155." GDCM_VERSION,
        0x0002, 0x0012, "UI");

    // Implementation Version Name
    std::string implementation_version_name = "ITK/GDCM ";
    implementation_version_name += gdcm::Util::GetVersion();
    document.InsertValEntry(implementation_version_name, 0x0002, 0x0013, "SH");

    std::ofstream stream(filename.c_str());
    document.WriteContent(&stream, gdcm::ExplicitVR);
    stream.close();
}
