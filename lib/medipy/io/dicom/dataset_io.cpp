/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include "dataset_io.h"

#include <fstream>
#include <stdexcept>
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

void write(PyObject* dataset, std::string const & filename)
{
    gdcm::Document document;
    {
        DataSetBridge data_set_bridge(PyObject_GetAttrString(dataset, "header"));
        data_set_bridge.to_gdcm(document);
    }
    {
        DataSetBridge data_set_bridge(dataset);
        data_set_bridge.to_gdcm(document);
    }

    // Set the File Meta Information

    // File Meta Information Version
    uint8_t file_meta_information_version[2] = {0, 1};
    document.InsertBinEntry(file_meta_information_version, 2, 0x0002, 0x0001, "OB");

    // Media Storage SOP Instance UID
    std::string sop = gdcm::Util::CreateUniqueUID();
    document.InsertValEntry(sop, 0x0002, 0x0003, "UI");

    // Transfer Syntax UID
    if(document.GetDocEntry(0x0002, 0x0010) == NULL)
    {
        std::string transferSyntax = "1.2.840.10008.1.2";
        std::vector<uint8_t> paddedTransferSyntax(
            (transferSyntax.size()%2==0)?transferSyntax.size():1+transferSyntax.size(),
            0);
        std::copy(transferSyntax.begin(), transferSyntax.end(), paddedTransferSyntax.begin());

        document.InsertBinEntry(
            &paddedTransferSyntax[0], paddedTransferSyntax.size(),
            0x0002, 0x0010, "UI");
    }

    // Media Storage Implementation Class UID
    document.InsertValEntry(
        gdcm::Util::GetRootUID() + ".147.144.143.155." GDCM_VERSION,
        0x0002, 0x0012, "UI");

    // Implementation Version Name
    std::string implementation_version_name = "ITK/GDCM ";
    implementation_version_name += gdcm::Util::GetVersion();
    document.InsertValEntry(implementation_version_name, 0x0002, 0x0013, "SH");

    // File Meta Information Group Length
    uint32_t groupLength = ComputeGroup0002Length(document);
    document.InsertBinEntry(
        reinterpret_cast<uint8_t*>(&groupLength), sizeof(groupLength),
        0x0002, 0x0000, "UL");

    std::ofstream stream(filename.c_str());

    std::string transferSyntax((char*)(document.GetBinEntry(0x0002, 0x0010)->GetBinArea()),
                               document.GetBinEntry(0x0002, 0x0010)->GetLength());
    // Remove the padding
    static std::string const whitespace(" \0", 2);
    std::string::size_type const last = transferSyntax.find_last_not_of(whitespace);
    std::string::size_type first = 0;
    if(last != std::string::npos)
    {
        first = transferSyntax.find_first_not_of(whitespace);

        if(first == std::string::npos)
        {
            first = 0;
        }
    }
    transferSyntax = transferSyntax.substr(first, last-first+1);

    if(transferSyntax == "1.2.840.10008.1.2.1")
    {
        document.WriteContent(&stream, gdcm::ExplicitVR);
    }
    else if(transferSyntax == "1.2.840.10008.1.2")
    {
        document.WriteContent(&stream, gdcm::ImplicitVR);
    }
    else
    {
        std::string message = "Cannot write with a transfer syntax of '";
        message += transferSyntax+"'";
        throw std::runtime_error(message);
    }
    stream.close();
}
