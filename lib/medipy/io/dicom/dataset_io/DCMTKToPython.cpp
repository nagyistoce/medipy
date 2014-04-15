/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include "DCMTKToPython.h"

#include <Python.h>

#include <map>
#include <stdexcept>
#include <string>

#include <dcmtk/config/osconfig.h>
#include <dcmtk/dcmdata/dctk.h>

DCMTKToPython
::DCMTKToPython()
: _medipy_io_dicom_vr(), _medipy_io_dicom_DataSet(NULL), 
  _medipy_io_dicom_Tag(NULL), _specific_character_set(""), _python_encoding("")
{
    PyObject* modules = PyImport_GetModuleDict(); // Borrowed reference
    PyObject* medipy = PyMapping_GetItemString(modules, const_cast<char*>("medipy")); // New reference
    PyObject* medipy_io = PyObject_GetAttrString(medipy, "io"); // New reference
    PyObject* medipy_io_dicom = PyObject_GetAttrString(medipy_io, "dicom"); // New reference
    
    this->_medipy_io_dicom_DataSet = PyObject_GetAttrString(medipy_io_dicom, "DataSet"); // New reference
    this->_medipy_io_dicom_Tag = PyObject_GetAttrString(medipy_io_dicom, "Tag"); // New reference
    
    // New references
#define ADD_VR_TO_MAP(vr) \
    this->_medipy_io_dicom_vr[EVR_##vr] = PyObject_GetAttrString(medipy_io_dicom, #vr)
    
    ADD_VR_TO_MAP(AE); ADD_VR_TO_MAP(AS); ADD_VR_TO_MAP(AT); ADD_VR_TO_MAP(CS);
    ADD_VR_TO_MAP(DA); ADD_VR_TO_MAP(DS); ADD_VR_TO_MAP(DT); ADD_VR_TO_MAP(FD);
    ADD_VR_TO_MAP(FL); ADD_VR_TO_MAP(IS); ADD_VR_TO_MAP(LO); ADD_VR_TO_MAP(LT);
    ADD_VR_TO_MAP(OB); ADD_VR_TO_MAP(OF); ADD_VR_TO_MAP(OW); ADD_VR_TO_MAP(PN);
    ADD_VR_TO_MAP(SH); ADD_VR_TO_MAP(SQ); ADD_VR_TO_MAP(SL); ADD_VR_TO_MAP(SS);
    ADD_VR_TO_MAP(ST); ADD_VR_TO_MAP(TM); ADD_VR_TO_MAP(UI); ADD_VR_TO_MAP(UL);
    ADD_VR_TO_MAP(UN); ADD_VR_TO_MAP(US); ADD_VR_TO_MAP(UT);

#undef ADD_VR_TO_MAP

    Py_DECREF(medipy_io_dicom);
    Py_DECREF(medipy_io);
    Py_DECREF(medipy);
}

DCMTKToPython
::DCMTKToPython(DCMTKToPython const & other)
: _medipy_io_dicom_vr(other._medipy_io_dicom_vr), 
  _medipy_io_dicom_DataSet(other._medipy_io_dicom_DataSet), 
  _medipy_io_dicom_Tag(other._medipy_io_dicom_Tag),
  _specific_character_set(""), _python_encoding("")
{
    for(std::map<DcmEVR, PyObject *>::iterator it=this->_medipy_io_dicom_vr.begin();
        it!=this->_medipy_io_dicom_vr.end(); ++it)
    {
        Py_XINCREF(it->second);
    }
    Py_XINCREF(this->_medipy_io_dicom_DataSet);
    Py_XINCREF(this->_medipy_io_dicom_Tag);
    
    this->set_specific_character_set(other.get_specific_character_set());
}

DCMTKToPython
::~DCMTKToPython()
{
    for(std::map<DcmEVR, PyObject *>::iterator it=this->_medipy_io_dicom_vr.begin();
        it!=this->_medipy_io_dicom_vr.end(); ++it)
    {
        Py_XDECREF(it->second);
    }
    Py_XDECREF(this->_medipy_io_dicom_DataSet);
    Py_XDECREF(this->_medipy_io_dicom_Tag);
}

DCMTKToPython &
DCMTKToPython
::operator=(DCMTKToPython const & other)
{
    if(this != &other)
    {
        for(std::map<DcmEVR, PyObject *>::iterator it=this->_medipy_io_dicom_vr.begin();
            it!=this->_medipy_io_dicom_vr.end(); ++it)
        {
            Py_XDECREF(it->second);
        }
        Py_XDECREF(this->_medipy_io_dicom_DataSet);
        Py_XDECREF(this->_medipy_io_dicom_Tag);
        
        this->_medipy_io_dicom_vr = other._medipy_io_dicom_vr;
        this->_medipy_io_dicom_DataSet = other._medipy_io_dicom_DataSet;
        this->_medipy_io_dicom_Tag = other._medipy_io_dicom_Tag;
        
        for(std::map<DcmEVR, PyObject *>::iterator it=this->_medipy_io_dicom_vr.begin();
            it!=this->_medipy_io_dicom_vr.end(); ++it)
        {
            Py_XINCREF(it->second);
        }
        Py_XINCREF(this->_medipy_io_dicom_DataSet);
        Py_XINCREF(this->_medipy_io_dicom_Tag);
        
        this->set_specific_character_set(other.get_specific_character_set());
    }
    return *this;
}

std::string const &
DCMTKToPython
::get_specific_character_set() const
{
    return this->_specific_character_set;
}

void
DCMTKToPython
::set_specific_character_set(std::string const & charset)
{
    if(charset != this->get_specific_character_set())
    {
        this->_specific_character_set = charset;
        
        if(charset == "") this->_python_encoding = "ascii";
        else if(charset == "ISO_IR 100") this->_python_encoding = "latin_1";
        else if(charset == "ISO_IR 101") this->_python_encoding = "iso8859_2";
        else if(charset == "ISO_IR 109") this->_python_encoding = "iso8859_3";
        else if(charset == "ISO_IR 110") this->_python_encoding = "iso8859_4";
        else if(charset == "ISO_IR 144") this->_python_encoding = "iso8859_5";
        else if(charset == "ISO_IR 127") this->_python_encoding = "iso8859_6";
        else if(charset == "ISO_IR 126") this->_python_encoding = "iso8859_7";
        else if(charset == "ISO_IR 138") this->_python_encoding = "iso8859_8";
        else if(charset == "ISO_IR 148") this->_python_encoding = "iso8859_9";
        else if(charset == "ISO_IR 13") this->_python_encoding = "iso2022_jp";
        // CP874 seems to be a superset of TIS-620/ISO-IR-166 (e.g. presence of 
        // the euro sign in the CP874 at an unassigned place of TIS-620), but we
        // should get away with it.
        else if(charset == "ISO_IR 166") this->_python_encoding = "cp874";
        // Single-byte character sets with code extensions (Table C.12-3)
        // ISO 2022 IR 6
        // ISO 2022 IR 100
        // ISO 2022 IR 101
        // ISO 2022 IR 109
        // ISO 2022 IR 110
        // ISO 2022 IR 144
        // ISO 2022 IR 127
        // ISO 2022 IR 126
        // ISO 2022 IR 138
        // ISO 2022 IR 148
        // ISO 2022 IR 113
        // ISO 2022 IR 166
        // Multi-byte character sets without code extensions (Table C.12-4)
        // ISO 2022 IR 87 // Kanji
        // ISO 2022 IR 159 // Supplementary Kanji
        // ISO 2022 IR 149 // Hangul, Hanja
        // Multi-byte character sets without code extensions (Table C.12-5)
        else if(charset == "ISO_IR 192") this->_python_encoding = "utf_8";
        else if(charset == "GB18030") this->_python_encoding = "gb18030";
    }
}

PyObject *
DCMTKToPython
::operator()(DcmObject * dataset)
{
    PyObject* python_dataset = PyObject_CallFunction(
        this->_medipy_io_dicom_DataSet, NULL);
    
    DcmObject * it = NULL;
    while(NULL != (it = dataset->nextInContainer(it)))
    {
        if(it->getTag() == DCM_SpecificCharacterSet)
        {
            // Specific Character Set: setup internal iconv converter
            DcmCodeString * specific_character_set = dynamic_cast<DcmCodeString*>(it);
            char* value;
            specific_character_set->getString(value);
            if(value != NULL)
            {
                this->set_specific_character_set(value);
            }
            else
            {
                this->set_specific_character_set("");
            }
        }
        
        if(it->getETag() == 0)
        {
            // Group length, do nothing
            continue;
        }
        else
        {
            this->_add_element(it, python_dataset);
        }
    }
    
    return python_dataset;
}

PyObject *
DCMTKToPython
::_to_python_tag(DcmTagKey const & dcmtk_tag) const
{
    PyObject* python_tag = PyObject_CallFunction(
        this->_medipy_io_dicom_Tag, const_cast<char*>("(II)"), 
        dcmtk_tag.getGroup(), dcmtk_tag.getElement());
    
    return python_tag;
}

/*******************************************************************************
 * Specializations of DCMTKToPython::_to_python for the different VRs.
 ******************************************************************************/

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_AE>(DcmObject * element) const
{
    return this->_to_python_text(dynamic_cast<DcmByteString*>(element), false);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_AS>(DcmObject * element) const
{
    return this->_to_python_text(dynamic_cast<DcmByteString*>(element), false);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_AT>(DcmObject * element) const
{
    DcmAttributeTag * attribute_tag = dynamic_cast<DcmAttributeTag*>(element);
    unsigned long const multiplicity = attribute_tag->getVM();
    
    PyObject * python_value = NULL;
    
    if(multiplicity > 1)
    {
        python_value = PyList_New(multiplicity);
        for(unsigned long position=0; position<multiplicity; ++position)
        {
            DcmTagKey dcmtk_tag;
            attribute_tag->getTagVal(dcmtk_tag, position);
            PyObject * python_tag = this->_to_python_tag(dcmtk_tag);
            PyList_SetItem(python_value, position, python_tag);
        }
    }
    else
    {
        DcmTagKey dcmtk_tag;
        attribute_tag->getTagVal(dcmtk_tag, 0);
        python_value = this->_to_python_tag(dcmtk_tag);
    }
    
    return python_value;
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_CS>(DcmObject * element) const
{
    return this->_to_python_text(dynamic_cast<DcmByteString*>(element), false);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_DA>(DcmObject * element) const
{
    return this->_to_python_text(dynamic_cast<DcmByteString*>(element), false);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_DT>(DcmObject * element) const
{
    return this->_to_python_text(dynamic_cast<DcmByteString*>(element), false);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_DS>(DcmObject * element) const
{
    PyObject * python_value = NULL;
    
    DcmDecimalString * ds = dynamic_cast<DcmDecimalString*>(element);
    unsigned long count = ds->getVM();

    if(count > 1)
    {
        python_value = PyList_New(count);
        
        for(unsigned long i=0; i<count; ++i)
        {
            Float64 value;
            ds->getFloat64(value, i);
            PyList_SetItem(python_value, i, PyFloat_FromDouble(value));
        }
    }
    else
    {
        Float64 value;
        ds->getFloat64(value, 0);
        python_value = PyFloat_FromDouble(value);
    }
    
    return python_value;
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_FD>(DcmObject * element) const
{
    return this->_to_python_number(dynamic_cast<DcmElement*>(element),
                                   &DcmElement::getFloat64);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_FL>(DcmObject * element) const
{
    return this->_to_python_number(dynamic_cast<DcmElement*>(element),
                                   &DcmElement::getFloat32);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_IS>(DcmObject * element) const
{
    PyObject * python_value = NULL;
    
    DcmIntegerString * is = dynamic_cast<DcmIntegerString*>(element);
    unsigned long count = is->getVM();

    if(count > 1)
    {
        python_value = PyList_New(count);
        
        for(unsigned long i=0; i<count; ++i)
        {
            Sint32 value;
            is->getSint32(value, i);
            PyList_SetItem(python_value, i, PyInt_FromLong(value));
        }
    }
    else
    {
        Sint32 value;
        is->getSint32(value, 0);
        python_value = PyInt_FromLong(value);
    }
    
    return python_value;
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_LO>(DcmObject * element) const
{
    return this->_to_python_text(dynamic_cast<DcmByteString*>(element), true);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_LT>(DcmObject * element) const
{
    return this->_to_python_text(dynamic_cast<DcmByteString*>(element), true);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_OB>(DcmObject * element) const
{
    return this->_to_python_binary(dynamic_cast<DcmElement*>(element));
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_OF>(DcmObject * element) const
{
    return this->_to_python_binary(dynamic_cast<DcmElement*>(element));
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_OW>(DcmObject * element) const
{
    return this->_to_python_binary(dynamic_cast<DcmElement*>(element));
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_PN>(DcmObject * element) const
{
    return this->_to_python_text(dynamic_cast<DcmByteString*>(element), true);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_SH>(DcmObject * element) const
{
    return this->_to_python_text(dynamic_cast<DcmByteString*>(element), true);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_SL>(DcmObject * element) const
{
    return this->_to_python_number(dynamic_cast<DcmElement*>(element),
                                   &DcmElement::getSint32);
}

// SQ is not processed here

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_SS>(DcmObject * element) const
{
    return this->_to_python_number(dynamic_cast<DcmElement*>(element),
                                   &DcmElement::getSint16);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_ST>(DcmObject * element) const
{
    return this->_to_python_text(dynamic_cast<DcmByteString*>(element), true);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_TM>(DcmObject * element) const
{
    return this->_to_python_text(dynamic_cast<DcmByteString*>(element), false);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_UI>(DcmObject * element) const
{
    return this->_to_python_text(dynamic_cast<DcmByteString*>(element), false);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_UL>(DcmObject * element) const
{
    return this->_to_python_number(dynamic_cast<DcmElement*>(element),
                                   &DcmElement::getUint32);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_UN>(DcmObject * element) const
{
    return this->_to_python_binary(dynamic_cast<DcmElement*>(element));
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_US>(DcmObject * element) const
{
    return this->_to_python_number(dynamic_cast<DcmElement*>(element),
                                   &DcmElement::getUint16);
}

template<>
PyObject *
DCMTKToPython
::_to_python<EVR_UT>(DcmObject * element) const
{
    return this->_to_python_text(dynamic_cast<DcmByteString*>(element), true);
}

/*******************************************************************************
 * End of specializations of DCMTKToPython::_to_python for the different VRs.
 ******************************************************************************/

#define TRY_DECODE(begin, size, encoding, python_value) \
    if(python_value==NULL) \
    { \
        /* Avoid propagating the error. */ \
        PyErr_Clear(); \
        python_value = PyUnicode_Decode(begin, size, encoding, "strict"); \
    }

PyObject *
DCMTKToPython
::_to_python_text(OFString const & value, DcmEVR const & vr, bool use_utf8) const
{
    // Explicit size is necessary due to the presence of NUL.
    // NUL is used to pad UI
    static OFString const whitespace(" \0", 2);
    
    PyObject * python_value = NULL;
    
    std::string::size_type first=0;
    std::string::size_type size=0;

    std::string::size_type const last = value.find_last_not_of(whitespace);
    
    if(last != std::string::npos)
    {
        if(vr != "LT" && vr != "PN" && vr != "ST" && vr != "UT")
        {
            // Leading spaces are significant for LT, PN, ST, and UT
            first = value.find_first_not_of(whitespace);
            if(first == std::string::npos)
            {
                first = 0;
            }
        }
        size = (last-first)+1;
    }
    
    char const * begin = value.c_str()+first;
    
    if(use_utf8)
    {
        TRY_DECODE(begin, size, this->_python_encoding.c_str(), python_value);
        TRY_DECODE(begin, size, "ascii", python_value);
        TRY_DECODE(begin, size, "latin_1", python_value);
        TRY_DECODE(begin, size, "iso8859_2", python_value);
        TRY_DECODE(begin, size, "iso8859_3", python_value);
        TRY_DECODE(begin, size, "iso8859_4", python_value);
        TRY_DECODE(begin, size, "iso8859_5", python_value);
        TRY_DECODE(begin, size, "iso8859_6", python_value);
        TRY_DECODE(begin, size, "iso8859_7", python_value);
        TRY_DECODE(begin, size, "iso8859_8", python_value);
        TRY_DECODE(begin, size, "iso8859_9", python_value);
        TRY_DECODE(begin, size, "iso2022_jp", python_value);
        TRY_DECODE(begin, size, "utf_8", python_value);
        TRY_DECODE(begin, size, "gb18030", python_value);
    }
    else
    {
        python_value = PyString_FromStringAndSize(begin, size);
    }
    
    return python_value;
}

#undef TRY_DECODE

PyObject *
DCMTKToPython
::_to_python_text(DcmByteString * element, bool use_utf8) const
{
    DcmEVR const vr = element->getVR();
    
    PyObject * python_value = NULL;
    
    unsigned long count = element->getVM();
    
    if(count>1)
    {
        python_value = PyList_New(count);
        for(unsigned long i=0; i<count; ++i)
        {
            OFString value;
            element->getOFString(value, i);
            PyList_SetItem(python_value, i, this->_to_python_text(value, vr, use_utf8));
        }
    }
    else
    {
        OFString value;
        element->getOFString(value, 0);
        python_value = this->_to_python_text(value, vr, use_utf8);
    }
    
    return python_value;
}

PyObject *
DCMTKToPython
::_to_python_binary(DcmElement * element) const
{
    PyObject * python_value = NULL;
    
    DcmOtherByteOtherWord* byte_string = dynamic_cast<DcmOtherByteOtherWord*>(element);
    if(byte_string != NULL)
    {
        char const * begin(NULL);
        if(element->getVR() != EVR_OW)
        {
            Uint8 * data;
            byte_string->getUint8Array(data);
            begin = reinterpret_cast<char const *>(data);
        }
        else
        {
            Uint16 * data;
            byte_string->getUint16Array(data);
            begin = reinterpret_cast<char const *>(data);
        }
        
        python_value = PyString_FromStringAndSize(begin, element->getLength());
    }
    else
    {
        throw std::runtime_error(
            std::string("Cannot handle conversion of ")+
                DcmVR(element->getVR()).getValidVRName());
    }
    
    return python_value;
}

void
DCMTKToPython
::_add_element(DcmObject * element, PyObject * python_dataset) const
{
    DcmEVR const dcmtk_vr(DcmVR(element->getVR()).getValidEVR());
    char const * const vr_name = DcmVR(dcmtk_vr).getValidVRName();

    PyObject * python_value = NULL;
    
    if(dcmtk_vr == EVR_SQ)
    {
        DcmSequenceOfItems * sequence = dynamic_cast<DcmSequenceOfItems*>(element);
        python_value = PyList_New(sequence->card());
        
        DcmObject * sequence_it = NULL;
        unsigned int index=0;
        while(NULL != (sequence_it = sequence->nextInContainer(sequence_it)))
        {
            DCMTKToPython converter(*this);
            PyObject * item_value = converter(sequence_it);
            PyList_SetItem(python_value, index, item_value);
            ++index;
        }
    }
    else if(element->getLength() == 0)
    {
        python_value = Py_None;
        Py_INCREF(python_value);
    }
    else
    {
        if(dcmtk_vr == EVR_AE) { python_value = this->_to_python<EVR_AE>(element); }
        if(dcmtk_vr == EVR_AS) { python_value = this->_to_python<EVR_AS>(element); }
        if(dcmtk_vr == EVR_AT) { python_value = this->_to_python<EVR_AT>(element); }
        if(dcmtk_vr == EVR_CS) { python_value = this->_to_python<EVR_CS>(element); }
        if(dcmtk_vr == EVR_DA) { python_value = this->_to_python<EVR_DA>(element); }
        if(dcmtk_vr == EVR_DS) { python_value = this->_to_python<EVR_DS>(element); }
        if(dcmtk_vr == EVR_DT) { python_value = this->_to_python<EVR_DT>(element); }
        if(dcmtk_vr == EVR_FD) { python_value = this->_to_python<EVR_FD>(element); }
        if(dcmtk_vr == EVR_FL) { python_value = this->_to_python<EVR_FL>(element); }
        if(dcmtk_vr == EVR_IS) { python_value = this->_to_python<EVR_IS>(element); }
        if(dcmtk_vr == EVR_LO) { python_value = this->_to_python<EVR_LO>(element); }
        if(dcmtk_vr == EVR_LT) { python_value = this->_to_python<EVR_LT>(element); }
        if(dcmtk_vr == EVR_OB) { python_value = this->_to_python<EVR_OB>(element); }
        if(dcmtk_vr == EVR_OF) { python_value = this->_to_python<EVR_OF>(element); }
        if(dcmtk_vr == EVR_OW) { python_value = this->_to_python<EVR_OW>(element); }
        if(dcmtk_vr == EVR_PN) { python_value = this->_to_python<EVR_PN>(element); }
        if(dcmtk_vr == EVR_SH) { python_value = this->_to_python<EVR_SH>(element); }
        // SQ is not processed here
        if(dcmtk_vr == EVR_SL) { python_value = this->_to_python<EVR_SL>(element); }
        if(dcmtk_vr == EVR_SS) { python_value = this->_to_python<EVR_SS>(element); }
        if(dcmtk_vr == EVR_ST) { python_value = this->_to_python<EVR_ST>(element); }
        if(dcmtk_vr == EVR_TM) { python_value = this->_to_python<EVR_TM>(element); }
        if(dcmtk_vr == EVR_UI) { python_value = this->_to_python<EVR_UI>(element); }
        if(dcmtk_vr == EVR_UL) { python_value = this->_to_python<EVR_UL>(element); }
        if(dcmtk_vr == EVR_UN) { python_value = this->_to_python<EVR_UN>(element); }
        if(dcmtk_vr == EVR_US) { python_value = this->_to_python<EVR_US>(element); }
        if(dcmtk_vr == EVR_UT) { python_value = this->_to_python<EVR_UT>(element); }
    }
    
    if(python_value == NULL)
    {
        throw std::runtime_error(std::string("Unhandled VR: ") + DcmVR(dcmtk_vr).getVRName());
    }

    // Build the tag
    PyObject* tag = PyObject_CallFunction(this->_medipy_io_dicom_Tag, 
        const_cast<char*>("(II)"), element->getGTag(), element->getETag());
    
    // Build the value
    std::map<DcmEVR, PyObject *>::const_iterator const vr_it = 
        this->_medipy_io_dicom_vr.find(dcmtk_vr);
    if(vr_it==this->_medipy_io_dicom_vr.end())
    {
        throw std::runtime_error(std::string("Unknown MediPy VR:") + vr_name);
    }
    PyObject* MediPyVR = vr_it->second; // Borrowed reference

    PyObject* typed_value = PyObject_CallFunction(MediPyVR, 
        const_cast<char*>("(O)"), python_value);
    Py_DECREF(python_value);

    // Update the dictionary
    PyDict_SetItem(python_dataset, tag, typed_value);
    Py_DECREF(tag);
    Py_DECREF(typed_value);
}
