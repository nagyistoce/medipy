#include "PythonToDCMTK.h"

#include <Python.h>

#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>

#include <dcmtk/config/osconfig.h>
#include <dcmtk/dcmdata/dctk.h>

PythonToDCMTK
::PythonToDCMTK()
{
    // Nothing else.
}

PythonToDCMTK
::PythonToDCMTK(PythonToDCMTK const & )
{
    // Nothing else.
}

PythonToDCMTK
::~PythonToDCMTK()
{
    // Nothing to do.
}

PythonToDCMTK &
PythonToDCMTK
::operator=(PythonToDCMTK const & )
{
    // Nothing else
    
    return *this;
}

DcmDataset
PythonToDCMTK
::operator()(PyObject * python_dataset)
{
    DcmDataset dataset;
    
    PyObject *tag, *value;
    Py_ssize_t pos = 0;
    while(PyDict_Next(python_dataset, &pos, &tag, &value))
    {
        this->_add_element(tag, value, dataset);
    }
    
    return dataset;
}

/*******************************************************************************
 * Specializations of PythonToDCMTK::_to_dcmtk for the different VRs.
 ******************************************************************************/

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_AE>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_text(python_element, false, ' ', dataset, tag);
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_AS>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_text(python_element, false, ' ', dataset, tag);
}

// TODO : AT

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_CS>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_text(python_element, false, ' ', dataset, tag);
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_DA>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_text(python_element, false, ' ', dataset, tag);
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_DS>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_number_string(python_element, dataset, tag);
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_DT>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_text(python_element, false, ' ', dataset, tag);
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_FD>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_binary(python_element, PyFloat_AsDouble, dataset, tag,
                     &DcmDataset::putAndInsertFloat64);
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_FL>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_binary(python_element, PyFloat_AsDouble, dataset, tag,
                     &DcmDataset::putAndInsertFloat32);
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_IS>(PyObject * python_element, DcmDataset & dataset,
                   DcmTag const & tag) const
{
    this->_to_number_string(python_element, dataset, tag);
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_LO>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_text(python_element, true, ' ', dataset, tag);
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_LT>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_text(python_element, true, ' ', dataset, tag);
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_OB>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_raw(python_element, dataset, tag);
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_OF>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_raw(python_element, dataset, tag);
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_OW>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_raw(python_element, dataset, tag);
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_PN>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_text(python_element, true, ' ', dataset, tag);
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_SH>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_text(python_element, true, ' ', dataset, tag);
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_SL>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_binary(python_element, PyInt_AsLong, dataset, tag,
                     &DcmDataset::putAndInsertSint32);
}

// SQ is not processed here

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_SS>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_binary(python_element, PyInt_AsLong,  dataset, tag,
                     &DcmDataset::putAndInsertSint16);
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_ST>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_text(python_element, true, ' ', dataset, tag);
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_TM>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_text(python_element, false, ' ', dataset, tag);
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_UI>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_text(python_element, false, '\0', dataset, tag);
}


template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_UL>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_binary(python_element, PyInt_AsLong, dataset, tag,
                     &DcmDataset::putAndInsertUint32);
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_UN>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_raw(python_element, dataset, tag);
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_US>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_binary(python_element, PyInt_AsLong, dataset, tag,
                     &DcmDataset::putAndInsertUint16);
}

template<>
void
PythonToDCMTK
::_to_dcmtk<EVR_UT>(PyObject * python_element, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    this->_to_text(python_element, true, ' ', dataset, tag);
}

/*******************************************************************************
 * End of specializations of BSONToDataSet::_to_dcmtk for the different VRs.
 ******************************************************************************/

void
PythonToDCMTK
::_add_element(PyObject * python_tag, PyObject * python_value, DcmDataset & dataset)
{
    // Get the tag
    unsigned long const numeric_tag = PyLong_AsLong(python_tag);
    DcmTag const tag(numeric_tag>>16, numeric_tag&0xffff);
    
    if(tag == DCM_SpecificCharacterSet)
    {
        // We encode everything to UTF-8, so we overwrite the Specific Character Set
        dataset.putAndInsertOFStringArray(tag, "ISO_IR 192");
        return;
    }
    
    // Get the VR
    PyObject* value_type = PyObject_GetAttrString(python_value, "__class__"); // New reference
    PyObject* python_vr = PyObject_GetAttrString(value_type, "__name__"); // New reference
    DcmVR const vr(PyString_AsString(python_vr));
    DcmEVR const evr(vr.getValidEVR());
    Py_DECREF(python_vr);
    Py_DECREF(value_type);

    // Get the value
    PyObject* nested_value = PyObject_GetAttrString(python_value, "value"); // New reference

    if(nested_value != Py_None)
    {
        if(evr == EVR_SQ)
        {
            for(Py_ssize_t index=0; index<PyList_Size(nested_value); ++index)
            {
                PythonToDCMTK converter(*this);
                DcmDataset * dcmtk_item = new DcmDataset(
                    converter(PyList_GetItem(nested_value, index)));
                dataset.insertSequenceItem(tag, dcmtk_item);
            }
        }
        else
        {
            if(evr == EVR_AE) this->_to_dcmtk<EVR_AE>(nested_value, dataset, tag);
            else if(evr == EVR_AS) this->_to_dcmtk<EVR_AS>(nested_value, dataset, tag);
            // else if(evr == EVR_AT) this->_to_dcmtk<EVR_AT>(nested_value, dataset, tag);
            else if(evr == EVR_CS) this->_to_dcmtk<EVR_CS>(nested_value, dataset, tag);
            else if(evr == EVR_DA) this->_to_dcmtk<EVR_DA>(nested_value, dataset, tag);
            else if(evr == EVR_DS) this->_to_dcmtk<EVR_DS>(nested_value, dataset, tag);
            else if(evr == EVR_DT) this->_to_dcmtk<EVR_DT>(nested_value, dataset, tag);
            else if(evr == EVR_FD) this->_to_dcmtk<EVR_FD>(nested_value, dataset, tag);
            else if(evr == EVR_FL) this->_to_dcmtk<EVR_FL>(nested_value, dataset, tag);
            else if(evr == EVR_IS) this->_to_dcmtk<EVR_IS>(nested_value, dataset, tag);
            else if(evr == EVR_LO) this->_to_dcmtk<EVR_LO>(nested_value, dataset, tag);
            else if(evr == EVR_LT) this->_to_dcmtk<EVR_LT>(nested_value, dataset, tag);
            else if(evr == EVR_OB) this->_to_dcmtk<EVR_OB>(nested_value, dataset, tag);
            else if(evr == EVR_OF) this->_to_dcmtk<EVR_OF>(nested_value, dataset, tag);
            else if(evr == EVR_OW) this->_to_dcmtk<EVR_OW>(nested_value, dataset, tag);
            else if(evr == EVR_PN) this->_to_dcmtk<EVR_PN>(nested_value, dataset, tag);
            else if(evr == EVR_SH) this->_to_dcmtk<EVR_SH>(nested_value, dataset, tag);
            // SQ is not processed here
            else if(evr == EVR_SL) this->_to_dcmtk<EVR_SL>(nested_value, dataset, tag);
            else if(evr == EVR_SS) this->_to_dcmtk<EVR_SS>(nested_value, dataset, tag);
            else if(evr == EVR_ST) this->_to_dcmtk<EVR_ST>(nested_value, dataset, tag);
            else if(evr == EVR_TM) this->_to_dcmtk<EVR_TM>(nested_value, dataset, tag);
            else if(evr == EVR_UI) this->_to_dcmtk<EVR_UI>(nested_value, dataset, tag);
            else if(evr == EVR_UL) this->_to_dcmtk<EVR_UL>(nested_value, dataset, tag);
            else if(evr == EVR_UN) this->_to_dcmtk<EVR_UN>(nested_value, dataset, tag);
            else if(evr == EVR_US) this->_to_dcmtk<EVR_US>(nested_value, dataset, tag);
            else if(evr == EVR_UT) this->_to_dcmtk<EVR_UT>(nested_value, dataset, tag);
            else
            {
                throw std::runtime_error(std::string("Unhandled VR: ") + vr.getValidVRName());
            }
        }
    }
    else
    {
        dataset.insertEmptyElement(tag);
    }

    Py_DECREF(nested_value);
}

void
PythonToDCMTK
::_to_text(PyObject * python_value, bool use_utf8, char padding,
           DcmDataset & dataset, DcmTag const & tag) const
{
    std::ostringstream stream;
    stream.imbue(std::locale("C"));
    
    if(PyList_Check(python_value))
    {
        Py_ssize_t const size = PyList_Size(python_value);
        for(Py_ssize_t index=0; index<size; ++index)
        {
            PyObject * python_item = PyList_GetItem(python_value, index);
            if(use_utf8)
            {
                PyObject * decoded_string = PyUnicode_AsEncodedString(
                    python_item,  "utf-8", "strict");
                if(decoded_string == NULL)
                {
                    throw std::runtime_error("Could not encode Unicode object to UTF-8");
                }
                stream << PyString_AsString(decoded_string);
                Py_DECREF(decoded_string);
            }
            else
            {
                stream << PyString_AsString(python_item);
            }
            
            if(index != size-1)
            {
                stream << "\\";
            }
        }
    }
    else
    {
        if(use_utf8)
        {
            PyObject * decoded_string = PyUnicode_AsEncodedString(
                python_value, "utf-8", "strict");
            if(decoded_string == NULL)
            {
                throw std::runtime_error("Could not encode Unicode object to UTF-8");
            }
            stream << PyString_AsString(decoded_string);
            Py_DECREF(decoded_string);
        }
        else
        {
            stream << PyString_AsString(python_value);
        }
    }
    
    OFString value(stream.str().c_str());

    if(value.size()%2!=0)
    {
        value += padding;
    }
    
    dataset.putAndInsertOFStringArray(tag, value);
}

template<typename TInserter, typename TGetter>
void
PythonToDCMTK
::_to_binary(PyObject * python_value, TGetter getter,
             DcmDataset & dataset, DcmTag const & tag, TInserter inserter) const
{
    if(PyList_Check(python_value))
    {
        for(Py_ssize_t index=0; index<PyList_Size(python_value); ++index)
        {
            (dataset.*inserter)(tag, getter(PyList_GetItem(python_value, index)), index, OFTrue);
        }
    }
    else
    {
        (dataset.*inserter)(tag, getter(python_value), 0, OFTrue);
    }
}

void
PythonToDCMTK
::_to_raw(PyObject * python_value, DcmDataset & dataset,
          DcmTag const & tag) const
{
    unsigned long size = PyString_Size(python_value);
    Uint8 const * data = reinterpret_cast<Uint8 const *>(PyString_AsString(python_value));
    dataset.putAndInsertUint8Array(tag, data, size);
}

void
PythonToDCMTK
::_to_number_string(PyObject * python_value, DcmDataset & dataset,
                    DcmTag const & tag) const
{
    std::ostringstream stream;
    stream.imbue(std::locale("C"));
    
    if(PyList_Check(python_value))
    {
        Py_ssize_t const size = PyList_Size(python_value);
        for(Py_ssize_t index=0; index<size; ++index)
        {
            PyObject * python_item = PyList_GetItem(python_value, index);
            stream << PyString_AsString(PyObject_Str(python_item));
            
            if(index != size-1)
            {
                stream << "\\";
            }
        }
    }
    else
    {
        stream << PyString_AsString(PyObject_Str(python_value));
    }
    
    OFString value(stream.str().c_str());

    if(value.size()%2!=0)
    {
        value += ' ';
    }
    
    dataset.putAndInsertOFStringArray(tag, value);
}
