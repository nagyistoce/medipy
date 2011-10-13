/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _7dbfbe3c_7822_11e0_a61c_0010188057be
#define _7dbfbe3c_7822_11e0_a61c_0010188057be

#include <Python.h>

#include <string>

/** @brief Parse a DICOM file.
 *  @return A new reference to a dictionary containing the DICOM elements.
 */
PyObject* parse_file(std::string const & filename);

#endif // _7dbfbe3c_7822_11e0_a61c_0010188057be