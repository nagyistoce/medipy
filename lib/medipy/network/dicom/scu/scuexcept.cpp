/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/
 
#include <string>
 
#include "scuexcept.h"

SCUException
::SCUException(std::string const & message)
:message(message)
{
    //Do nothing
}

const char*
SCUException
::what() const throw()
{
    return this->message.c_str();
}
