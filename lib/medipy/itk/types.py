##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import itk

# Convert ImageIOBase::IOComponentType to its corresponding ITK type
io_component_type_to_type = {
    itk.ImageIOBase.UCHAR : itk.UC,
    itk.ImageIOBase.CHAR : itk.SC,
    itk.ImageIOBase.USHORT : itk.US,
    itk.ImageIOBase.SHORT : itk.SS,
    itk.ImageIOBase.UINT : itk.UI,
    itk.ImageIOBase.INT : itk.SI,
    itk.ImageIOBase.ULONG : itk.UL,
    itk.ImageIOBase.LONG : itk.SL,
    itk.ImageIOBase.FLOAT : itk.F,
    itk.ImageIOBase.DOUBLE : itk.D,
}

# Map an ITK type to its immediately smaller ITK type
smaller_type = {
    itk.US : itk.UC,
    itk.UL : itk.US,
    itk.UI : itk.US,
    itk.SS : itk.SC,
    itk.SL : itk.SS,
    itk.SI : itk.SS,
    itk.D : itk.F
}

# Map an ITK type to its immediately larger ITK type
larger_type = {
    itk.UC : itk.US,
    itk.US : itk.UL,
    itk.UI : itk.UL,
    itk.SC : itk.SS,
    itk.SS : itk.SL,
    itk.SI : itk.SL,
    itk.F : itk.D
}
