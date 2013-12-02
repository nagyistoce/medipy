find_package(DCMTK REQUIRED)

add_definitions(-DHAVE_CONFIG_H)

find_path(OFLOG "dcmtk/oflog/oflog.h")

if(NOT ${OFLOG} STREQUAL "OFLOG-NOTFOUND")
    swig_python_module(scu 
        scu.i scuexcept.cpp scu.cpp DcmSCU.cpp MedipyDcmSCU.cpp
        echoscu.cpp findscu.cpp getscu.cpp
        ../../../io/dicom/dataset_io/DCMTKToPython.cpp
        ../../../io/dicom/dataset_io/PythonToDCMTK.cpp
    )
endif()
    
set_directory_properties(PROPERTIES
    INCLUDE_DIRECTORIES "${DCMTK_INCLUDE_DIRS};${PYTHON_INCLUDE_DIRS}")

if("${DCMTK_LIBRARIES}" STREQUAL "")
    set(DCMTK_LIBRARIES dcmdata dcmimage dcmimgle dcmjpeg dcmnet dcmpstat 
        dcmqrdb dcmsr dcmtls ijg12 ijg16 ijg8 ofstd)
endif()
find_library(oflog oflog)
if(NOT ${oflog} STREQUAL "oflog-NOTFOUND")
    set(DCMTK_LIBRARIES ${DCMTK_LIBRARIES} oflog)
endif()

set(DCMTK_LIBRARIES ${DCMTK_LIBRARIES} dcmdata)

if(NOT ${OFLOG} STREQUAL "OFLOG-NOTFOUND")
    swig_link_libraries(scu
        dcmdsig ${DCMTK_LIBRARIES} wrap pthread)
    install_swig_python_module(scu)
endif()
