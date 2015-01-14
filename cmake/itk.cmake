find_package(WrapITK REQUIRED)

macro(itk_end_wrap_module_swig_interface)

  # Loop over the extra swig input files and copy them to the Typedefs directory
  foreach(source ${WRAPPER_LIBRARY_SWIG_INPUTS})
    file(COPY "${source}"
         DESTINATION "${WRAPPER_MASTER_INDEX_OUTPUT_DIR}")
    get_filename_component(basename ${source} NAME)
    set(dest "${WRAPPER_MASTER_INDEX_OUTPUT_DIR}/${basename}")
    WRAP_ITK_INSTALL("/Configuration/Typedefs" "${dest}")
#    set(SWIG_INTERFACE_MODULE_CONTENT "${SWIG_INTERFACE_MODULE_CONTENT}%import ${basename}\n")
  endforeach()

  # the list of .i files generated for the module
  set(SWIG_INTERFACE_FILES )

  # prepare dependencies
  set(DEPS )
  foreach(dep ${WRAPPER_LIBRARY_DEPENDS})
    set(DEPS ${DEPS} ${${dep}IdxFiles} ${${dep}SwigFiles})
# MEDIPY FIX: do not ITK dependencies to the .mdx, since igenerator looks for 
# them in the local directory
#    set(SWIG_INTERFACE_MDX_CONTENT "${dep}.mdx\n${SWIG_INTERFACE_MDX_CONTENT}")
  endforeach()

  # add some libs required by this module
  set(swig_libs )
  foreach(swig_lib ${WRAPPER_SWIG_LIBRARY_FILES})
    get_filename_component(basename ${swig_lib} NAME)
    set(swig_libs ${swig_libs} --swig-include ${basename})
    file(COPY "${swig_lib}"
      DESTINATION "${WRAPPER_MASTER_INDEX_OUTPUT_DIR}")
    set(dest "${WRAPPER_MASTER_INDEX_OUTPUT_DIR}/${basename}")
    WRAP_ITK_INSTALL("/Configuration/Typedefs" "${dest}")
  endforeach()

  # accumulate the idx files already generated before this module to generate usable depenencies
  set(idx_files )
  foreach(module ${SWIG_INTERFACE_MODULES})
    # create the swig interface
    set(interface_file "${WRAPPER_MASTER_INDEX_OUTPUT_DIR}/${module}.i")
    set(xml_file "${WRAPPER_LIBRARY_OUTPUT_DIR}/${module}.xml")
    set(idx_file "${WRAPPER_MASTER_INDEX_OUTPUT_DIR}/${module}.idx")
    set(idx_files ${idx_files} ${idx_file})
#    set(includes_file "${WRAPPER_MASTER_INDEX_OUTPUT_DIR}/${WRAPPER_LIBRARY_NAME}.includes")
    set(module_interface_file "${WRAPPER_MASTER_INDEX_OUTPUT_DIR}/${WRAPPER_LIBRARY_NAME}.i")
    set(typedef_in_file "${WRAPPER_LIBRARY_OUTPUT_DIR}/${module}SwigInterface.h.in")
    set(typedef_file "${WRAPPER_MASTER_INDEX_OUTPUT_DIR}/${module}SwigInterface.h")
    # prepare the options
    set(opts )
    foreach(dep ${WRAPPER_LIBRARY_DEPENDS})
      set(opts ${opts} --mdx "${WRAP_ITK_TYPEDEFS_DIRECTORY}/${dep}.mdx")
#      set(opts ${opts} --include "${dep}.includes")
#      set(opts ${opts} --import "${dep}.i")
    endforeach()
    # import the interface files previously defined instead of importing all the files defined
    # in the module to avoid many warnings when running swig
#    foreach(i ${SWIG_INTERFACE_FILES})
#      get_filename_component(bi "${i}" NAME)
#      set(opts ${opts} --import "${bi}")
#    endforeach()

    if(ITK_WRAP_EXPLICIT)
      set(opts ${opts} --include "${WRAPPER_LIBRARY_NAME}Explicit.h")
    endif()

    add_custom_command(
      OUTPUT ${interface_file} ${typedef_file}
      COMMAND ${PYTHON_EXECUTABLE} ${IGENERATOR}
        ${opts}
        --swig-include itk.i
        ${swig_libs}
        --mdx ${mdx_file}
#        --include ${WRAPPER_LIBRARY_NAME}.includes
#        --import ${module_interface_file}
        --swig-include ${module}_ext.i
        -w1 -w3 -w51 -w52 -w53 -w54
        -A protected -A private
        -p ${PYGCCXML_DIR}
        -g ${GCCXML}
        --typedef-input ${typedef_in_file}
        --typedef-output ${typedef_file}
        --include ${module}SwigInterface.h
        ${xml_file}
        ${interface_file}
      DEPENDS ${DEPS} ${idx_files} ${IGENERATOR} # ${SWIG_INTERFACE_IDX_FILES} ${SWIG_INTERFACE_FILES} ${typedef_in_file} # ${includes_file}
    )
  #   add_custom_target(${module}Swig DEPENDS ${interface_file})
  #   add_dependencies(${module}Swig ${WRAPPER_LIBRARY_NAME}Idx)
  #   add_dependencies(${WRAPPER_LIBRARY_NAME}Swig ${module}Swig)

    set(SWIG_INTERFACE_FILES ${SWIG_INTERFACE_FILES} ${interface_file})

    WRAP_ITK_INSTALL("/Configuration/Typedefs" "${interface_file}")

  endforeach()


  # the mdx file
  set(mdx_file "${WRAPPER_MASTER_INDEX_OUTPUT_DIR}/${WRAPPER_LIBRARY_NAME}.mdx")
  set(CONFIG_INDEX_FILE_CONTENT "${SWIG_INTERFACE_MDX_CONTENT}")
  configure_file("${ITK_WRAP_SWIGINTERFACE_SOURCE_DIR}/Master.mdx.in" "${mdx_file}"
     @ONLY)
  WRAP_ITK_INSTALL("/Configuration/Typedefs" "${mdx_file}")

  set(module_interface_file "${WRAPPER_MASTER_INDEX_OUTPUT_DIR}/${WRAPPER_LIBRARY_NAME}.i")
  set(module_interface_ext_file "${WRAPPER_MASTER_INDEX_OUTPUT_DIR}/${WRAPPER_LIBRARY_NAME}_ext.i")
  set(deps_imports )
  set(deps_includes )
  foreach(dep ${WRAPPER_LIBRARY_DEPENDS})
    set(deps_imports "${deps_imports}%import ${dep}.i\n")
#    set(deps_includes "${deps_includes}#include \"${dep}.includes\"\n")
  endforeach()

#  set(SWIG_INTERFACE_INCLUDES "${deps_includes}#include \"${WRAPPER_LIBRARY_NAME}.includes\"")
  set(CONFIG_MODULE_INTERFACE_CONTENT ) #"${deps_imports}${SWIG_INTERFACE_MODULE_CONTENT}")
  configure_file("${ITK_WRAP_SWIGINTERFACE_SOURCE_DIR}/module.i.in" "${module_interface_file}"
    @ONLY)
  WRAP_ITK_INSTALL("/Configuration/Typedefs/" "${module_interface_file}")

#  set(WRAP_ITK_FILE_CONTENT )
#  configure_file("${WRAP_ITK_CONFIG_DIR}/empty.in" "${module_interface_ext_file}"
#    @ONLY)
#  install(FILES "${module_interface_ext_file}"
#    DESTINATION "${WRAP_ITK_INSTALL_PREFIX}/Configuration/Typedefs/"
#  )

  add_custom_target(${WRAPPER_LIBRARY_NAME}Idx DEPENDS ${SWIG_INTERFACE_IDX_FILES})
  add_dependencies(${WRAPPER_LIBRARY_NAME}Idx ${WRAPPER_LIBRARY_NAME}GccXML)
  set(${WRAPPER_LIBRARY_NAME}IdxFiles ${SWIG_INTERFACE_IDX_FILES} CACHE INTERNAL "Internal ${WRAPPER_LIBRARY_NAME}Idx file list.")

  add_custom_target(${WRAPPER_LIBRARY_NAME}Swig DEPENDS ${SWIG_INTERFACE_FILES})
  set(${WRAPPER_LIBRARY_NAME}SwigFiles ${SWIG_INTERFACE_FILES} CACHE INTERNAL "Internal ${WRAPPER_LIBRARY_NAME}Swig file list.")
  add_dependencies(${WRAPPER_LIBRARY_NAME}Swig ${WRAPPER_LIBRARY_NAME}Idx)
  if(NOT EXTERNAL_WRAP_ITK_PROJECT)
    # don't depends on the targets from wrapitk in external projects
    foreach(dep ${WRAPPER_LIBRARY_DEPENDS})
      add_dependencies(${WRAPPER_LIBRARY_NAME}Swig ${dep}Swig ${dep}Idx)
    endforeach()
  endif()

endmacro()

macro(itk_end_wrap_module_python)

  # Loop over the extra swig input files and add them to the generated files
  # lists. Guess that the generated cxx output will have the same name as
  # the .i input file.
  set(ITK_WRAP_PYTHON_PROCCESS_SWIG_INPUTS ON)
  foreach(source ${WRAPPER_LIBRARY_SWIG_INPUTS})
    get_filename_component(base_name ${source} NAME_WE)
    itk_wrap_submodule_python("${base_name}")
    itk_end_wrap_submodule_python("${base_name}")
  endforeach()
  set(ITK_WRAP_PYTHON_PROCCESS_SWIG_INPUTS OFF)

  # create the python config file
  # this file store all the name - type association and a dependencies list for the modules
  #
  # first build the dependency list
  set(ITK_WRAP_PYTHON_CONFIGURATION_DEPENDS "")

  foreach(dep ${WRAPPER_LIBRARY_DEPENDS})
    set(ITK_WRAP_PYTHON_CONFIGURATION_DEPENDS "'${dep}', ${ITK_WRAP_PYTHON_CONFIGURATION_DEPENDS}")
    set(ITK_WRAP_PYTHON_LIBRARY_IMPORTS "import ${dep}Python\n${ITK_WRAP_PYTHON_LIBRARY_IMPORTS}")
  endforeach()

  # ITKPyBase is always included, excepted ITKPyBase itself
  if(NOT "${WRAPPER_LIBRARY_NAME}" STREQUAL "ITKPyBase")
    set(ITK_WRAP_PYTHON_CONFIGURATION_DEPENDS "'ITKPyBase', ${ITK_WRAP_PYTHON_CONFIGURATION_DEPENDS}")
    set(ITK_WRAP_PYTHON_LIBRARY_IMPORTS "import ITKPyBasePython\n${ITK_WRAP_PYTHON_LIBRARY_IMPORTS}")
  endif()
  
  # MEDIPY FIX: Generate all wrappers in the same location so that sys.path 
  # configuration is easier when running from a build tree
  
  # and create the file, with the var ITK_WRAP_PYTHON_CONFIGURATION_TEMPLATES and
  # ITK_WRAP_PYTHON_CONFIGURATION_DEPENDS created earlier
  
  configure_file("${ITK_WRAP_PYTHON_SOURCE_DIR}/ModuleConfig.py.in"
    "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${WRAPPER_LIBRARY_NAME}Config.py"
    @ONLY)
  WRAP_ITK_BINDINGS_INSTALL("/Python/Configuration"
    "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${WRAPPER_LIBRARY_NAME}Config.py"
  )

  # create the advanced lib module python file
  # this file let the final user _not_ use the itk module, but rather
  # something like "import Base"
  set(CONFIG_LIBRARY_NAME "${WRAPPER_LIBRARY_NAME}")
  configure_file("${ITK_WRAP_PYTHON_SOURCE_DIR}/ModuleLoader.py.in"
    "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${WRAPPER_LIBRARY_NAME}.py"
    @ONLY)
  WRAP_ITK_BINDINGS_INSTALL("/Python" "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${WRAPPER_LIBRARY_NAME}.py")
  
  # END MEDIPY FIX

  # create the python customization stuff in the main module
  # it allow to group the pythons module in a single shared lib, by loading the int
  # functions of the module. I also import the objects from the submodules in the
  # main module.
  #
  # It uses ITK_WRAP_PYTHON_LIBRARY_DECLS, ITK_WRAP_PYTHON_LIBRARY_CALLS and
  # ITK_WRAP_PYTHON_LIBRARY_IMPORTS
  configure_file("${ITK_WRAP_PYTHON_SOURCE_DIR}/main_module_ext.i.in"
    "${WRAPPER_MASTER_INDEX_OUTPUT_DIR}/python/${WRAPPER_LIBRARY_NAME}_ext.i"
    @ONLY)
  WRAP_ITK_INSTALL("/Configuration/Typedefs/python"
    "${WRAPPER_MASTER_INDEX_OUTPUT_DIR}/python/${WRAPPER_LIBRARY_NAME}_ext.i"
  )

  # set some var reused later
  set(interface_file "${WRAPPER_MASTER_INDEX_OUTPUT_DIR}/${WRAPPER_LIBRARY_NAME}.i")
  set(lib ${WRAPPER_LIBRARY_NAME}Python)
  set(python_file "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${lib}.py")
  set(cpp_file "${CMAKE_CURRENT_BINARY_DIR}/${WRAPPER_LIBRARY_NAME}Python.cpp")

  # if this is for an external library, let the user add extra swig args
  if(EXTERNAL_WRAP_ITK_PROJECT)
    set(WRAP_ITK_SWIG_ARGS_PYTHON "" CACHE STRING "Extra user-defined swig arguments to be to the swig executable.")
  endif()

  # Run swig to produce the .so and the .py file
  itk_setup_swig_python("module" ${base_name} ${interface_file} ${python_file} ${cpp_file} "")

  # build all the c++ files from this module in a common lib
  add_library(${lib} MODULE ${cpp_file} ${ITK_WRAP_PYTHON_CXX_FILES} ${WRAPPER_LIBRARY_CXX_SOURCES})
  set_target_properties(${lib} PROPERTIES PREFIX "_")
  # gcc 4.4 complains a lot without this flag when building in release mode
  if(CMAKE_COMPILER_IS_GNUCC)
    set_target_properties(${lib} PROPERTIES COMPILE_FLAGS "-fno-strict-aliasing -Wno-shadow")
  endif()
  # work around linkage problem on some solaris systems
  if(CMAKE_SYSTEM MATCHES "SunOS-." AND CMAKE_COMPILER_IS_GNUCXX AND CMAKE_COMPILER_IS_GNUCC)
    target_link_libraries(${lib} stdc++)
  endif()
  # extension is not the same on windows
  if(WIN32)
    set_target_properties(${lib} PROPERTIES SUFFIX .pyd)
  endif()
  target_link_libraries(${lib} ${WRAPPER_LIBRARY_LINK_LIBRARIES} ${PYTHON_LIBRARY})
  add_dependencies(${lib} ${WRAPPER_LIBRARY_NAME}Swig)
  if(ITK_WRAP_DOC)
    add_dependencies(${lib} ${WRAPPER_LIBRARY_NAME}Doxygen)
  endif()
  if(ITK_WRAP_EXPLICIT AND NOT ${WRAPPER_LIBRARY_NAME} STREQUAL ITKPyBase)
    target_link_libraries(${lib} ${WRAPPER_LIBRARY_NAME}Explicit)
    add_dependencies(${lib} ${WRAPPER_LIBRARY_NAME}Explicit)
  endif()
  # 64 bit binaries are not installed in the same directories on solaris
  set(sun64 )
  if(CMAKE_SYSTEM MATCHES "SunOS.*" AND CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(sun64 "/64")
  endif()
  install(TARGETS "${lib}" DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/ITK-${ITK_VERSION_MAJOR}.${ITK_VERSION_MINOR}/Python${sun64}")
  if(NOT EXTERNAL_WRAP_ITK_PROJECT)
    # don't depends on the targets from wrapitk in external projects
    foreach(dep ${WRAPPER_LIBRARY_DEPENDS})
      add_dependencies(${lib} ${dep}Swig)
      if(ITK_WRAP_DOC)
        add_dependencies(${lib} ${dep}Doxygen)
      endif()
    endforeach()
  endif()

endmacro()

macro(find_swig_library_files)
    foreach(module ${ARGN})
        set(WRAPPER_SWIG_LIBRARY_FILES 
            ${WRAPPER_SWIG_LIBRARY_FILES} "${CMAKE_CURRENT_SOURCE_DIR}/itk${module}.swg")
    endforeach()
endmacro()

wrap_type("itk::Matrix" "M")
    unique(types "${WRAP_ITK_SCALAR};D")
    foreach(d ${ITK_WRAP_DIMS})
        foreach(type ${types})
            add_template("${ITKM_${type}}${d}${d}"  "${ITKT_${type}},${d},${d}")
        endforeach()
    endforeach()
end_wrap_type()
set(itk_Wrap_Matrix ${WRAPPER_TEMPLATES})
