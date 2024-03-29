##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# If uuid module is loaded too late, TLS might be uninitialized and calls to
# uuid functions will cause a segfault on some Linux versions. See
# http://www.sourceware.org/bugzilla/show_bug.cgi?id=12453
import uuid

# If expat is loaded too late, WrapITK causes crashes when an exception is 
# raised. See https://issues.itk.org/jira/browse/HISTITK-834
import xml.parsers.expat

import ConfigParser
import os
import sys
import warnings

# Hide some WrapITK warnings
warnings.filterwarnings(
    "ignore", ".*Unknown parameter.* in template", module="itkTemplate")

def configure_plugins_path() :
    plugins_path = []
    if "MEDIPY_PLUGINS_PATH" in os.environ :
        plugins_path = [
            x for x in os.environ["MEDIPY_PLUGINS_PATH"].split(os.pathsep)
            if os.path.isdir(x)
        ]
    elif sys.platform == "linux2" :
        """ The the global plugins path on Linux platform, based on the contents
            of /etc/medipy.cfg
        """
        
        plugins_path = []
        config_file = "/etc/medipy.cfg"
        if os.path.isfile(config_file) :
            config = ConfigParser.ConfigParser()
            config.read(config_file)
            if config.has_section("plugins") and config.has_option("plugins", "path") :
                plugins_path = [
                    x for x in config.get("plugins", "path").split(os.pathsep)
                    if os.path.isdir(x)
                ]
    for path in plugins_path :
        if path not in __path__ :
            __path__.append(path)

if not hasattr(sys, "frozen") :
    configure_plugins_path()
