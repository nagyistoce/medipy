import docutils
import docutils.parsers.rst
import docutils.utils

from parameter_gui import ParameterGUI

def parseRSTAsParameters(docstring) :
    """ Create a list of ParameterGUI from a string 
    """
    
    parser = docutils.parsers.rst.Parser()
    settings = docutils.frontend.OptionParser(components=(docutils.parsers.rst.Parser,)).get_default_values()
    document = docutils.utils.new_document("/", settings)
    parser.parse(docstring, document)
    dom = document.asdom()
    root = dom.documentElement
    
    # Find the field element that has a field_name child with value gui
    fields = dom.getElementsByTagName("field")
    guiNode = None
    for field in fields :
        fieldName = field.getElementsByTagName("field_name")[0]
        if fieldName.firstChild.data == "gui" :
            guiNode = field.getElementsByTagName("definition_list")[0]
    
    # Find each element
    parameters = []
    for item in guiNode.childNodes :
        name = item.getElementsByTagName("term")[0].firstChild.data
        
        classifiers = item.getElementsByTagName("classifier")
        type = classifiers[0].firstChild.data
        
        initializer = None
        if len(classifiers) >= 2 :
            initializer = classifiers[1].firstChild.data
        
        definitions = item.getElementsByTagName("definition")[0].childNodes
        label = definitions[0].firstChild.data
        tooltip = len(definitions)>1 and definitions[1].firstChild.data or None
        parameter = ParameterGUI(name, type, initializer, label, tooltip)
        
        parameters.append(parameter)
    
    return parameters