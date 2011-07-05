import xml.dom.minidom

def parse_docstring(docstring):
    """ Parse the given docstring to a list of dictionaries, describing GUI
        parameters.
    """
    
    gui_description_begin = docstring.rfind("<gui>")
    gui_description_end = docstring.rfind("</gui>")
    gui_description = docstring[gui_description_begin:gui_description_end+len("</gui>")]
    
    parameters = []
    
    if not gui_description :
        return parameters
    
    document = xml.dom.minidom.parseString(gui_description)
    for node in document.childNodes[0].childNodes :
        if node.nodeType != node.ELEMENT_NODE or node.tagName != "item":
            continue
            
        parameters.append(dict(node.attributes.items()))

    return parameters