import sys
import os.path

def find_resource(where) :
    """ Return the absolute file name of a resource inside the project.
        This function is necessary as bbfreeze cannot embed arbitrary files.
    """
    
    # Try first as a local application using MediPy
    command_name = sys.argv[0]
    abs_path = os.path.abspath(command_name)
    dirname = os.path.dirname(abs_path)
    elements = dirname.split(os.path.sep)
    
    resource_path = None
    while len(elements) > 0 :
        path = os.path.sep.join(elements)
        if os.path.exists(os.path.join(path, where)) :
            resource_path = os.path.join(path, where)
            break
        else :
            elements.pop()
    if resource_path is not None :
        return resource_path
    else :
        # Try in MediPy
        import medipy
        dirname = os.path.dirname(medipy.__file__)
        resource_path = os.path.join(dirname, where)
        if os.path.exists(resource_path) :
            return resource_path
        else :
            raise Exception("Cannot find resource %s"%where)