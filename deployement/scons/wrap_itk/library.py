class Library(object):
    def __init__(self, name, dependencies=None, linked_libraries=None,
                 modules=None):
        self.name = name
        self.dependencies = dependencies or []
        self.linked_libraries = linked_libraries or []
        self.modules = modules or []