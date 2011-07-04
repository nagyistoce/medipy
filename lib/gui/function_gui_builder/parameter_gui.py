class ParameterGUI :
    """Simple class to store the GUI informations of a function parameter"""
    name = None
    type = None
    initializer = None
    label = None
    tooltip = None
    
    def __init__(self, name, type, initializer, label, tooltip) :
        self.name = name
        self.type = type
        self.initializer = initializer
        self.label = label
        self.tooltip = tooltip

