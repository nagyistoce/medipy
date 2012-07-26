class Module(object):
    def __init__(self, name, classes=None):
        self.name = name
        self.classes = classes or []