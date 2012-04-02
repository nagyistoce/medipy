import re

class ItkTypeWrapper(object):

    def __init__(self):
        self.templates = None
        self.includes = []
        
        self.ITKT = {
            "UC" : "unsigned char", "US" : "unsigned short", 
            "UI" : "unsigned int", "UL" : "unsigned long",
            "SC" : "signed char", "SS" : "signed short", 
            "SI" : "signed int", "SL" : "signed long",
            "F" : "float", "D" : "double", "LD" : "long double",
            "B" : "bool"}
        self.ITKM = dict([(key,key) for key in self.ITKT])
        self.ITKN = {}
        
        self._itk_Wrap_Class = None
        self._itk_Wrap_Prefix = None

    def wrap_type(self, class_, prefix):
            self.templates = []
            self._itk_Wrap_Prefix = prefix
            self._itk_Wrap_Class = class_
            
            if "itk::" in class_ :
                includeFileName = re.sub(r"itk.*::(.*)", r"itk\1", class_)
                self.includes.append("{0}.h".format(includeFileName))
    
    def add_template(self, name, types):
        self.templates.append((name, types))    
    
    def end_wrap_type(self):
        for wrapTpl, wrapType in self.templates :
            if "::" in self._itk_Wrap_Class :
                elements = self._itk_Wrap_Class.split("::")
                top_namespace, base_name = elements[0], elements[-1]
                swig_name = top_namespace+base_name
            else :
                swig_name = self._itk_Wrap_Class
            
            key = "{0}{1}".format(self._itk_Wrap_Prefix, wrapTpl)
            self.ITKT[key] = "{0}< {1} >".format(self._itk_Wrap_Class, wrapType)
            self.ITKM[key] = "{0}{1}".format(self._itk_Wrap_Prefix, wrapTpl)
            self.ITKN[key] = "{0}{1}".format(swig_name, wrapTpl)
