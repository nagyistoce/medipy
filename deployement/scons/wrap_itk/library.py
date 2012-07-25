class Library(object):
    def __init__(self, name, dependencies=None, linked_libraries=None,
                 modules=None):
        self.name = name
        self.dependencies = dependencies or []
        self.linked_libraries = linked_libraries or []
        self.modules = modules or []
    
    def _get_includes(self):
        return [
            "itkCommand.h",
            # itkStatisticsLabelObject.h does not exist in ITK 3.14
            #"itkStatisticsLabelObject.h"
            "itkOffset.h", "itkVector.h", "itkCovariantVector.h", "itkContinuousIndex.h", 
            "itkArray.h", "itkFixedArray.h", "itkRGBPixel.h", "itkRGBAPixel.h", "complex",
            "itkSymmetricSecondRankTensor.h", "itkImage.h", "itkVectorImage.h", 
            "itkVariableLengthVector.h", "itkPoint.h", "itkLevelSetNode.h", 
            "itkFlatStructuringElement.h", "itkSpatialObject.h", "itkHistogram.h", 
            "itkLabelMap.h"]
    
    includes = property(_get_includes)
    