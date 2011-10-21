import copy

from itk_types import ItkTypeWrapper

wrapitk_root = "/usr/lib/InsightToolkit/WrapITK/"

WRAP_ITK_DIMS = [2,3]
    
WRAP = {}
WRAP["unsigned char"] = True
WRAP["unsigned long"] = False
WRAP["unsigned short"] = False
WRAP["signed char"] = False
WRAP["signed long"] = False
WRAP["signed short"] = True
WRAP["float"] = True
WRAP["double"] = True

WRAP["complex_float"] = True
WRAP["complex_double"] = False

WRAP["vector_double"] = False
WRAP["vector_float"] = True
WRAP["covariant_vector_double"] = False
WRAP["covariant_vector_float"] = True

WRAP["rgb_unsigned_char"] = True
WRAP["rgb_unsigned_short"] = False
WRAP["rgba_unsigned_char"] = True
WRAP["rgba_unsigned_short"] = False

unsigned_integers = [
    t for t in ["unsigned char", "unsigned long", "unsigned short"] if WRAP[t]]
signed_integers = [
    t for t in ["signed char", "signed long", "signed short"] if WRAP[t]]
reals = [t for t in ["float", "double"] if WRAP[t]]
rgb = []#t for t in ["rgb_unsigned_char", "rgb_unsigned_short"] if WRAP[t]]
rgba = []#t for t in ["rgba_unsigned_char", "rgba_unsigned_short"] if WRAP[t]]
vector_real = []#t for t in ["vector_double", "vector_float"] if WRAP[t]]
covariant_vector_real = [
    ]#t for t in ["covariant_vector_double", "covariant_vector_float"] if WRAP[t]]

integers = unsigned_integers+signed_integers
scalars = integers+reals
vectors = vector_real+covariant_vector_real
colors = rgb+rgba
all_types = integers+scalars+vectors+colors
dimensions = [1,2,3,4]

mapping = {
    "unsigned char" : "UC", "unsigned long" : "UL", "unsigned short" : "US",
    "signed char" : "SC", "signed long" : "SL", "signed short" : "SS",
    "float" : "F", "double" : "D",
    "rgb_unsigned_char" : "RGBUC", "rgb_unsigned_short" : "RGBUS",
    "rgba_unsigned_char" : "RGBAUC", "rgba_unsigned_short" : "RGBAUS",
    "vector_double" : "VD", "vector_float" : "VF",
    "covariant_vector_double" : "CVD", "covariant_vector_float" : "CVF",
    "complex_double" : "CD", "complex_float" : "CF",
}
reverse_mapping = dict([(v,k) for k,v in mapping.items()])

WRAP_ITK_USIGN_INT = [
    mapping[t] for t in ["unsigned char", "unsigned long", "unsigned short"] if WRAP[t]
]
WRAP_ITK_SIGN_INT = [
    mapping[t] for t in ["signed char", "signed long", "signed short"] if WRAP[t]
]
WRAP_ITK_REAL = [
    mapping[t] for t in ["float", "double"] if WRAP[t]
]
WRAP_ITK_RGB = [
    mapping[t] for t in ["rgb_unsigned_char", "rgb_unsigned_short"] if WRAP[t]
]
WRAP_ITK_RGBA = [
    mapping[t] for t in ["rgba_unsigned_char", "rgba_unsigned_short"] if WRAP[t]
]
WRAP_ITK_VECTOR_REAL = [
    mapping[t] for t in ["vector_double", "vector_float"] if WRAP[t]
]
WRAP_ITK_COV_VECTOR_REAL = [
    mapping[t] for t in ["covariant_vector_double", "covariant_vector_float"] if WRAP[t]
]
WRAP_ITK_COMPLEX_REAL = [
    mapping[t] for t in ["complex_double", "complex_float"] if WRAP[t]
]

WRAP_ITK_INT = WRAP_ITK_SIGN_INT+WRAP_ITK_USIGN_INT
WRAP_ITK_SCALAR = WRAP_ITK_INT+WRAP_ITK_REAL
WRAP_ITK_VECTOR = WRAP_ITK_VECTOR_REAL+WRAP_ITK_COV_VECTOR_REAL
WRAP_ITK_COLOR = WRAP_ITK_RGB+WRAP_ITK_RGBA
WRAP_ITK_ALL_TYPES = WRAP_ITK_RGB+WRAP_ITK_RGBA+WRAP_ITK_VECTOR+WRAP_ITK_SCALAR+WRAP_ITK_COMPLEX_REAL

SMALLER_THAN_D = set(WRAP_ITK_SCALAR).intersection(["F","UL","US","UC","SL","SS","SC"])
SMALLER_THAN_F = set(WRAP_ITK_SCALAR).intersection(["UL","US","UC","SL","SS","SC"])
SMALLER_THAN_UL = set(WRAP_ITK_INT).intersection(["US","UC","SL","SS","SC"])
SMALLER_THAN_US = set(WRAP_ITK_INT).intersection(["UC","SC"])
SMALLER_THAN_SL = set(WRAP_ITK_INT).intersection(["US","UC","SS","SC"])
SMALLER_THAN_SS = set(WRAP_ITK_INT).intersection(["UC","SC"])

type_wrapper = None

if type_wrapper is None :

    type_wrapper = ItkTypeWrapper()
    
    type_wrapper.wrap_type("itk::Offset", "O")
    for d in set(WRAP_ITK_DIMS+[1,2]) :
        type_wrapper.add_template(d, d)
    type_wrapper.end_wrap_type()
    
    type_wrapper.wrap_type("itk::Vector", "V")
    for d in set(WRAP_ITK_DIMS+[1,6]) :
        for t in set(WRAP_ITK_SCALAR+["UC","F","D"]) :
            type_wrapper.add_template("{0}{1}".format(type_wrapper.ITKM[t], d),
                                      "{0},{1}".format(type_wrapper.ITKT[t], d))
    type_wrapper.end_wrap_type()
    
    type_wrapper.wrap_type("itk::CovariantVector", "CV")
    for d in set(WRAP_ITK_DIMS) :
        type_wrapper.add_template("{0}{1}".format(type_wrapper.ITKM["F"], d),
                                  "{0},{1}".format(type_wrapper.ITKT["F"], d))
        type_wrapper.add_template("{0}{1}".format(type_wrapper.ITKM["D"], d),
                                  "{0},{1}".format(type_wrapper.ITKT["D"], d))
    type_wrapper.end_wrap_type()
    
    type_wrapper.wrap_type("itk::ContinuousIndex", "CI")
    for d in set(WRAP_ITK_DIMS) :
        type_wrapper.add_template("{0}{1}".format(type_wrapper.ITKM["F"], d),
                                  "{0},{1}".format(type_wrapper.ITKT["F"], d))
        type_wrapper.add_template("{0}{1}".format(type_wrapper.ITKM["D"], d),
                                  "{0},{1}".format(type_wrapper.ITKT["D"], d))
    type_wrapper.end_wrap_type()
    
    type_wrapper.wrap_type("itk::Array", "A")
    for t in ["D", "F", "UL", "SL"] :
        type_wrapper.add_template(type_wrapper.ITKM[t], type_wrapper.ITKT[t])
    type_wrapper.end_wrap_type()
    
    type_wrapper.wrap_type("itk::FixedArray", "FA")
    dims = copy.copy(WRAP_ITK_DIMS)
    for d in WRAP_ITK_DIMS :
        dims.append(d*2)
        dims.append(d*(d+1)/2)
    for d in set(dims+range(1,7)) :
        for t in ["D", "F", "UL", "US", "UC", "UI", "SL", "SS", "SC", "B"] :
            type_wrapper.add_template("{0}{1}".format(type_wrapper.ITKM[t], d),
                                      "{0},{1}".format(type_wrapper.ITKT[t], d))
    type_wrapper.end_wrap_type()
    
    type_wrapper.wrap_type("itk::RGBPixel", "RGB")
    type_wrapper.add_template(type_wrapper.ITKM["UC"], type_wrapper.ITKT["UC"])
    if WRAP["rgb_unsigned_short"] :
        type_wrapper.add_template(type_wrapper.ITKM["US"], type_wrapper.ITKT["US"])
    type_wrapper.end_wrap_type()
    
    type_wrapper.wrap_type("itk::RGBAPixel", "RGBA")
    type_wrapper.add_template(type_wrapper.ITKM["UC"], type_wrapper.ITKT["UC"])
    if WRAP["rgba_unsigned_short"] :
        type_wrapper.add_template(type_wrapper.ITKM["US"], type_wrapper.ITKT["US"])
    type_wrapper.add_template(type_wrapper.ITKM["F"], type_wrapper.ITKT["F"])
    type_wrapper.end_wrap_type()
    
    type_wrapper.wrap_type("std::complex", "C")
    if WRAP["complex_float"] :
        type_wrapper.add_template(type_wrapper.ITKM["F"], type_wrapper.ITKT["F"])
    if WRAP["complex_double"] :
        type_wrapper.add_template(type_wrapper.ITKM["D"], type_wrapper.ITKT["D"])
    type_wrapper.end_wrap_type()
    
    type_wrapper.wrap_type("itk::SymmetricSecondRankTensor", "SSRT")
    for d in WRAP_ITK_DIMS :
        for t in ["D", "F"] :
            type_wrapper.add_template("{0}{1}".format(type_wrapper.ITKM[t], d),
                                      "{0},{1}".format(type_wrapper.ITKT[t], d))
    type_wrapper.end_wrap_type()
    
    type_wrapper.wrap_type("itk::Image", "I")
    for d in WRAP_ITK_DIMS :
        for t in set(WRAP_ITK_ALL_TYPES+["D","UC","UL","RGBUC","RGBAUC"]) :
            if t in WRAP_ITK_VECTOR :
                t = "{0}{1}".format(t,d)
    
            type_wrapper.add_template("{0}{1}".format(type_wrapper.ITKM[t], d),
                                      "{0},{1}".format(type_wrapper.ITKT[t], d))
        type_wrapper.add_template("{0}{1}".format(type_wrapper.ITKM["FAF{0}".format(d)], d),
                                  "{0},{1}".format(type_wrapper.ITKT["FAF{0}".format(d)], d))
        type_wrapper.add_template("{0}{1}".format(type_wrapper.ITKM["O{0}".format(d)], d),
                                  "{0},{1}".format(type_wrapper.ITKT["O{0}".format(d)], d))
        type_wrapper.add_template("{0}{1}".format(type_wrapper.ITKM["SSRT{0}{1}".format("D", d)], d),
                                  "{0},{1}".format(type_wrapper.ITKT["SSRT{0}{1}".format("D", d)], d))
    type_wrapper.end_wrap_type()
    
    type_wrapper.wrap_type("itk::VectorImage", "VI")
    for d in WRAP_ITK_DIMS :
        for t in set(WRAP_ITK_SCALAR+["UC"]) :
            type_wrapper.add_template("{0}{1}".format(type_wrapper.ITKM[t], d),
                                      "{0},{1}".format(type_wrapper.ITKT[t], d))
    type_wrapper.end_wrap_type()
    
    type_wrapper.wrap_type("itk::VariableLengthVector", "VLV")
    for t in set(WRAP_ITK_SCALAR+["UC"]) :
        type_wrapper.add_template(type_wrapper.ITKM[t], type_wrapper.ITKT[t])
    type_wrapper.end_wrap_type()
    
    type_wrapper.wrap_type("itk::Point", "P")
    for d in WRAP_ITK_DIMS :
        for t in set(["F", "D"]) :
            type_wrapper.add_template("{0}{1}".format(type_wrapper.ITKM[t], d),
                                      "{0},{1}".format(type_wrapper.ITKT[t], d))
    type_wrapper.end_wrap_type()
    
    type_wrapper.wrap_type("itk::LevelSetNode", "LSN")
    for d in WRAP_ITK_DIMS :
        for t in WRAP_ITK_SCALAR :
            type_wrapper.add_template("{0}{1}".format(type_wrapper.ITKM[t], d),
                                      "{0},{1}".format(type_wrapper.ITKT[t], d))
    type_wrapper.end_wrap_type()
    
    type_wrapper.wrap_type("itk::FlatStructuringElement", "SE")
    for d in WRAP_ITK_DIMS :
        type_wrapper.add_template(d, d)
    type_wrapper.end_wrap_type()
    
    type_wrapper.wrap_type("itk::SpatialObject", "SO")
    for d in WRAP_ITK_DIMS :
        type_wrapper.add_template(d, d)
    type_wrapper.end_wrap_type()
    
    type_wrapper.wrap_type("itk::Statistics::Histogram", "H")
    for t in ["F", "D"] :
        type_wrapper.add_template(type_wrapper.ITKM[t], type_wrapper.ITKT[t])
    type_wrapper.end_wrap_type()
    
    type_wrapper.wrap_type("itk::LabelMap", "LM")
    for d in WRAP_ITK_DIMS :
        type_wrapper.add_template(d, "itk::StatisticsLabelObject< {0}, {1} >".format(type_wrapper.ITKT["UL"], d))
    type_wrapper.end_wrap_type()

