##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

class VR(object):
    """ Base class for all VR classes.
    """
    
    def __init__(self, value):
        self.value = value
    
    def __eq__(self, other):
        return self.value == other
    
    def __ne__(self, other):
        return self.value != other

class AE(VR): pass
class AS(VR): pass
class AT(VR): pass
class CS(VR): pass
class DA(VR): pass
class DS(VR): pass
class DT(VR): pass
class FL(VR): pass
class FD(VR): pass
class IS(VR): pass

class LO(VR): pass
class LT(VR): pass
class OB(VR): pass
class OF(VR): pass
class OW(VR): pass
class PN(VR): pass
class SH(VR): pass
class SL(VR): pass
class SQ(VR): pass
class SS(VR): pass

class UI(VR): pass
class TM(VR): pass
class ST(VR): pass
class UL(VR): pass
class UN(VR): pass
class US(VR): pass
class UT(VR): pass
