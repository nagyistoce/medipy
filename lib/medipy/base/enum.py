##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

def enum(name, *members, **named):
    """ Create an enumerated type called ``name`` with given ``members``: ::
        
            >>> Colors = medipy.base.enum("Colors", "red", "green", "blue")
            >>> color = Colors.red
        
        By default, the members are integers, but their value can be specified
        using named parameters. In this case, the values must be storable in a
        dictionary: ::
        
            >>> Colors = medipy.base.enum("Colors", "red", green="green", blue=3.14)
            >>> Colors.red
            0
            >>> Colors.green
            'green'
            >>> Colors.blue
            3.14
        
        The name of members can be retrieved either in bulk or using their 
        values: ::
        
            >>> Colors.members
            ['red', 'green', 'blue']
            >>> Colors[0]
            'red'
            >>> Colors["green"]
            'green'
            >>> Colors[3.14]
            'blue'
        
        The enumerated type is derived of :class:`~medipy.base.Enum`: ::
        
            >>> isinstance(Colors, medipy.base.Enum)
            True
    """
    
    enums = dict(zip(members, range(len(members))), **named)
    reverse = dict((value, key) for key, value in enums.iteritems())
    new_type = Enum(name, (), enums)
    new_type._reverse = reverse
    return new_type

class Enum(type) :
    """ Base class for all enumerated types. This class should not be used
        directly
    """
    
    def __init__(cls, name, bases, dct):
        type.__init__(cls, name, bases, dct)
    
    @property
    def members(self) :
        return self._reverse.values()
    
    def __getitem__(self, item) :
        """ Return the name of a member given its value.
        """
        
        return self._reverse[item]
