""" This module contains a DICOM routing engine which applies :class:`rules <Rule>`
    to :class:`datasets <medipy.io.dicom.DataSet>`. Each rule contains 
    :class:`conditions <Condition>` and :class:`actions <Action>`: if all 
    conditions are met, then the actions are applied on the dataset.
    
    The following conditions are available :
        * :class:`AlwaysTrue`
        * :class:`AlwaysFalse`
        * :class:`And`
        * :class:`Or`
        * :class:`Not`
        * :class:`ElementMatch`
    
    The following actions are available :
        * :class:`SetElement`
        * :class:`DeleteElement`
        * :class:`EmptyElement`
        * :class:`ModifyDataSet`
        * :class:`SaveDataSet`
    
    The following example modifies datasets if they match the criteria, and
    saves the modified datasets to a directory. ::
    
        # Get some data sets (e.g. from a PACS)
        datasets = []
        
        condition = medipy.io.dicom.routing.ElementMatch("patients_name", "Doe^John")
        
        actions = [
            medipy.io.dicom.routing.EmptyElement("patients_name"),
            medipy.io.dicom.routing.SetElement("patient_id", "12345"),
            medipy.io.dicom.routing.EmptyElement("patients_birth_date"),
            medipy.io.dicom.routing.SetElement("study_id", "FOO"),
            medipy.io.dicom.routing.SaveDataSet("/some/where", "hierarchical"),
        ]
        
        rule = medipy.io.dicom.routing.Rule(condition, actions)
        
        for dataset in datasets :
            rule(dataset)
"""

from actions import (Action, SetElement, DeleteElement, EmptyElement, 
                     ModifyDataSet, SaveDataSet)
from conditions import (Condition, AlwaysTrue, AlwaysFalse, And, Or, Not,
                        ElementMatch)
from Rule import Rule

__all__ = [
    "Action", "SetElement", "DeleteElement", "EmptyElement", "ModifyDataSet", 
    "SaveDataSet", "Condition", "AlwaysTrue", "AlwaysFalse", "And", "Or", "Not",
    "ElementMatch", "Rule"
]