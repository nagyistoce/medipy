##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import wx
import wx.xrc as xrc
import wx.gizmos as gizmos

class TreeListCtrlXmlHandler(xrc.XmlResourceHandler):
    def __init__(self):
        xrc.XmlResourceHandler.__init__(self)
        # Standard styles
        self.AddWindowStyles()
        # Custom styles
        self.AddStyle('wxDEFAULT_COL_WIDTH', gizmos.DEFAULT_COL_WIDTH)
        self.AddStyle('wxTL_MODE_NAV_FULLTREE', gizmos.TL_MODE_NAV_FULLTREE)
        self.AddStyle('wxTL_MODE_NAV_EXPANDED', gizmos.TL_MODE_NAV_EXPANDED)
        self.AddStyle('wxTL_MODE_NAV_VISIBLE', gizmos.TL_MODE_NAV_VISIBLE)
        self.AddStyle('wxTL_MODE_NAV_LEVEL', gizmos.TL_MODE_NAV_LEVEL)
        self.AddStyle('wxTL_MODE_FIND_EXACT', gizmos.TL_MODE_FIND_EXACT)
        self.AddStyle('wxTL_MODE_FIND_PARTIAL', gizmos.TL_MODE_FIND_PARTIAL)
        self.AddStyle('wxTL_MODE_FIND_NOCASE', gizmos.TL_MODE_FIND_NOCASE)
        self.AddStyle('wxTREE_HITTEST_ONITEMCOLUMN', gizmos.TREE_HITTEST_ONITEMCOLUMN)
        self.AddStyle('wxTR_COLUMN_LINES', gizmos.TR_COLUMN_LINES)
        self.AddStyle('wxTR_VIRTUAL', gizmos.TR_VIRTUAL)
        self.AddStyle('wxTL_ALIGN_LEFT  ', wx.ALIGN_LEFT)
        self.AddStyle('wxTL_ALIGN_RIGHT ', wx.ALIGN_RIGHT)
        self.AddStyle('wxTL_ALIGN_CENTER', wx.ALIGN_CENTER)

        self.AddStyle('wxTL_SEARCH_VISIBLE', gizmos.TL_MODE_NAV_VISIBLE)
        self.AddStyle('wxTL_SEARCH_LEVEL  ', gizmos.TL_MODE_NAV_LEVEL)
        self.AddStyle('wxTL_SEARCH_FULL   ', gizmos.TL_MODE_FIND_EXACT)
        self.AddStyle('wxTL_SEARCH_PARTIAL', gizmos.TL_MODE_FIND_PARTIAL)
        self.AddStyle('wxTL_SEARCH_NOCASE ', gizmos.TL_MODE_FIND_NOCASE)

        self.AddStyle('wxTR_DONT_ADJUST_MAC', gizmos.TR_DONT_ADJUST_MAC)
        self.AddStyle('wxTR_DEFAULT_STYLE', wx.TR_DEFAULT_STYLE)
        
        self.AddStyle('wxTR_HIDE_ROOT', wx.TR_HIDE_ROOT)
        self.AddStyle('wxTR_FULL_ROW_HIGHLIGHT',wx.TR_FULL_ROW_HIGHLIGHT)
        self.AddStyle('wxTR_EXTENDED',wx.TR_EXTENDED)
        self.AddStyle('wxTR_MULTIPLE',wx.TR_MULTIPLE)

    def CanHandle(self, node):
        return self.IsOfClass(node, 'TreeListCtrl')

    # Process XML parameters and create the object
    def DoCreateResource(self):
        assert self.GetInstance() is None

        w = gizmos.TreeListCtrl(self.GetParentAsWindow(),
                                self.GetID(),
                                style=self.GetStyle(),
                                name=self.GetName())
        return w
