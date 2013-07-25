##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import weakref
from vtk import (vtkCellLocator, vtkClipPolyData, vtkLODActor, vtkMaskPoints, 
                 vtkPointLocator, vtkPolyData, vtkPolyDataMapper)

class Object3D(object):
    """ Represent a 3D object, based on a vtkPolyData
    """
    def __init__(self, dataset=None, name=None, image=None, gui_image=None):
        ####################
        # Public interface #
        ####################
        self.name = name
        
        #####################
        # Private interface #
        #####################
        self._dataset = None
        self._actor = vtkLODActor()
        self._image = None
        self._gui_image = None
        
        self._point_locator = vtkPointLocator()
        self._cell_locator = vtkCellLocator()
        self._point_locator_dirty = True
        self._cell_locator_dirty = True
        
        self._mapper = vtkPolyDataMapper()
        self._clipping_mapper = vtkPolyDataMapper()
        self._inside_out = True
        self._clipping_functions = None
        self._clippers = []
        self._mask_filter = vtkMaskPoints()
        
        ##################
        # Initialization #
        ##################
        self._mask_filter.SetOnRatio(1)
        self._mask_filter.SetOffset(0)
        self._mask_filter.GenerateVerticesOn()
        
        self._actor.SetMapper(self._mapper)
        self._set_dataset(dataset)
        self.color_by_material()
        self._set_image(image)
        self._set_gui_image(image)
            
        
    ####################
    # Public interface #
    ####################
#        def get_point_array(self, array):
#            return self._dataset.GetPointData().GetArray(array)
#        
#        def get_cell_array(self, array):
#            return self._dataset.GetCellData().GetArray(array)
#        
#        def get_field_array(self, array):
#            return self._dataset.GetFieldData().GetArray(array)
#        
#        def set_scalar_range(self, range):
#            self._mapper.SetScalarRange(*range)
        
    def color_by_material(self):
        self._mapper.ScalarVisibilityOff()
        self._clipping_mapper.ScalarVisibilityOff()
        
            
    def color_by_scalars(self, lut, min, max):
        for mapper in (self._mapper, self._clipping_mapper) :
            mapper.SetLookupTable(lut)
            mapper.SetColorModeToMapScalars()        
            mapper.SetScalarRange(min, max)
            mapper.ScalarVisibilityOn()

    def color_by_point_scalars(self,array_name):
        self._mapper.SetScalarModeToUsePointData()
        self._dataset.GetPointData().SetActiveScalars(array_name)
        self._mapper.ScalarVisibilityOn()

    def point_arrays(self):
        result = []
        for i in range(self._dataset.GetPointData().GetNumberOfArrays()) :
            result.append(self._dataset.GetPointData().GetArray(i).GetName())
        return result
        
        
           
#            self._mapper.SelectColorArray(1)
#            self._mapper.InterpolateScalarsBeforeMappingOff()
#        
#        def color_by_cell_scalars(self, array):
#            self._mapper.SetScalarModeToUseCellData()
#            self._dataset.GetCellData().SetActiveScalars(array)
#            self._mapper.ScalarVisibilityOn()
#            
#            self._mapper.SelectColorArray(1)
#            self._mapper.SetScalarRange(0, 6)
#            self._mapper.InterpolateScalarsBeforeMappingOff()
#        
#        
#        def cell_arrays(self):
#            result = []
#            for i in range(self._dataset.GetCellData().GetNumberOfArrays()) :
#                result.append(self._dataset.GetCellData().GetArray(i).GetName())
#            return result
#        
#        def field_arrays(self):
#            result = []
#            for i in range(self._dataset.GetFieldData().GetNumberOfArrays()) :
#                result.append(self._dataset.GetFieldData().GetArray(i).GetName())
#            return result
        
#        def show_vertices(self):
#            self._mapper.SetInputConnection(self._mask_filter.GetOutputPort())
#            self._mapper.Modified()
#        
#        def hide_vertices(self):
#            self._mapper.SetInput(self._dataset)
#            self._mapper.Modified()
        
    ##############
    # Properties #
    ##############
    def _set_dataset(self, dataset):
        if dataset is not None :
            self._dataset = dataset
        else : 
            self._dataset = vtkPolyData()
        self._point_locator.SetDataSet(dataset)
        self._cell_locator.SetDataSet(dataset)
        
        self._point_locator_dirty = True
        self._cell_locator_dirty = True
        
        self._mask_filter.SetInput(self._dataset)
        self._mapper.SetInput(self._dataset)
    
    def _get_color(self):
        return self._actor.GetProperty().GetColor()
    
    def _set_color(self, color):
        return self._actor.GetProperty().SetColor(color)
    
    def _get_diffuse_color(self):
        return self._actor.GetProperty().GetDiffuseColor()
    
    def _set_diffuse_color(self, color):
        self._actor.GetProperty().SetDiffuseColor(color)
    
    def _get_ambient_color(self):
        return self._actor.GetProperty().GetAmbientColor()
    
    def _set_ambient_color(self, color):
        self._actor.GetProperty().SetAmbientColor(color)
    
    def _get_specular_color(self):
        return self._actor.GetProperty().GetSpecularColor()
    
    def _set_specular_color(self, color):
        self._actor.GetProperty().SetSpecularColor(color)
        
    def _get_opacity(self):
        return self._actor.GetProperty().GetOpacity()
    
    def _set_opacity(self, opacity):
        self._actor.GetProperty().SetOpacity(opacity)
    
    def _get_visibility(self):
        return self._actor.GetVisibility()
    
    def _set_visibility(self, visibility):
        self._actor.SetVisibility(visibility)
        
    def _get_representation(self):
        return self._actor.GetProperty().GetRepresentationAsString()
    
    def _set_representation(self, representation):
        getattr(self._actor.GetProperty(), "SetRepresentationTo" + representation)()
        
    def _get_points_size(self):
        return self._actor.GetProperty().GetPointSize()
    
    def _set_points_size(self, points_size):
        self._actor.GetProperty().SetPointSize(points_size)    
        
    def _get_shading(self):
        return self._actor.GetProperty().GetInterpolationAsString()
    
    def _set_shading(self, shading):
        getattr(self._actor.GetProperty(), "SetInterpolationTo" + shading)()
    
    def _get_image(self):
        if self._image is not None :
            return self._image()
        else :
            return None
    
    def _set_image(self, image):
        if image is not None :
            self._image = weakref.ref(image)
        else : 
            self._image = None
    
    def _get_gui_image(self):
        if self._gui_image is not None :
            return self._gui_image()
        else :
            return None
    
    def _set_gui_image(self, image):
        if image is not None :
            self._gui_image = weakref.ref(image)
        else : 
            self._gui_image = None
    
    def _set_inside_out(self, flag):
        self._inside_out = flag
        for clipper in self._clippers :
            clipper.SetInsideOut(self._inside_out)
    
    def _set_clipping_functions(self, functions):
        if functions :
            self._build_clippers(functions)
            self._actor.SetMapper(self._clipping_mapper)
        else :
            self._actor.SetMapper(self._mapper)
        self._clipping_functions = functions
    
    def _get_point_locator(self):
        if self._point_locator_dirty :
            self._point_locator.BuildLocator()
            self._point_locator_dirty = False
        return self._point_locator
    
    def _get_cell_locator(self):
        if self._cell_locator_dirty :
            self._cell_locator.BuildLocator()
            self._cell_locator_dirty = False
        return self._cell_locator
        
    actor = property(lambda x:x._actor)
    mapper = property(lambda x:x._actor.GetMapper())
    dataset = property(lambda x:x._dataset, _set_dataset)
    diffuse_color = property(_get_diffuse_color, _set_diffuse_color)
    ambient_color = property(_get_ambient_color, _set_ambient_color)
    specular_color = property(_get_specular_color, _set_specular_color)
    color = property(_get_color, _set_color)
    opacity = property(_get_opacity, _set_opacity)
    visibility = property(_get_visibility, _set_visibility)
    representation = property(_get_representation, _set_representation)
    points_size = property(_get_points_size, _set_points_size)
    shading = property(_get_shading, _set_shading)
    inside_out = property(lambda x:x._inside_out, _set_inside_out)
    clipping_functions = property(lambda x:x._clipping_functions, _set_clipping_functions)
    image = property(_get_image, _set_image)
    gui_image = property(_get_gui_image, _set_gui_image)
    point_locator = property(_get_point_locator)
    cell_locator = property(_get_cell_locator)

    #####################
    # Private interface #
    #####################
    
    def _build_clippers(self, functions):
        input = self.dataset
        
        for function in functions : 
            clipper = vtkClipPolyData()
            clipper.SetClipFunction(function)
            clipper.GenerateClipScalarsOff()
            clipper.GenerateClippedOutputOff()
            clipper.SetInsideOut(self._inside_out)
            clipper.SetValue(0)
            clipper.SetInput(input)
            self._clippers.append(clipper)
            input = clipper.GetOutput()

        self._clipping_mapper.SetInput(input)
