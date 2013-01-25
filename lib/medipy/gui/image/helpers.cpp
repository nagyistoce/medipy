#include <Python.h>

#include <iostream>
#include <PyVTKObject.h>
#include <vtkImageChangeInformation.h>
#include <vtkImageData.h>
#include <vtkImageReslice.h>
#include <vtkMath.h>
#include <vtkMatrix4x4.h>

void sequence_to_array(PyObject* sequence, double* array, bool reverse=false)
{
    Py_ssize_t const size = PySequence_Size(sequence);
    for(unsigned int i=0; i<size; ++i)
    {
        PyObject* object = PySequence_GetItem(sequence, i);
        double const value = PyFloat_AsDouble(object);
        Py_DECREF(object);

        unsigned int const index = reverse?(size-i-1):i;
        array[index] = value;

    }
}

void index_to_physical(double* index, double* origin, double* spacing, double* physical)
{
    for(unsigned int i=0; i<3; ++i)
    {
        physical[i] = index[i]*spacing[i]+origin[i];
    }
}

void physical_to_index(double* physical, double* origin, double* spacing, double* index)
{
    for(unsigned int i=0; i<3; ++i)
    {
        index[i] = (physical[i]-origin[i])/spacing[i];
    }
}

PyObject* layer_index_to_world(PyObject * layer, PyObject * index_python)
{
    // Get VTK objects for the Layer class
    vtkImageChangeInformation * change_information =
        reinterpret_cast<vtkImageChangeInformation*>(reinterpret_cast<PyVTKObject*>(
            PyObject_GetAttrString(layer, "_change_information"))->vtk_ptr);
    vtkImageReslice * reslicer =
                reinterpret_cast<vtkImageReslice*>(reinterpret_cast<PyVTKObject*>(
                    PyObject_GetAttrString(layer, "_reslicer"))->vtk_ptr);
    vtkImageData * vtk_image =
        reinterpret_cast<vtkImageData*>(reinterpret_cast<PyVTKObject*>(
            PyObject_GetAttrString(layer, "_vtk_image"))->vtk_ptr);
    vtkMatrix4x4 * reslicer_axes_inverse =
        reinterpret_cast<vtkMatrix4x4*>(reinterpret_cast<PyVTKObject*>(
            PyObject_GetAttrString(layer, "_reslicer_axes_inverse"))->vtk_ptr);

    // Get the index coordinates as a C object, converting to VTK order
    double index[3];
    sequence_to_array(index_python, index, true);

    // Make sure the pipeline is up-to-date
    change_information->Update();

    // Convert from the non-sliced image point coordinates (VTK order)
    double physical[4];
    index_to_physical(index,
        vtk_image->GetOrigin(), vtk_image->GetSpacing(), physical);

    // Apply the inverse reslicer transform (homogeneous coordinates, VTK
    // order), converting to the sliced image
    physical[3] = 1;
    double transformed[4];
    reslicer_axes_inverse->MultiplyPoint(physical, transformed);
    vtkMath::MultiplyScalar(transformed, 1./transformed[3]);

    // Convert to index coordinate in resliced image (VTK order)
    double index_resliced[3];
    physical_to_index(transformed,
        reslicer->GetOutput()->GetOrigin(), reslicer->GetOutput()->GetSpacing(),
        index_resliced);

    // Convert to world coordinates
    double world[3];
    index_to_physical(index_resliced,
        change_information->GetOutputOrigin(), change_information->GetOutputSpacing(),
        world);

    return Py_BuildValue("(d,d,d)", world[0], world[1], world[2]);
}

PyObject* layer_world_to_index(PyObject * layer, PyObject * world_python)
{
    // Get VTK objects for the Layer class
    vtkImageChangeInformation * change_information =
        reinterpret_cast<vtkImageChangeInformation*>(reinterpret_cast<PyVTKObject*>(
            PyObject_GetAttrString(layer, "_change_information"))->vtk_ptr);
    vtkImageReslice * reslicer =
            reinterpret_cast<vtkImageReslice*>(reinterpret_cast<PyVTKObject*>(
                PyObject_GetAttrString(layer, "_reslicer"))->vtk_ptr);
    vtkImageData * vtk_image =
        reinterpret_cast<vtkImageData*>(reinterpret_cast<PyVTKObject*>(
            PyObject_GetAttrString(layer, "_vtk_image"))->vtk_ptr);

    // Get the world coordinates as a C object
    double world[3];
    sequence_to_array(world_python, world);

    // Make sure the pipeline is up-to-date
    change_information->Update();

    // Convert to index coordinate in resliced image (VTK order)
    double index[3];
    physical_to_index(world,
        change_information->GetOutputOrigin(), change_information->GetOutputSpacing(),
        index);

    // Set height to 0, since the picked value will depend on the position of the actor
    index[2] = 0;

    // Apply the reslicer transform (homogeneous coordinates, VTK order),
    // converting to the non-sliced image
    double physical[4];
    index_to_physical(index,
        reslicer->GetOutput()->GetOrigin(), reslicer->GetOutput()->GetSpacing(),
        physical);
    physical[3] = 1;

    double transformed[4];
    reslicer->GetResliceAxes()->MultiplyPoint(physical, transformed);
    vtkMath::MultiplyScalar(transformed, 1./transformed[3]);

    // Convert to index coordinate in the non-sliced image (VTK order)
    physical_to_index(transformed,
        vtk_image->GetOrigin(), vtk_image->GetSpacing(), index);

    // VTK order -> NumPy order
    return Py_BuildValue("(d,d,d)", index[2], index[1], index[0]);
}
