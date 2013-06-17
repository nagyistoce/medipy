#ifndef _ed3c3514_93e1_46f7_993b_84fd6bbafd9f
#define _ed3c3514_93e1_46f7_993b_84fd6bbafd9f

/**
 * @brief Convert a VTK world coordinate (VTK order) to the corresponding image
 *        index (numpy order).
 *
 * @param layer an instance of medipy.gui.image.Layer
 * @param index a sequence of 3 floats, in numpy (i.e. z,y,x) order
 */
PyObject* layer_index_to_world(PyObject * layer, PyObject * index);

/**
 * @brief Convert a VTK world coordinate (VTK order) to the corresponding image
 *        index (numpy order).
 *
 * @param layer an instance of medipy.gui.image.Layer
 * @param world a sequence of 3 floats, in VTK (i.e. x,y,z) order
 */
PyObject* layer_world_to_index(PyObject * layer, PyObject * world);

#endif // _ed3c3514_93e1_46f7_993b_84fd6bbafd9f
