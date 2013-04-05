import unittest

import numpy
import medipy.base
import medipy.io.image_annotation

class TestImageAnnotation(unittest.TestCase):
    def test_xml(self) :
        annotation = medipy.base.ImageAnnotation(
            (1, 2.718, 3.14), "Something", 
            medipy.base.ImageAnnotation.Shape.cube, 10, (0.5,0.2,0.3), True, 
            "lorem\nipsum")
        data = medipy.io.image_annotation.annotation_to_xml(annotation)
        other_annotation = medipy.io.image_annotation.annotation_from_xml(data)
        
        numpy.testing.assert_equal(annotation.position, other_annotation.position)
        self.assertEqual(annotation.label, other_annotation.label)
        self.assertEqual(annotation.shape, other_annotation.shape)
        self.assertEqual(annotation.size, other_annotation.size)
        numpy.testing.assert_equal(annotation.color, other_annotation.color)
        self.assertEqual(annotation.filled, other_annotation.filled)
        self.assertEqual(annotation.comment, other_annotation.comment)
    
    def test_list_xml(self) :
        annotation1 = medipy.base.ImageAnnotation(
            (1, 2.718, 3.14), "Something", 
            medipy.base.ImageAnnotation.Shape.cube, 10, (0.5,0.2,0.3), True, 
            "lorem\nipsum")
        annotation2 = medipy.base.ImageAnnotation(
            (42, 1.2345, 6.789), "Nothing", 
            medipy.base.ImageAnnotation.Shape.cross, 20, (0.1,0.9,0.7), False, 
            "ipsum\nlorem")
        
        annotations = [annotation1, annotation2]
        data = medipy.io.image_annotation.annotations_to_xml(annotations)
        other_annotations = medipy.io.image_annotation.annotations_from_xml(data)
        
        self.assertEqual(len(annotations), len(other_annotations))
        for index, annotation in enumerate(annotations) :
            other_annotation = other_annotations[index] 
            numpy.testing.assert_equal(annotation.position, other_annotation.position)
            self.assertEqual(annotation.label, other_annotation.label)
            self.assertEqual(annotation.shape, other_annotation.shape)
            self.assertEqual(annotation.size, other_annotation.size)
            numpy.testing.assert_equal(annotation.color, other_annotation.color)
            self.assertEqual(annotation.filled, other_annotation.filled)
            self.assertEqual(annotation.comment, other_annotation.comment)

if __name__ == '__main__':
    unittest.main()