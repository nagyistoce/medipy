import os
import unittest

import numpy

import medipy.base
import medipy.diffusion

class TestTractography(unittest.TestCase) :
    
    def setUp(self) :
        self.tensors = medipy.base.Image(shape=(11,11,11), value=0, dti="tensor_2")
        for i in xrange(11) :
            for j in xrange(-1, 2) :
                for k in xrange(-1, 2) :
                    self.tensors[i, 5+j, 5+k] = [0.2, 0, 0, 0.2, 0, 1]
                    self.tensors[5+j, i, 5+k] = [0.2, 0, 0, 1, 0, 0.2]
                    self.tensors[5+j, 5+k, i] = [1, 0, 0, 0.2, 0, 0.2]
        self.tensors.data[4:7,4:7,4:7] = 0
        self.tensors.data *= 1e-3
        
        self.fibers = [
            [[1.0, 4.0, 4.0], [1.5, 4.0, 4.0], [2.0, 4.0, 4.0], [2.5, 4.0, 4.0], 
             [3.0, 4.0, 4.0], [3.5, 4.0, 4.0]], 
            [[1.0, 4.0, 6.0], [1.5, 4.0, 6.0], [2.0, 4.0, 6.0], [2.5, 4.0, 6.0], 
             [3.0, 4.0, 6.0], [3.5, 4.0, 6.0]], 
            [[1.0, 6.0, 4.0], [1.5, 6.0, 4.0], [2.0, 6.0, 4.0], [2.5, 6.0, 4.0], 
             [3.0, 6.0, 4.0], [3.5, 6.0, 4.0]], 
            [[1.0, 6.0, 6.0], [1.5, 6.0, 6.0], [2.0, 6.0, 6.0], [2.5, 6.0, 6.0], 
             [3.0, 6.0, 6.0], [3.5, 6.0, 6.0]], 
            [[4.0, 1.0, 4.0], [4.0, 1.5, 4.0], [4.0, 2.0, 4.0], [4.0, 2.5, 4.0], 
             [4.0, 3.0, 4.0], [4.0, 3.5, 4.0]], 
            [[4.0, 1.0, 6.0], [4.0, 1.5, 6.0], [4.0, 2.0, 6.0], [4.0, 2.5, 6.0], 
             [4.0, 3.0, 6.0], [4.0, 3.5, 6.0]], 
            [[4.0, 4.0, 1.0], [4.0, 4.0, 1.5], [4.0, 4.0, 2.0], [4.0, 4.0, 2.5], 
             [4.0, 4.0, 3.0], [4.0, 4.0, 3.5]], 
            [[4.0, 4.0, 6.5], [4.0, 4.0, 7.0], [4.0, 4.0, 7.5], [4.0, 4.0, 8.0], 
             [4.0, 4.0, 8.5], [4.0, 4.0, 9.0], [4.0, 4.0, 9.5]], 
            [[4.0, 6.0, 1.0], [4.0, 6.0, 1.5], [4.0, 6.0, 2.0], [4.0, 6.0, 2.5], 
             [4.0, 6.0, 3.0], [4.0, 6.0, 3.5]], 
            [[4.0, 6.0, 6.5], [4.0, 6.0, 7.0], [4.0, 6.0, 7.5], [4.0, 6.0, 8.0], 
             [4.0, 6.0, 8.5], [4.0, 6.0, 9.0], [4.0, 6.0, 9.5]], 
            [[4.0, 6.5, 4.0], [4.0, 7.0, 4.0], [4.0, 7.5, 4.0], [4.0, 8.0, 4.0], 
             [4.0, 8.5, 4.0], [4.0, 9.0, 4.0], [4.0, 9.5, 4.0]], 
            [[4.0, 6.5, 6.0], [4.0, 7.0, 6.0], [4.0, 7.5, 6.0], [4.0, 8.0, 6.0], 
             [4.0, 8.5, 6.0], [4.0, 9.0, 6.0], [4.0, 9.5, 6.0]], 
            [[6.0, 1.0, 4.0], [6.0, 1.5, 4.0], [6.0, 2.0, 4.0], [6.0, 2.5, 4.0], 
             [6.0, 3.0, 4.0], [6.0, 3.5, 4.0]], 
            [[6.0, 1.0, 6.0], [6.0, 1.5, 6.0], [6.0, 2.0, 6.0], [6.0, 2.5, 6.0], 
             [6.0, 3.0, 6.0], [6.0, 3.5, 6.0]], 
            [[6.0, 4.0, 1.0], [6.0, 4.0, 1.5], [6.0, 4.0, 2.0], [6.0, 4.0, 2.5], 
             [6.0, 4.0, 3.0], [6.0, 4.0, 3.5]], 
            [[6.0, 4.0, 6.5], [6.0, 4.0, 7.0], [6.0, 4.0, 7.5], [6.0, 4.0, 8.0], 
             [6.0, 4.0, 8.5], [6.0, 4.0, 9.0], [6.0, 4.0, 9.5]], 
            [[6.0, 6.0, 1.0], [6.0, 6.0, 1.5], [6.0, 6.0, 2.0], [6.0, 6.0, 2.5], 
             [6.0, 6.0, 3.0], [6.0, 6.0, 3.5]], 
            [[6.0, 6.0, 6.5], [6.0, 6.0, 7.0], [6.0, 6.0, 7.5], [6.0, 6.0, 8.0], 
             [6.0, 6.0, 8.5], [6.0, 6.0, 9.0], [6.0, 6.0, 9.5]], 
            [[6.0, 6.5, 4.0], [6.0, 7.0, 4.0], [6.0, 7.5, 4.0], [6.0, 8.0, 4.0], 
             [6.0, 8.5, 4.0], [6.0, 9.0, 4.0], [6.0, 9.5, 4.0]], 
            [[6.0, 6.5, 6.0], [6.0, 7.0, 6.0], [6.0, 7.5, 6.0], [6.0, 8.0, 6.0], 
             [6.0, 8.5, 6.0], [6.0, 9.0, 6.0], [6.0, 9.5, 6.0]], 
            [[6.5, 4.0, 4.0], [7.0, 4.0, 4.0], [7.5, 4.0, 4.0], [8.0, 4.0, 4.0], 
             [8.5, 4.0, 4.0], [9.0, 4.0, 4.0], [9.5, 4.0, 4.0]], 
            [[6.5, 4.0, 6.0], [7.0, 4.0, 6.0], [7.5, 4.0, 6.0], [8.0, 4.0, 6.0], 
             [8.5, 4.0, 6.0], [9.0, 4.0, 6.0], [9.5, 4.0, 6.0]], 
            [[6.5, 6.0, 4.0], [7.0, 6.0, 4.0], [7.5, 6.0, 4.0], [8.0, 6.0, 4.0], 
             [8.5, 6.0, 4.0], [9.0, 6.0, 4.0], [9.5, 6.0, 4.0]], 
            [[6.5, 6.0, 6.0], [7.0, 6.0, 6.0], [7.5, 6.0, 6.0], [8.0, 6.0, 6.0], 
             [8.5, 6.0, 6.0], [9.0, 6.0, 6.0], [9.5, 6.0, 6.0]]]

        self.clusters = {
            0: {'indices': [0, 20], 'hidden': numpy.asarray(
                [[7.5, 8., 8.], [10.25, 8., 8.], [13., 8., 8.]]), 'N': 2,
                'skel': 0}, 
            1: {'indices': [1, 21], 'hidden': numpy.asarray(
                [[7.5, 8., 12.], [10.25, 8., 12.], [13., 8., 12.]]), 'N': 2,
                'skel': 1}, 
            2: {'indices': [2, 22], 'hidden': numpy.asarray(
                [[7.5, 12., 8.], [10.25, 12.,  8.], [13., 12., 8.]]), 'N': 2,
                'skel': 2}, 
            3: {'indices': [3, 23], 'hidden': numpy.asarray(
                [[7.5, 12., 12.], [10.25, 12., 12.], [13., 12., 12.]]), 'N': 2,
                'skel': 3}, 
            4: {'indices': [4, 12], 'hidden': numpy.asarray(
                [[10., 2., 8.], [10.,  4.5,  8. ], [10. ,  7., 8.]]), 'N': 2,
                'skel': 4}, 
            5: {'indices': [5, 13], 'hidden': numpy.asarray(
                [[10., 2., 12.], [10.,  4.5, 12. ], [10. ,  7., 12.]]), 'N': 2,
                'skel': 5}, 
            6: {'indices': [6, 14], 'hidden': numpy.asarray(
                [[10., 8., 2.], [10.,  8. ,  4.5], [10. ,  8., 7.]]), 'N': 2,
                'skel': 6}, 
            7: {'indices': [7, 15], 'hidden': numpy.asarray(
                [[10., 8., 13.], [10.,  8., 16.], [10.,  8., 19.]]), 'N': 2,
                'skel': 7}, 
            8: {'indices': [8, 16], 'hidden': numpy.asarray(
                [[10., 12., 2.], [10., 12. ,  4.5], [10. , 12., 7.]]), 'N': 2,
                'skel': 8}, 
            9: {'indices': [9, 17], 'hidden': numpy.asarray(
                [[10., 12., 13.], [10., 12., 16.], [10., 12., 19.]]), 'N': 2,
                'skel': 9}, 
            10: {'indices': [10, 18], 'hidden': numpy.asarray(
                [[10., 13., 8.], [10., 16.,  8.], [10., 19., 8.]]), 'N': 2,
                'skel': 10}, 
            11: {'indices': [11, 19], 'hidden': numpy.asarray(
                [[10., 13., 12.], [10., 16., 12.], [10., 19., 12.]]), 'N': 2,
                'skel': 11}
        }

    
    def test_streamline(self) :
        fibers = medipy.diffusion.tractography.streamline(self.tensors,
            0.5, 0.2, numpy.pi/3, 1, "Euler", None, None)
        
        distinct_fibers = []
        for fiber in fibers :
            if fiber.tolist() not in distinct_fibers :
                distinct_fibers.append(fiber.tolist())
        
        self.assertEqual(self.fibers, distinct_fibers)
    
    def test_local_skeleton_clustering(self) :
        clusters = medipy.diffusion.tractography.local_skeleton_clustering(
            self.fibers, 1)
        
        for key, cluster in self.clusters.items() :
            self.assertTrue(key in clusters)
            self.assertEqual(cluster["indices"], clusters[key]["indices"])
            numpy.testing.assert_array_equal(
                cluster["hidden"], clusters[key]["hidden"])
            self.assertEqual(cluster["N"], clusters[key]["N"])
        for key in clusters :
            self.assertTrue(key in self.clusters)
    
    def test_skeleton(self) :
        clusters = medipy.diffusion.tractography.local_skeleton_clustering(
            self.fibers, 1)
        medipy.diffusion.tractography.skeleton(self.fibers, clusters)
        
        for key, cluster in self.clusters.items() :
            self.assertEqual(cluster["skel"], clusters[key]["skel"])
        for key in clusters :
            self.assertTrue(key in self.clusters)
    
    def test_create_object_3d_without_clusters(self) :
        object_3d = medipy.diffusion.tractography.create_object_3d(
            [numpy.asarray(f) for f in self.fibers])
        self.assertEqual(object_3d.point_arrays(), [])
    
    def test_create_object_3d_with_clusters(self) :
        clusters = medipy.diffusion.tractography.local_skeleton_clustering(
            self.fibers, 1)
        medipy.diffusion.tractography.skeleton(self.fibers, clusters)
        
        object_3d = medipy.diffusion.tractography.create_object_3d(
            [numpy.asarray(f) for f in self.fibers], clusters)
        self.assertEqual(object_3d.point_arrays(), ["cluster", "skeleton"])

if __name__ == "__main__" :
    unittest.main()
