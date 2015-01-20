import os
import unittest

import numpy

import medipy.base
import medipy.diffusion

class TestTractography(unittest.TestCase) :
    
    def setUp(self) :
        # Generate a cross pattern in the xy plane
        self.tensors = medipy.base.Image(shape=(3,11,11), value=0, dti="tensor_2")
        
        for i in range(-1, 2):
            for j in range(-4, 5):
                self.tensors[1, 5+i, 5+j] = [1, 0, 0, 0.2, 0, 0.2]
                self.tensors[1, 5+j, 5+i] = [0.2, 0, 0, 1, 0, 0.2]

        self.tensors.data[1,4:7,4:7] = 0
        self.tensors.data *= 1e-3
        
        self.fibers = [
            [[1.0, 1.0, 4.0], [1.0, 1.5, 4.0], [1.0, 2.0, 4.0], [1.0, 2.5, 4.0], 
                [1.0, 3.0, 4.0], [1.0, 3.5, 4.0]],
            [[1.0, 1.0, 5.0], [1.0, 1.5, 5.0], [1.0, 2.0, 5.0], [1.0, 2.5, 5.0], 
                [1.0, 3.0, 5.0], [1.0, 3.5, 5.0]],
            [[1.0, 1.0, 6.0], [1.0, 1.5, 6.0], [1.0, 2.0, 6.0], [1.0, 2.5, 6.0], 
                [1.0, 3.0, 6.0], [1.0, 3.5, 6.0]],
            [[1.0, 4.0, 1.0], [1.0, 4.0, 1.5], [1.0, 4.0, 2.0], [1.0, 4.0, 2.5], 
                [1.0, 4.0, 3.0], [1.0, 4.0, 3.5]],
            [[1.0, 4.0, 7.0], [1.0, 4.0, 7.5], [1.0, 4.0, 8.0], [1.0, 4.0, 8.5], 
                [1.0, 4.0, 9.0]],
            [[1.0, 5.0, 1.0], [1.0, 5.0, 1.5], [1.0, 5.0, 2.0], [1.0, 5.0, 2.5], 
                [1.0, 5.0, 3.0], [1.0, 5.0, 3.5]],
            [[1.0, 5.0, 7.0], [1.0, 5.0, 7.5], [1.0, 5.0, 8.0], [1.0, 5.0, 8.5], 
                [1.0, 5.0, 9.0]],
            [[1.0, 6.0, 1.0], [1.0, 6.0, 1.5], [1.0, 6.0, 2.0], [1.0, 6.0, 2.5], 
                [1.0, 6.0, 3.0], [1.0, 6.0, 3.5]],
            [[1.0, 6.0, 7.0], [1.0, 6.0, 7.5], [1.0, 6.0, 8.0], [1.0, 6.0, 8.5], 
                [1.0, 6.0, 9.0]],
            [[1.0, 7.0, 4.0], [1.0, 7.5, 4.0], [1.0, 8.0, 4.0], [1.0, 8.5, 4.0], 
                [1.0, 9.0, 4.0]],
            [[1.0, 7.0, 5.0], [1.0, 7.5, 5.0], [1.0, 8.0, 5.0], [1.0, 8.5, 5.0], 
                [1.0, 9.0, 5.0]],
            [[1.0, 7.0, 6.0], [1.0, 7.5, 6.0], [1.0, 8.0, 6.0], [1.0, 8.5, 6.0], 
                [1.0, 9.0, 6.0]],
        ]

        self.clusters = {
            0: {'indices': [0, 1, 2], 'hidden': numpy.asarray(
                [[3., 3., 15.], [3., 6.75, 15.], [3., 10.5, 15.]]), 'N': 3,
                'skel': 1}, 
            1: {'indices': [3, 5, 7], 'hidden': numpy.asarray(
                [[3., 15., 3.], [3., 15., 6.75], [3., 15., 10.5]]), 'N': 3,
                'skel': 5}, 
            2: {'indices': [4, 6, 8], 'hidden': numpy.asarray(
                [[3., 15., 21.], [3., 15., 24.], [3., 15., 27.]]), 'N': 3,
                'skel': 6}, 
            3: {'indices': [9, 10, 11], 'hidden': numpy.asarray(
                [[3., 21., 15.], [3., 24., 15.], [3., 27., 15.]]), 'N': 3,
                'skel': 10}, 
        }

    
    def test_streamline(self) :
        medipy.io.save(self.tensors, "/home/lamy/tmp/tensors.nii.gz")
        fibers = medipy.diffusion.tractography.streamline(self.tensors,
            0.5, 0.2, numpy.pi/3, 2, "Euler", self.tensors.spacing, None)
        
        distinct_fibers = []
        for fiber in fibers :
            if fiber.tolist() not in distinct_fibers :
                distinct_fibers.append(fiber.tolist())
        
        self.assertEqual(self.fibers, distinct_fibers)
    
    def test_local_skeleton_clustering(self) :
        clusters = medipy.diffusion.tractography.local_skeleton_clustering(
            self.fibers, 2)
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
            self.fibers, 2)
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
