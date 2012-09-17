# -*- coding: utf8 -*-
import unittest
import numpy
import medipy.base
import medipy.gui.base

class TestPreferences(unittest.TestCase) :
    
    application_name = "MediPy-Tests"
    company_name = u"Université de Strasbourg"
    
    test_keys = [
        "basic/boolean", "basic/int", "basic/float", "basic/str", "basic/unicode", 
        "container/list", "container/tuple", "container/dict",
        "numpy/float64", "numpy/array"
    ]
    
    def __init__(self, *args, **kwargs) :
        preferences = medipy.gui.base.Preferences(self.application_name, self.company_name)
        for key in self.test_keys :
            if key in preferences :
                print ("Preference {0} already exists: "
                       "test will not be representative".format(key))
        
        unittest.TestCase.__init__(self, *args, **kwargs)
    
    def _test(self, key, value) :
        preferences = medipy.gui.base.Preferences(self.application_name, self.company_name)
        
        preferences.set(key, value)
        self.assertTrue(key in preferences)
        
        self.assertEqual(type(preferences.get(key)), type(value))
        if(isinstance(value, numpy.ndarray)) :
            numpy.testing.assert_equal(preferences.get(key), value)
        else :
            self.assertEqual(preferences.get(key), value)
        
        preferences.delete(key)
        self.assertFalse(key in preferences)
    
    def test_boolean(self) :
        self._test("basic/boolean", True)
    
    def test_int(self) :
        self._test("basic/int", 1234)
    
    def test_float(self) :
        self._test("basic/float", 12.34)
    
    def test_str(self) :
        self._test("basic/str", "abcd")
    
    def test_unicode(self) :
        self._test("basic/unicode", u"é€œ")
    
    def test_list(self) :
        self._test("container/list", [1,2,3,4])
    
    def test_tuple(self) :
        self._test("container/tuple", (1,2,3,4))
    
    def test_dict(self) :
        self._test("container/dict", {"first": 1, "second": 2})
    
    def test_float64(self) :
        self._test("numpy/float64", numpy.float64(12.34))
    
    def test_array(self) :
        self._test("numpy/array", numpy.asarray([[12.34, 34.56], [56.78], [78.91]]))

if __name__ == "__main__" :
    unittest.main()
