import unittest
import tofimaging.ReductionTools as rt

class TestReductionTools(unittest.TestCase):
    def test_Ang2Mev(self):
        self.assertAlmostEqual(rt.Ang2MeV(1),81.82)
        
    def test_MeV2Ang(self):
        self.assertAlmostEqual(rt.MeV2Ang(81.82),1)
        
    # def test_tof2l(self):