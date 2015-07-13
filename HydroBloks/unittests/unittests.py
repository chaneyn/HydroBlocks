import unittest
import sys

class DynamicTopmodel(unittest.TestCase):

  def test_celerity_subsurface(self):
    from pyDTopmodel import dynamic_topmodel
    #Initialize the model
    nhru = 4
    nhru_outlet = 0
    dt = dynamic_topmodel.Dynamic_Topmodel(nhru,nhru_outlet)
    self.assertEqual('hello'.upper(),'HELLO')

sys.path.append('../')
suite = unittest.TestLoader().loadTestsFromTestCase(DynamicTopmodel)
unittest.TextTestRunner(verbosity=2).run(suite)
