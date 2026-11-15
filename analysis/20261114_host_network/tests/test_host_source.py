import importlib.util,pathlib,unittest
P=pathlib.Path(__file__).parents[1]/'scripts'/'host_source.py'
class Gates(unittest.TestCase):
 def test_implementation_exists(self): self.assertTrue(P.exists(), 'Source gates not implemented')
 def test_gates(self):
  self.assertTrue(P.exists(), 'Source gates not implemented')
  spec=importlib.util.spec_from_file_location('hs',P);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
  self.assertEqual(m.binomial('Urtica dioica'),'urtica dioica')
  self.assertIsNone(m.binomial('Urtica sp.'))
  self.assertIsNone(m.binomial('Urtica dioica subsp. dioica'))
  j={'matchType':'EXACT','rank':'SPECIES','kingdom':'Plantae','canonicalName':'Urtica dioica','status':'ACCEPTED','usageKey':123}
  self.assertEqual(m.plant_key('Urtica dioica',j),'123')
  self.assertIsNone(m.plant_key('Urtica dioica',dict(j,matchType='FUZZY')))
  self.assertIsNone(m.plant_key('Urtica dioica',dict(j,rank='GENUS')))
  self.assertIsNone(m.plant_key('Urtica dioica',dict(j,canonicalName='Urtica urens')))
  self.assertEqual(m.plant_key('Urtica dioica',dict(j,status='SYNONYM',acceptedUsageKey=456)),'456')
if __name__=='__main__':unittest.main()
