#!/usr/bin/env python3

import os, sys
import unittest
import json

src_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
sys.path.append(src_path)

from lib.nanuq  import Nanuq

# class TestConfigurator(unittest.TestCase):
#     def setUp(self):
#         self.configurator = Configurator()

#     def test_default_config_file(self):
#         expected = os.path.expanduser("~/.illumina/gapp_conf.json")
#         self.assertEqual(self.configurator.file, expected)
    
#     def test_init_attributes(self):
#         self.assertIsNotNone(self.configurator.instance)
#         self.assertIsNotNone(self.configurator.domain)
#         self.assertIsNotNone(self.configurator.token)
#         self.assertIsNotNone(self.configurator.testDefinitionId)
#         self.assertIsNotNone(self.configurator.phenotips_server)
#         self.assertIsNotNone(self.configurator.phenotips_auth)
#         self.assertIsNotNone(self.configurator.phenotips_secret)
    
#     def tearDown(self):
#         pass
    

if __name__ == '__main__':
    unittest.main()
    print('\nDone\n')