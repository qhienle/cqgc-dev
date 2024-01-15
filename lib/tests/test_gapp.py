#!/usr/bin/env python3
# USAGE: python -m unittest -v __file__

import os, sys
import unittest

LIB_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(LIB_DIR))
from gapp import Phenotips


class TestPhenotips(unittest.TestCase):
    def setUp(self):
        self.pho  = Phenotips()
    
    def test_init_attributes(self):
        self.assertIsNotNone(self.pho.server)
        self.assertIsNotNone(self.pho.auth)
        self.assertIsNotNone(self.pho.headers['Authorization'])
        self.assertIsNotNone(self.pho.headers['X-Gene42-Secret'])


    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
