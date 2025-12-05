#!/usr/bin/env python3
# USAGE: python -m unittest -v __file__

import os, sys
import unittest

LIB_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(LIB_DIR))
from qlin import qlin


class TestQlin(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_load_url(self):
        q = qlin('https://qlin-me-hybrid.cqgc.hsj.rtss.qc.ca')
        self.assertIsNotNone(q.url, 'Server URL should not be NoneType')
        self.assertIsInstance(q.url, str, 'qlin.url should be of type str')
        
    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
