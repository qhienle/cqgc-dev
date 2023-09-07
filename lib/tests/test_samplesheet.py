#!/usr/bin/env python3
# USAGE: python -m unittest __file__

import os, sys
import unittest

LIB_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(LIB_DIR))
from samplesheet import SampleSheet


class TestSampleSheet(unittest.TestCase):
    def setUp(self):
        workdir = os.path.dirname(os.path.abspath(__file__))
        os.chdir(workdir)
        self.sheet = SampleSheet('SampleSheet_v2.csv')
    

    def test_init_attributes(self):
        self.assertIsNotNone(self.sheet.file)
        self.assertIsNotNone(self.sheet.version)
        self.assertIsNotNone(self.sheet.sections)
    

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
