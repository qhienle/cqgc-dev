#!/usr/bin/env python3
# USAGE: python -m unittest -v __file__

import os, sys
import unittest

LIB_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(LIB_DIR))
from gapp import Configurator
from gapp import Phenotips
from gapp import REDCap


class TestConfigurator(unittest.TestCase):
    def setUp(self):
        self.configs = Configurator()
    
    def test_load_configs_file(self):
        self.assertIsNotNone(self.configs.file, 'Could not load configuration file')

    def tearDown(self):
        pass


class TestPhenotips(unittest.TestCase):
    def setUp(self):
        self.pho  = Phenotips()
    
    def test_init_phenotips(self):
        self.assertIsNotNone(self.pho.server)
        self.assertIsNotNone(self.pho.auth)
        self.assertIsNotNone(self.pho.headers['Authorization'])
        self.assertIsNotNone(self.pho.headers['X-Gene42-Secret'])

    def test_get_hpo_old(self):
        pass

    def tearDown(self):
        pass


class TestREdCap(unittest.TestCase):
    def setUp(self):
        self.red  = REDCap()
    
    def test_init_redcap(self):
        self.assertIsNotNone(self.red.server)
        self.assertIsNotNone(self.red.token)

    def test_get_patient(self):
        foo = self.red.get_record_id('Q1K_HSJ_10050_P')
        self.assertEqual(foo, '50', 'REDCap record_id for patient "Q1K_HSJ_10050_P" should be "50"')
        self.assertIsInstance(foo, str, 'record_id should be an instance of `str`')

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
