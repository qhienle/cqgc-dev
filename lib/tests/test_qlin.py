#!/usr/bin/env python3
# USAGE: python -m unittest -v __file__

import os, sys
import unittest

LIB_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(LIB_DIR))
from qlin import qlin


class TestQlin(unittest.TestCase):
    def setUp(self):
        self.q = qlin('https://qlin-me-hybrid.cqgc.hsj.rtss.qc.ca')
    
    def test_load_url(self):
        # q = qlin('https://qlin-me-hybrid.cqgc.hsj.rtss.qc.ca')
        self.assertIsNotNone(self.q.url, 'Server URL should not be NoneType')
        self.assertIsInstance(self.q.url, str, 'qlin.url should be of type str')
        
    def test_attr_authenticatedHeaders(self):
        self.assertIsNotNone(self.q.authenticatedHeaders, 'Server URL should not be NoneType')
        self.assertIsInstance(self.q.authenticatedHeaders, dict, 'qlin.authenticatedHeaders should be of type dict')
        self.assertIsInstance(self.q.authenticatedHeaders['Authorization'], str, 'qlin.authenticatedHeaders["Authorization"] should be of type str')
 
    def test_search_analysis_by_aliquot(self):
        analyses = self.q.search_analysis(aliquot='40250')
        self.assertIsInstance(analyses, list, 'qlin.search_analysis() should return a list')
        with self.assertRaises(Exception):
            print('qlin.search_analysis() without params returns an error')
            self.q.search_analysis()

    def test_search_analysis_by_mrn(self):
        analyses = self.q.search_analysis(mrn='3554393')
        self.assertIsInstance(analyses, list, 'qlin.search_analysis() should return a list')

    def test_get_an_analysis(self):
        analysis_id = '822034'
        analyses    = self.q.get_an_analysis(analysis_id)
        self.assertIsInstance(analyses, dict, 'Expected a dict')
        
    def test_get_hpo(self):
        mrns = ['3554393', '3554229', '3552343', '3555439']
        self.assertIsInstance(self.q.get_hpo(mrns[0]), list, 'Expected a list')
        self.assertEqual(self.q.get_hpo(mrns[1])[0], 'HP:0001644', 'Expected HP:0001644')
        self.assertEqual(len(self.q.get_hpo(mrns[2])), 2, 'Expected 2')

    def test_get_hpo_mrn_notfound(self):
        self.assertIsNone(self.q.get_hpo('inexistant_mrn'), 'Expected NoneType')

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
