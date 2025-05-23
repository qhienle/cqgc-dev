#!/usr/bin/env python3
# USAGE: python -m unittest -v __file__

import os, sys
import unittest

LIB_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(LIB_DIR))
from gapp import Configurator
from gapp import Phenotips
from gapp import REDCap
from gapp import Emedgene


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
        self.sample = 'Q1K_HSJ_10050_P'
    
    def test_init_redcap(self):
        self.assertIsNotNone(self.red.server)
        self.assertIsNotNone(self.red.token)

    def test_get_record_id(self):
        foo = self.red.get_record_id(self.sample)
        self.assertEqual(foo, '50', 'REDCap record_id for patient "Q1K_HSJ_10050_P" should be "50"')
        self.assertIsInstance(foo, str, 'record_id should be an instance of `str`')

    def test_get_hpo_is_str(self):
        self.assertIsInstance(self.red.get_hpo(self.sample), str, 'record_id should be an instance of `str`')
    
    def test_get_hpo_equals_10(self):
        hpos_str = self.red.get_hpo('Q1K_HSJ_100123_P')
        num_hpos = len(hpos_str.split(';'))
        self.assertEqual(num_hpos, 10, 'Patient "Q1K_HSJ_10050_P" should have 10 HPO terms!')

    def tearDown(self):
        pass


class TestEmedgene(unittest.TestCase):
    def setUp(self):
        self.emg = Emedgene()

    def test_init_emedgene(self):
        self.assertIsNotNone(self.emg.username)
        self.assertIsNotNone(self.emg.password)
        self.assertIsNotNone(self.emg.prag_server)
        self.assertIsNotNone(self.emg.eval_server)

    def test_init_emedgene_case(self):
        self.assertIsNotNone(self.emg.case)
        self.assertIsInstance(self.emg.case, dict, "Emedgene Case should be a dict")
        self.assertIsInstance(self.emg.case['test_data'], dict, "Case test_data should be a dict")

    def test_authenticate_emedgene(self):
        auth = self.emg.authenticate()
        self.assertIsInstance(auth, str, '`auth key must be instance of str`')
        self.assertTrue(auth.startswith('Bearer '))

    def test_get_emg_id_GM221763(self):
        self.assertEqual(self.emg.get_emg_id('GM221763'), 'EMG398184424')

    def test_create_test_data(self):
        self.emg.create_test_data(presets='Genome-v1.1')
        self.assertIsInstance(self.emg.case['test_data'], dict, "Case test_data should be a dict")
        self.assertEqual(self.emg.case['test_data']['type'], 'Whole Genome')
        self.assertEqual(self.emg.case['test_data']['selected_preset_set'], 'Genome-v1.1')

    def test_create_test_data_proband(self):
        hpos = ['HP:0001508', 'HP:0001510', 'HP:0002190']
        proband = self.emg.create_patient(sample_name='25-04410-T1', biosample='35676', gender='Male', relation='Test Subject', phenotypes=hpos)
        self.assertIsInstance(proband, dict, "Case test_data should be a dict")
        self.assertIsInstance(proband['fastq_sample'], str, "Case test_data should be a dict")
        self.assertIsInstance(proband['phenotypes'], list)
        self.assertEqual(proband['id'], '25-04410-T1')
        self.assertEqual(proband['gender'], 'Male')
        self.assertEqual(proband['relationship'], 'Test Subject')
        self.assertEqual(len(proband['phenotypes']), 3)

    def test_submit_emg_case(self):
        #self.assertSomething()
        print(self.emg)

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
