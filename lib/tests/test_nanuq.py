#!/usr/bin/env python3
# USAGE: python -m unittest -v __file__

import os, sys
import unittest

LIB_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(LIB_DIR))
from nanuq import Nanuq


class TestNanuq(unittest.TestCase):
    def setUp(self):
        self.nanuq  = Nanuq()
        self.run_id = "A00516_0428"
    
    def test_init_attributes(self):
        self.assertIsNotNone(self.nanuq.config_file)
        self.assertIsNotNone(self.nanuq.server)
        self.assertIsNotNone(self.nanuq.username)
        self.assertIsNotNone(self.nanuq.password)
        self.assertIsNotNone(self.nanuq.auth_data)

    def test_get_auth(self):
        auth = self.nanuq.get_auth()
        self.assertTrue(isinstance(auth, dict))

    def test_pass_credentials_as_args(self):
        self.nanuq = Nanuq(username="foo", password="bar")
        self.assertEqual(self.nanuq.username, 'foo')
        self.assertEqual(self.nanuq.password, 'bar')
        self.assertEqual(self.nanuq.auth_data['j_username'], 'foo')
        self.assertEqual(self.nanuq.auth_data['j_password'], 'bar')

    def test_get_api(self):
        url = self.nanuq.server
        self.assertRaises(Exception, self.nanuq.get_api(url))
        response = self.nanuq.get_api(url)
        self.assertEqual(response.status_code, 200)

    def test_get_clinical_sample(self):
        sample = self.nanuq.get_clinical_sample('22762')
        self.assertTrue(isinstance(sample, str))
    
    def test_get_sample(self):
        sample = self.nanuq.get_clinical_sample('32994')
        self.assertTrue(isinstance(sample, str))

    def test_get_sample_bad_request(self):
        self.assertRaises(Exception, self.nanuq.get_sample(''), 'Error 400: Bad request')

    def test_get_sample_not_found(self):
        self.assertRaises(Exception, self.nanuq.get_sample(21310), 'Error 404')
        self.assertRaises(Exception, self.nanuq.get_sample('00001'), 'Error 404')

    def test_get_sample_is_string(self):
        sample = self.nanuq.get_sample('21057')
        self.assertTrue(isinstance(sample, str))
        sample = self.nanuq.get_sample(21057)
        self.assertTrue(isinstance(sample, str))

    def test_check_downloaded_file(self):
        pass

    def test_check_run_name_long(self):
        fc_short = self.nanuq.check_run_name('200302_A00516_0106_BHNKHFDMXX')
        self.assertEqual(fc_short, 'A00516_0106', 'Convert RunID to short form')

    def test_check_run_name_is_ok(self):
        fc_short = self.nanuq.check_run_name('A00516_0106')
        self.assertEqual(fc_short, 'A00516_0106', 'Convert RunID to short form')

    # def test_check_run_name_is_bs(self):
    #     self.assertRaises(ValueError, self.nanuq.check_run_name('nimportequoi'), 'RunID is not the expected format')

    # def test_get_samplesheet(self):
    #     file = 'SampleSheet.csv'
    #     response = self.nanuq.get_samplesheet(self.run_id)
    #     self.assertEqual(response.status_code, 200, 'HTTP response code should be 200')
    #     self.assertFalse(os.stat(file).st_size == 0, 'File size should be greater than zero')
    #     os.remove(file)

    # def test_get_samplenames(self):
    #     file = 'SampleNames.txt'
    #     response = self.nanuq.get_samplenames(self.run_id)
    #     self.assertEqual(response.status_code, 200, 'HTTP response code should be 200')
    #     self.assertFalse(os.stat(file).st_size == 0, 'File size should be greater than zero')
    #     os.remove(file)

    # def test_get_samplepools(self):
    #     file = 'SamplePools.csv'
    #     response = self.nanuq.get_samplepools(self.run_id)
    #     self.assertEqual(response.status_code, 200, 'HTTP response code should be 200')
    #     self.assertFalse(os.stat(file).st_size == 0, 'File size should be greater than zero')
    #     os.remove(file)

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
