#!/usr/bin/env python3

import os, sys
import unittest
import json

src_path = os.path.dirname(os.path.dirname((os.path.realpath(__file__))))
sys.path.append(src_path)
import emg_make_batch_from_nanuq

class Test_emg_make_batch_from_nanuq:
    def test_tests(self):
        self.assertEqual(emg_make_batch_from_nanuq.tests(), 1)

if __name__ == '__main__':
    unittest.main()
