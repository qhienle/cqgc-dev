#!/usr/bin/env python3

import os, sys
import unittest
import json

src_path = os.path.dirname(os.path.dirname((os.path.realpath(__file__))))
sys.path.append(src_path)
import emg_make_batch

class Test_emg_make_batch:
    def test_tests(self):
        self.assertEqual(emg_make_batch_from_nanuq.tests(), 1)

if __name__ == '__main__':
    unittest.main()

"""

# batch_manifest.csv
# (See Emedgene Help for specifications to the current file format). 
# For the example, "Files Names" truncated after the first two.

[Data],,,,,,,,,,,,,,,,,,,,,
Family Id,Case Type,Files Names,Sample Type,BioSample Name,Visualization Files,Storage Provider Id,Default Project,Execute_now,Relation,Gender,Phenotypes,Phenotypes Id,Date Of Birth,Boost Genes,Gene List Id,Kit Id,Selected Preset,Label Id,Clinical Notes,Due Date,Opt In
18-4493-T1,Whole Genome,/projects/3703703/biosamples/4407403/datasets/ds.4458d03cfc314c6a8d87ea9545a92226/sequenced files/297927630;/projects/3703703/biosamples/4407403/datasets/ds.4458d03cfc314c6a8d87ea9545a92226/sequenced files/297927631;,FASTQ,18-4493-T1,,10126,PRAGMatIQ_CHUS,False,proband,M,,HP:0001258;HP:0001264,2013-09-08,,,,Default,12,P0000468,,
18-4493-T1,Whole Genome,/projects/3703703/biosamples/4407404/datasets/ds.08fa2204208a409191b54c23829d5826/sequenced files/297927646;/projects/3703703/biosamples/4407404/datasets/ds.08fa2204208a409191b54c23829d5826/sequenced files/297927647;,FASTQ,24-07358-T1,,10126,PRAGMatIQ_CHUS,False,mother,F,Healthy,,1985-04-11,,,,Default,12,,,
18-4493-T1,Whole Genome,/projects/3703703/biosamples/4407405/datasets/ds.c55aa956bd9342b1a0eb1388dd91be5b/sequenced files/297928186;/projects/3703703/biosamples/4407405/datasets/ds.c55aa956bd9342b1a0eb1388dd91be5b/sequenced files/297928187;,FASTQ,24-07574-T1,,10126,PRAGMatIQ_CHUS,False,father,M,Healthy,,1975-07-04,,,,Default,12,,,
"""