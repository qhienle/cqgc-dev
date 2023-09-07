#!/usr/bin/env python3

"""
SampleSheet Class

Convert SampleSheet.csv files:

[Section_1]
key1,value1
key2,value2
key3,value3
...
[Section_2]
Col1,Col2,Col3,Col4
Val1,Val2,Val3,Val4
Val5,Val6,Val7,Val8
...

Into Python object:

self.sections = {
    'Header': [[key1,val1], [key2,val2], [key3,val3],...],
    'BCLConver_Settings': [[key1,val1], [key2, val2], [key3,val3],...],
    'BCLConvert_Data': [['col1','col2',...], [val,val,...], [val,val,...],
    'Cloud_Data': [['col1','col2',...], [val,val,...], [val,val,...],...],
    ...
}

Two versions of SampleSheets (indicated in the first section, [Header]).

- SampleSheet v1 ("IEMFileVersion,5" for some unknown reason) has 3 sections:
    - [Header] and [Settings]: comma-separated key-value pairs;
    - [Data]:  Lines (rows) of 11 comma-delimited fields, column names on the first line.
- SampleSheet v2 similar structure to v1, has different and additional sections:
    - [Header]
    - [Reads]
    - [BCLConvert_Settings]
    - [BCLConvert_Data]
    - [Cloud_Settings]
    - [Cloud_Data]
    - [CQGC_Data], a custom section for internal controls
"""

import os
import sys
import re
import copy
import time

__version__ = "0.1"

def now(format='full'):
    """
    Return Date-Time string for logging (timestamp).
    - `format`: can be 'full' (date-time) or 'time' only. Default='full'.
    """
    # import time
    if format == 'full':
        return(time.strftime('[%Y-%m-%d@%H:%M:%S]'))
    elif format == 'time':
        return(time.strftime('[@%H:%M:%S]'))
    else:
        print('WARNING: `now()` accepts "full" or "time". Default="full"')
        return(time.strftime('[%Y-%m-%d@%H:%M:%S]'))

class SampleSheet:
    """
    Initialize attributes with the function _init_() and populate from
     `file`, depending on the `version`. If file is not provided, an 
    object can be instanciated with empty attributes. Set `version` based on
    the value in the [Header] section. 

    self.sections = {
        'Header': [[key1,val1], [key2,val2], [key3,val3],...],
        'BCLConver_Settings': [[key1,val1], [key2, val2], [key3,val3],...],
        'BCLConvert_Data': [['col1','col2',...], [val,val,...], [val,val,...],
        'Cloud_Data': [['col1','col2',...], [val,val,...], [val,val,...],...],
        ...
    }
    """
    def __init__(self, file=None):
        self.file     = file
        self.version  = None
        self.sections = {'Header': []}
        
        # Parse the [Header] section, which we assume is the first one and
        # determine which file version to load (set self.version).
        #
        try:
            with open(file, 'r') as fh:
                lines = fh.readlines()
        except FileNotFoundError as error_fnf:
            print(f"{error_fnf}: {file} not found.")
            exit # or return(None)?

        for line in lines:
            section_re = re.search('\[(.+)\]', line)
            if section_re:
                section = section_re.group(1)
                if section in self.sections:
                    pass
                else:
                    self.sections[section] = []
            else:
                # We're on a csv-formatted line: split into list and append
                # non-empty lines (len() > 1) to `self.sections[{section}]`.
                #
                cols = line.rstrip().split(',')
                if section == 'Header':
                    if len(cols) == 2:
                        self.sections['Header'].append(cols)
                        # V1: self.header['IEMFileVersion'] == 5, key is not in v2
                        # V2: self.header['FileFormatVersion'] == 2, key is not in v1
                        if cols[0] == 'IEMFileVersion':
                            self.version = 1
                        elif cols[0] == 'FileFormatVersion':
                            self.version = 2
                else:
                    if len(cols) > 1:
                        self.sections[section].append(cols)

    def filter_samples(self, index_size=10):
        """
        Filter samples by "index_size" in a SampleSheet section [Data] (v1) or
        [BCLConvert_Data] (v2).
        - index_size: size of index (default=10)
        - Returns   : Object, modified list of samples, filtered by size of index.
        """

        # Determine SampleSheet version to access the sections [Data] (for v1)
        # or [BCLConvert_Data] (for v2), to populate list of samples.
        # Use var 'data' as a proxy to reference either of the two sections.
        # Var 'col' references the column with the index sequence in section.
        # SampleID at 2nd column in both versions. No need to reference, use '1'
        #
        if self.version == 1:
            data = 'Data'
            #col  = 7 # index at columns 7 and 9
        elif self.version == 2:
            data = 'BCLConvert_Data'
            #col  = 2 # index at columns 2 and 3

        samples = []
        samples.append(self.sections[data][0])  # Add the column headers
        samples_kept = [] # Track SampleIDs of given size of index
        for sample in self.sections[data]:
            if len(sample[2]) == index_size:
                samples.append(sample)
                samples_kept.append(sample[1])

        # Replace the original list of samples with the one we just filtered.
        #
        new_self = copy.deepcopy(self)
        new_self.sections[data] = samples

        # For v2 only: Filter-out 
        #
        if self.version == 2:
            for section in ['Cloud_Data', 'CQGC_Data']:
                samples_subset = []
                
                # Add the columns headers for this section, which is the first
                # item in this list
                #
                samples_subset.append(self.sections[section][0])
                for sample in self.sections[section]:
                    if sample[0] in samples_kept:
                        samples_subset.append(sample)
                new_self.sections[section] = samples_subset

        return(new_self)
        
    # TODO: Filter by library preparation kit (for Chromium samples)

    def remove_samples_with_one_index(self):
        """
        DEPRECATED: We now receive 10X libraries with 10-bp index
        Remove samples with only one index of 8 bases (which presumably 
        correspond to 10X Single-Cell)
        TODO: NOT IDEAL SOLUTION. Filter based on Info from LibraryPrepKit
        along with filtering samples based on index size
        """
        if self.version == 1:
            pass
        elif self.version == 2:
            data = 'BCLConvert_Data'
            col  = 2 # index at columns 2 and 3


    def add_base_mask(self, base_mask='Y101;I8N2;I8N2;Y101'):
        """
        Add base-mask to Settings (v1) or BCLConvert_Settings
        - SureSelect: OverrideCycles,Y101;I8N2;I8N2;Y101
        - 10XChromium: OverrideCycles,Y26;I8;Y101
        """
        if self.version == 1:
            settings = 'Settings'
            # Setting exists in v1, replace the base-mask value
            i = 0
            for setting in self.sections[settings]:
                if setting[0] == 'OverrideCycles':
                    self.sections[settings][i][1] = base_mask
                i += 1
        elif self.version == 2:
            settings = 'BCLConvert_Settings'
            self.sections[settings].append(['OverrideCycles', base_mask])
        return(self.sections[settings])

    def reverse_complement(seq):
        """
        Returns a reverse-complement of 'seq', as a string.
        """
        seq = seq.upper()
        seq = seq.replace('a', 'x')
        seq = seq.replace('t', 'a')
        seq = seq.replace('x', 't')
        seq = seq.replace('c', 'x')
        seq = seq.replace('g', 'c')
        seq = seq.replace('x', 'g')
        return(seq[::-1])

    def to_csv(self, file=None, version=2):
        """
        Write output to String content. Content can be written to screen,
        or to 'file', if specified.
        """
        content = ''
        if version == 2:
            order = ['Header', 'Reads', 'BCLConvert_Settings', 'BCLConvert_Data',
                     'Cloud_Settings', 'Cloud_Data', 'CQGC_Data']
            for section in order:
                content += f"[{section}]\n"
                for line in self.sections[section]:
                    content += ','.join(line) + "\n"
        elif version == 1:
            print(f"Sorry, cannot print to_csv() for SampleSheet v1, yet.")
            return()
        # TODO: Convert from v1 to v2 and vice-versa
        if file:
            with open(file, 'w') as fh:
                fh.write(content)
        else:
            print(content)

    def to_json():
        """Return self as a JSON string"""
        #TODO
        pass

    def to_dataframe():
        """Return self as a pandas dataframe"""
        # TODO
        pass

if __name__ == '__main__':

    if len(sys.argv) > 1:
        file = sys.argv[1]
    else:
        file = 'tests' + os.sep + 'Seq_S2XP_DS_20210602.csv'

    # Create a SampleSheet v1 instance
    sheet = SampleSheet(file)
    
    # Access attributes
    #
    print(f"File contains:\n{sheet.file}")
    print(f"Samplesheet version is {sheet.version}")
    if sheet.version == 1:
        print(f"Settings contains:\n{sheet.sections['Settings']}")
        print(f"Data contains:\n{sheet.sections['Data']}")
    elif sheet.version == 2: 
        print(f"Reads contains\n{sheet.sections['Reads']}")
        print(f"BCLConvert_Settings contains\n{sheet.sections['BCLConvert_Settings']}")
        print(f"BCLConvert_Data contains\n{sheet.sections['BCLConvert_Data']}")
        print(f"Cloud_Settings contains\n{sheet.sections['Cloud_Settings']}")
        print(f"Cloud_Data contains\n{sheet.sections['Cloud_Data']}")
        print(f"CQGC_Data contains\n{sheet.sections['CQGC_Data']}")
    
    # Add atrributes
    #
    result = sheet.add_base_mask()
    print(f"\nResult from add_base_mask():\n{result}")

    # Split the samplesheet
    #
    for index_size in [8, 10]:
        if sheet.version == 1:
            pass
        else:
            samples = sheet.filter_samples(index_size)
            print(f"\nSplit data contains:\n")
            for row in samples.sections['BCLConvert_Data']:
                print(row)
            if sheet.version == 2:
                print(f"\nSamples:\n{samples.sections}")

    # Output to csv
    #
    print(f"Testing output to_csv()\n")
    sheet.to_csv()
    print(sheet.sections)

    # Test reverse_complement
    #TODO

    print("\nDone.\n")