#!/usr/bin/env python3
"""
Foo Template

Template for Python developments.

USAGE: foo.py --help # And the rest is handled by `argparse`
"""

import os
import argparse

__version__ = "0.1"

class Class:
    """
    Define Class here
    """
    def __init__(self):
        """
        Initialize Class
        """
        pass


def parse_args():
    """
    Parse command-line options
    """
    parser = argparse.ArgumentParser(description="Template for Python developments")
    parser.add_argument('arg', help="Mandatory argument [REQUIRED]")
    parser.add_argument('-o', '--optional', help="Optional argument")
    parser.add_argument('-f', '--flag', action="store_true", help="Optional flag")
    return(parser.parse_args())


def _test(arg, opt="."):
    """
    Define Function1, _e.g._ for testing stuff here
    - arguments:
    - returns:
    """
    print(f"Required command-line argument is: {arg}")
    return(os.stat(opt))


if __name__ == '__main__':
    args = parse_args()
    _test(args.arg)
    print("\nDone.\n")