"""test_rvd27.py: Test fixtures for RVD2.7 model."""

import rvd27
import numpy as np
from nose.tools import assert_true, assert_false, assert_equal

def test_pileup2dc(pileupFileName):
    print pileupFileName
    rvd27.make_depth(pileupFileName)

if __name__ == '__main__':
    pileupFileName="../../../../../freeze/baker_yeast/GSY1135/test/gen007_test_4.pileup"
    test_pileup2dc(pileupFileName)

    