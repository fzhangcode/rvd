"""test_rvd27.py: Test fixtures for RVD2.7 model."""

import rvd27
import numpy as np
from nose.tools import assert_true, assert_false, assert_equal

def test_gamma_mle_empty():
    x = np.array([])
    (a,b) = rvd27.gamma_mle(x)
    assert_true(np.isnan(a) and np.isnan(b))