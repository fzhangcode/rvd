"""test_rvd26.py: Test fixtures for RVD2.6 model."""

import rvd26
import numpy as np
from nose.tools import assert_true, assert_false, assert_equal

def test_gamma_mle_empty():
    x = np.array([])
    (a,b) = rvd26.gamma_mle(x)
    assert_true(np.isnan(a) and np.isnan(b))
    
def test_generate_sample():

    phi = {'mu0':0.25, 'M0':2E3, 'a':1000, 'b':1}
    r, theta, mu, M = rvd26.generate_sample(phi, n=100, N=3, J=100, seedint=10241978)
    print r
    assert_equal(r[0,0],26)
    
if __name__ == '__main__':
    test_generate_sample()