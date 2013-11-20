""" test_rvd30.py: test fixtures for RVD30 Model."""

import rvd30

def test_generate_sample():
    
    r, theta, mu, M = rvd30.generate_sample(N=3)
    
    
    phi_hat, mu_hat, theta_hat, M_hat = rvd30.estimate_mom(r, 100)
    
    print phi_hat
    
if __name__ == '__main__':
    test_generate_sample()