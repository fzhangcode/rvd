""" test_rvd28.py: test fixtures for RVD28 Model."""

import rvd28

def test_generate_sample():
    
    r, theta, mu, M = rvd28.generate_sample(N=30)
    
    print r
    print mu
    print M
    
    phi_hat, mu_hat, theta_hat, M_hat = rvd28.estimate_mom(r, 100)
    
    print M_hat
    print phi_hat
    
if __name__ == '__main__':
    test_generate_sample()