""" test_rvd29.py: test fixtures for RVD29 Model."""

import rvd29

def test_generate_sample():
    
	r, theta, mu, M = rvd29.generate_sample(N=3)
	
	phi_hat, mu_hat, theta_hat, M_hat = rvd29.estimate_mom(r, 100)
	
	print phi_hat
	
if __name__ == '__main__':
    test_generate_sample()