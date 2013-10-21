#!/usr/bin/env python

"""\
:mod:'test_CLI_rvd27.py'--Test Command Line Interface
/***************************************************
Author:Fan Zhang
Time:10/20/2013
***************************************************/
"""
__license__="""Copyright (C) 2013 Fan Zhang <fzhang@wpi.edu>
Permission to use, copy, modify, and distribute this file for any
purpose is hereby granted
"""
import cli.tests
import unittest
import CLI_rvd27_v2
import rvd27
import numpy as np

from cli.tests import AppTest
from CLI_rvd27_v2 import RVD_CLI as CLI

class test_CLI_case(unittest.TestCase):
	def setUp(self):
		self.app = rvd27
		
	def test_gen_gibbs_case(self):
		seed = 10241978
		J_pre = 100
		M_pre = [100 for i in range(J_pre)]
		n_pre = 100
		N_pre = 3
		M0_pre = 10
		mu0_pre = 0.0025
		ns = np.zeros((N_pre, J_pre))
		
		for n_index in range(N_pre):
			for J_index in range(J_pre):
				ns[n_index, J_index] = n_pre
				
		r, theta, mu, M = self.app.generate_sample(M=M_pre, M0=M0_pre, mu0=mu0_pre, n=n_pre, N=N_pre, J=J_pre, seedint = seed)
		
		(phi, theta_s, mu_s) = self.app.mh_sample(r, ns)
		self.app.save_model('test1.hdf5', phi, mu=mu_s, theta=theta_s, r=r, n=ns, loc=None, refb=None)
		
		# test conditions
		self.assertEqual(mu, mu_s)
		self.assertEqual(theta, theta_s)
		self.assertEqual(phi[M0], M0_pre)
		self.assertEqual(phi[M], M_pre)
		self.assertEqual(phi[mu0], mu0_pre)

# create the test suit
def suite():
	suite = unittest.TestSuite()
	suite.addTest(test_CLI_case("test_gen_gibbs_case"))
	return suite
# test
if __name__ == "__main__":
	unittest.main(defaultTest = 'suite')
	
	
		
		