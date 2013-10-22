#!/usr/bin/env python

"""\
:mod:'CLI_rvd27.py'--Command Line Interface
/***************************************************
Author:Fan Zhang
Time:10/8/2013
***************************************************/
"""
__license__="""Copyright (C) 2013 Fan Zhang <fzhang@wpi.edu>
Permission to use, copy, modify, and distribute this file for any
purpose is hereby granted
"""

import cli.app
import cli.log
import os
import sys
import numpy as np

from cli.app import Application
from cli.app import CommandLineMixin
from cli.util import StringIO
from cli.log import LoggingMixin


import rvd27
# Append the parent folder to the python path
sys.path.append(os.path.join(os.path.dirname(__file__), '../'))

class RVD_CLI(LoggingMixin, CommandLineMixin, Application):

	# Populate our options, -h/--help is already there.
	description = """
                    RVD is a hierarchical bayesian model for identifying
                    rare variants from short-read sequence data. """
					
	#Temp_Directroy = "pileup"
					
	#add paramters prog='rvd', description= description
	def __init__(self, main=None, **kwargs):
		LoggingMixin.__init__(self, **kwargs)
		CommandLineMixin.__init__(self, **kwargs)
		Application.__init__(self, main, **kwargs)
		
	
	def configure(self, **kargs):
			for name, value in kargs.items():
				isExist = hasattr(self, name)
				if isExist:
					setattr(self, name, value)
				else:
					print("ABORT EXECUTION: NO exists function ", name)
					sys.exit(1)
	def main(self):
		""" 
		manipulator, running in the run function of Application.
		"""
		
		args = self.argparser.parse_args()
			
		#self.log.basicConfig(level=log_level, format='%(levelname)s:%(module)s:%(message)s')
		
		# Tell where is the pileup temp directory 
		self.log.info('Info: The pileup files are: ./pileup/')
		
		# Create pileup directory if it doesn't exist
		if not os.path.isdir(rvd27.Temp_Directroy):
			os.makedirs(rvd27.Temp_Directroy)
			
		args.func(args)
		pass
	
	def func_gibbs(self, args):
		rvd27.gibbs(args)
		pass
		
	def func_test(self, args):
		rvd27.test(args.controlHDF5Name, args.caseHDF5Name, args.T, args.N, args.outputFile)
		pass
		
	def func_gen(self, args):
		"""Returns a sample with n reads, N replicates, and J locations. 
		The parameters of the model are in the structure phi.
		and the output file is stored in the hdf5file."""
		
		r, theta, mu, M= rvd27.generate_sample(args.M, args.M0, args.mu0, args.n, args.N, args.J, args.seedint)
		
		phi = {'mu0': args.mu0, 'M0': args.M0, 'M': M}
		
		ns = np.zeros((args.N, args.J))
		for n_index in range(args.N):
			for J_index in range(args.J):
				ns[n_index, J_index] = args.n
		
		rvd27.save_model(args.outputFile, phi, mu = mu, theta = theta, r=r, n=ns, loc = None, refb = None)
		
	def func_Proc(self, args):
		
		# Check that BAM file exists
		assert os.path.isfile(args.bamFileName), "BAM file does not exist: %s" % bamFileName
		
		# Check that FASTA reference file exists
		assert os.path.isfile(args.fastaFileName), "FASTA file does not exist: %s" % fastaFileName
		
		pileupFileName = rvd27.make_pileup(args.bamFileName, args.fastaFileName, args.region)
		dcFileName = rvd27.make_depth(pileupFileName)
		rvd27.gibbs(args.outputfile, dcFileName)
		
		"""test(controlHDF5Name, caseHDF5Name)"""
		pass
		
	def setup(self):
		Application.setup(self)
		CommandLineMixin.setup(self)
		LoggingMixin.setup(self)
		#create the top-level parser
		self.add_param('--version', action='version', version='%(prog)s 2.7')
		"""self.add_param('-v', '--verbose', dest='verbose', action='count', 
		                help="increase verbosity (specify multiple times for more)")"""
		
		subparsers = self.argparser.add_subparsers(help='sub-command help')
		
		#create subparser for Gibbs fitting
		argpgibbs = subparsers.add_parser('gibbs', help='fit the RVD model using Gibbs sampling')
		argpgibbs.add_argument('dcfile', nargs='+', help='depth chart file name')
		argpgibbs.add_argument('-o', '--outputFile', dest='outputFile', default='output.hdf5', help='output HDF5 file name')
		argpgibbs.add_argument('-p', '--pool', type=int, default=None, help='number of workers in multithread pool')
		argpgibbs.set_defaults(func=self.func_gibbs)
		
		#create subparser to compare two model files
		argptest = subparsers.add_parser('test', help='test if case error rate is greater than control by T')
		argptest.add_argument('controlHDF5Name', help='control model file (HDF5)')
		argptest.add_argument('caseHDF5Name', help='case model file (HDF5)')
		argptest.add_argument('-T', type=float, default=0.005, help='threshold for computing variant probability (default = 0.005)')
		argptest.add_argument('-N', type=int, default=1000, help='Monte Carlo Sample size (default=1000)')
		argptest.add_argument('-o', '--output', dest='outputFile', nargs='?', default='test.hdf5')
		argptest.set_defaults(func=self.func_test)
		
		#create subparser to sample the model
		argpgen = subparsers.add_parser('gen', help='Sample Data from the RVD model')
		argpgen.add_argument('-M', default=None, help='precision(default=100)')
		argpgen.add_argument('-M0', type=int, default=10, help='beta prior distribution(default=10)')
		argpgen.add_argument('-mu0', type=float, default=0.0025, help='beta prior distribution(default=.0025)')
		argpgen.add_argument('-n', type=int, default=100, help='reads(default=100)')
		argpgen.add_argument('-N', type=int, default=3, help='replicates(default=3)')
		argpgen.add_argument('-J', type=int, default=100, help='locations(default=100)')
		argpgen.add_argument('-seedint', type=int, default=None, help='(default=None)')
		argpgen.add_argument('-o', '--ouput', dest='outputFile', nargs='?', default='Output.hdf5')
		argpgen.set_defaults(func=self.func_gen)
		
		pass
		
	def pre_run(self):
		Application.pre_run(self)
		CommandLineMixin.pre_run(self)
		LoggingMixin.pre_run(self)
		
		#print("""

#######################################
# Welcome to RVD CLI Python Prototype #
#######################################
   		
		#""")
		#prompt = "\n" + self.PromptStr
		#sys.stdout.write(prompt)
		#sys.stdout.flush()
		
		#Parse the arguments(defaults to parsing sys.argv)
		#self.parse_args()
		pass
	
		
	
		
if __name__ == '__main__':
	RVD_CLI().run()



