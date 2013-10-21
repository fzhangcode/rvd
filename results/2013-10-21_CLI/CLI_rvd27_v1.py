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
from cli.app import Application
from cli.app import CommandLineMixin
from cli.util import StringIO
from cli.log import LoggingMixin


import rvd27
# Append the parent folder to the python path
sys.path.append(os.path.join(os.path.dirname(__file__), '../'))

class RVD_CLI(LoggingMixin, CommandLineMixin, Application):

	#add paramters prog='rvd', description= description
	def __init__(self, main=None, **kwargs):
		CommandLineMixin.__init__(self, **kwargs)
		LoggingMixin.__init__(self, **kwargs)
		Application.__init__(self, main, **kwargs)
		
		
		
		#self.PromptStr = "RVD_CLI>"
			
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
		args.func(args)
		pass
	
	def func_Gibbs(self, args):
		rvd27.gibbs(args)
		pass
		
	def func_Test(self, args):
		rvd27.test(args.controlHDF5Name, args.caseHDF5Name, args.T, args.N, args.outputFile)
		pass
		
	def func_Gen(self, args):
		"""Returns a sample with n reads, N replicates, and J locations. 
		The parameters of the model are in the structure phi.
		and the output file is stored in the hdf5file."""
		r, theta, mu, M = rvd27.generate_sample(args.phi, args.n, args.N, args.J, args.seedint)
		rvd27.save_model(args.outputfile, args.phi, mu, theta, r, args.n, args.J)
		pass
		
	def func_Proc(self, args):
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
		argpGibbs = subparsers.add_parser('Gibbs', help='fit the RVD model using Gibbs sampling')
		argpGibbs.add_argument('dcfile', nargs='+', help='depth chart file name')
		argpGibbs.add_argument('-o', '--outputFile', dest='outputFile', default='output.hdf5', help='output HDF5 file name')
		argpGibbs.add_argument('-p', '--pool', type=int, default=None, help='number of workers in multithread pool')
		argpGibbs.set_defaults(func=self.func_Gibbs)
		
		#create subparser to compare two model files
		argpTest = subparsers.add_parser('Test', help='test if case error rate is greater than control by T')
		argpTest.add_argument('controlHDF5Name', help='control model file (HDF5)')
		argpTest.add_argument('caseHDF5Name', help='case model file (HDF5)')
		argpTest.add_argument('-T', type=float, default=0.005, help='threshold for computing variant probability (default = 0.005)')
		argpTest.add_argument('-N', type=int, default=1000, help='Monte Carlo Sample size (default=1000)')
		argpTest.add_argument('-o', '--output', dest='outputFile', nargs='?', default='test.hdf5')
		argpTest.set_defaults(func=self.func_Test)
		
		#create subparser to sample the model
		argpGen = subparsers.add_parser('Gen', help='Sample Data from the RVD model')
		argpGen.add_argument('input', nargs='+')
		argpGen.add_argument('-o', '--ouput', dest='outputFile', nargs='?', default='Output.hdf5')
		argpGen.set_defaults(func=self.func_Gen)
		pass
		
	def pre_run(self):
		Application.pre_run(self)
		CommandLineMixin.pre_run(self)
		LoggingMixin.pre_run(self)
		#LoggingMixin.pre_run(self)
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



