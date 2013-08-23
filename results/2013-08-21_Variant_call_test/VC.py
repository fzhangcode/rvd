import os
import sys
# Insert the src/python/rvd27 directory at front of the path
rvddir = os.path.join('../../src/python/rvd27/')
sys.path.insert(0, rvddir)
import rvd27

os.system("python ../../src/python/rvd27/rvd27.py -v test_main -T 0.00005 -N 1000 -o 'test' Control.hdf5 Case1_0.hdf5")
