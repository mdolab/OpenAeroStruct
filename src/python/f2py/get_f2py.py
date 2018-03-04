from __future__ import print_function
#------------- VERY IMPORTANT ------------
import os,sys
tmp = sys.stdout
sys.stdout = sys.stderr
import numpy.f2py
sys.stdout = tmp
print (os.path.dirname(os.path.abspath(numpy.f2py.__file__)))
