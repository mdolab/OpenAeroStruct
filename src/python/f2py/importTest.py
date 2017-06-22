
from __future__ import print_function
import sys

print('Testing if module OpenAeroStruct can be imported...')
try:
    import OAS_API
except:
    print('Error importing OAS_API.so')
    import traceback
    traceback.print_exc()
    sys.exit(1)
# end try

print('Module OpenAeroStruct was successfully imported.')
