from numpy import f2py
import numpy as np

with open("rhs-collection.f90") as sourcefile:
    sourcecode = sourcefile.read()

f2py.compile(sourcecode, modulename='rhs_collection', extension=".f90")

try:
    import rhs_collection 
    print("Successful compilation!")
    print(rhs_collection.__doc__)
    import os, glob
    compiled_module = glob.glob("*cpython*")[0]
    os.remove(compiled_module)
except Exception as ex: 
    print("Failed to compile.")
    print(f"{ex}")
