from numpy import f2py
import numpy as np

with open("rhs-collection.f90") as sourcefile:
    sourcecode = sourcefile.read()

f2py.compile(sourcecode, modulename='rhs_collection', extension=".f90")

try:
    import rhs_collection 
    print("Sucessful compilation!")
    print(rhs_collection.__doc__)
except Exception as ex:
    print("Failed to compile.")
