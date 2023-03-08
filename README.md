# Purpose

Wrap the righ-hand-side of the bunch1s models one by one with f2py.

# Usage

```shell
python3 python_test.py
```

Will compile all the subroutines in `rhs-collections.f90` and try to import the compiled module.
If it's successful it would print the compiled module documentation:

```shell
$ python3 python_test.py
Sucessful compilation!
This module 'rhs_collection' is auto-generated with f2py (version:1.24.1).
Functions:
    dy = gpmm2(y,m,bdef,p1,p2)
    dy = g1smm(y,m,bdef,p1,p2)
    dy = g1slw(y,m,par)
```

# Concept notes

We re-format the right-hand-sides to be compatible with the .f90 (Fortan-90) standart as it's more flexible with formatting
and automated formatters and linters (fprettify) exist that make the source code more readable and maintable.