# Purpose

Wrap the righ-hand-side of the bunch1s models one by one with f2py.

# Usage

Download the latest build library built for your platform (the code is automatically compiled when pushed to github):

[Windows](https://nightly.link/vasilvas99/bunch1s_f2py/workflows/compile-bunch1s/main/build_windows.zip?h=5e973cebe0192e1902ec70ffb6946e18b78a24b2)

[Linux](https://nightly.link/vasilvas99/bunch1s_f2py/workflows/compile-bunch1s/main/build_unix.zip?h=5e973cebe0192e1902ec70ffb6946e18b78a24b2)

Then include it your python script with:

```python
import bunch1s
```

Currently provides the following equation right hand sides:

```shell
python3 -c "import bunch1s; print(bunch1s.__doc__);"
...
This module 'bunch1s' is auto-generated with f2py (version:1.24.0).
Functions:
    dy = gpmm2(y,bdef,p1,p2,m=shape(y, 0))
    dy = g1smm(y,bdef,p1,p2,m=shape(y, 0))
    dy = g1slw(y,par,m=shape(y, 0))
    dy = gkrug(y,par,m=shape(y, 0))
    dy = g_pk2(y,par,m=shape(y, 0))
    dy = g_mm0(y,par,m=shape(y, 0))
    dy = g_mm1(y,par,m=shape(y, 0))
    dy = gise2(y,par,m=shape(y, 0))
```

# Concept notes

We re-format the right-hand-sides to be compatible with the .f90 (Fortan-90) standart as it's more flexible with formatting
and automated formatters and linters (fprettify) exist that make the source code more readable and maintainable.

# Compiling manually

To compile manually clone this repository and install: cmake, python, gfortran, g++, gcc, make (provided by mingw in windows).

Open a terminal in the root of the repository and run:

```shell
$ mkdir build
$ cmake -S . -B build/ -G "Unix Makefiles"
$ cd build
$ make
```

CMake will autoconfigure all libraries, dependencies and compilers and emit valid unix makefiles which upon invoking `make` inside the build folder will compile the project successfully.

# Initial results

Integrating the wrapped g1smm with LSODA from SciPy for parameters:

```
bdef = 2.0
p1 = 3.0
p2 = 2.0
```

Yields nice looking trajectories and grouping:

![g1smm](initial_results/g1smm.png)


## pybunch1s and integrate_trajectories

pybunch1s provides a high-level python interface for all the wrapped functions in order to enable type-checking and IDE suggestions.

integrate_trajectories is now a library that provides all the basic tools needed to integrate a model.

## surface_statistics

Provides lstat and hstat