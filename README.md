### jointfit_nova -- Toy implementation of NOvA fitting tools for NOvA-T2K

----------------------------------------------------------------------

  *Package author*:               C. Backhouse <c.backhouse@ucl.ac.uk>
  
  *Build system*:                 J. Wolcott <jwolcott@fnal.gov>

----------------------------------------------------------------------

This document briefly describes the build process and usage of the demonstration NOvA fitting tools contained in this package.

## Prerequisites
* A compiler that understands C++14 (tested with gcc)
* ROOT v6.00.00 or above

## Build

The toy fitting tools contained in this package (a simplified version of NOvA's full CAFAna framework) should be built using CMake.

CMake builds are 'out-of-source' --- that is, they use a dedicated build directory that you can put anywhere.  In this example, we'll use a build directory that is a subdirectory of the source root (that is, the directory containing this file).

By default, the generated libraries and test executables are installed back into the source tree at `lib/` and `bin/`, respecively.  If you'd like to put them somewhere else, pass `-DCMAKE_INSTALL_PREFIX=/path/to/dir` (where `/path/to/dir` is where you'd like them to be installed) in the `cmake` command below.  

```shell script
$ mkdir build
$ cd build
$ cmake ..
$ make install
```

## Usage

### Testing installation
A test script is shipped with the package in order to test that it built and installed correctly.  To use it, you must set up your environment.  Assuming that `$INSTALL_DIR` is the installation directory (the root of the source tree if you didn't specify an installation directory to `cmake`, or the directory you specified if you did):
  
  * Set `$JOINTFIT_INC` to point to `$INSTALL_DIR/inc`
  * Append `$INSTALL_DIR/lib` to `$LD_LIBRARY_PATH`
  
Once the environment is prepared, change directory to the `CAFAna` directory inside the source tree and run the macro:

```shell script
$ cd CAFAna
$ root -b -q load_libs.C test.C+
```

(the `load_libs.C` macro loads the libraries required for the test macro to function)

### Integrating into other software

Probably you want to build an executable against the libraries generated in the build step.  The installation directory provides all the needed ingredients:
* Include files: `$INSTALL_DIR/inc`
* Libraries:     `$INSTALL_DIR/lib`
* Binaries (just tests): `$INSTALL_DIR/bin` 