[![Build Status](https://travis-ci.org/d-meiser/mbo.png?branch=master)](https://travis-ci.org/d-meiser/mbo)
[![Coverage Status](https://coveralls.io/repos/d-meiser/mbo/badge.png?branch=master)] (https://coveralls.io/r/d-meiser/mbo)

mbo
===

The many body operators (mbo) library provides functions for the
numerical solution of quantum mechanical many body problems.  It is
designed to be efficient both in terms of space usage and execution
time.

The latest version of mbo can be obtained from the github repository at

(https://github.com/d-meiser/mbo)


Build instructions
==================

mbo uses [cmake](http://www.cmake.org/) for cross platform
configuration.

Linux
-----

On linux, mbo can be downloaded, configured, and installed with the
following steps:

```
git clone https://github.com/d-meiser/mbo
cd mbo
mkdir build
cd build
cmake ..
make -j4
sudo make install
```
Note that the last step typically requires super user privileges.

Windows
-------

On windows mbo can be built with the microsoft visual studio compilers.
The free version (visual studio express) is sufficient.
