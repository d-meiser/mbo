[![Build Status](https://travis-ci.org/d-meiser/mbo.png?branch=master)](https://travis-ci.org/d-meiser/mbo)
[![Build status](https://ci.appveyor.com/api/projects/status/ti3765jb3ohy2q5t/branch/master?svg=true)](https://ci.appveyor.com/project/d-meiser/mbo/branch/master)
[![Coverage Status](https://coveralls.io/repos/d-meiser/mbo/badge.png?branch=master)] (https://coveralls.io/r/d-meiser/mbo)

mbo
===

The many body operators (mbo) library provides functions for the
numerical solution of quantum mechanical many body problems.  It is
designed to be efficient both in terms of space usage and execution
time.

The latest version of mbo can be obtained from the github repository at
[https://github.com/d-meiser/mbo](https://github.com/d-meiser/mbo).

The mbo website is here: [mbo website](https://d-meiser.github.io/mbo).

The documentation can be found here: [mbo documentation](http://d-meiser.github.io/mbo/docs/html/index.html).

mbo uses [cmake](http://www.cmake.org/) for cross platform
configuration.  On linux, mbo can be downloaded, configured, and
installed with the following steps:

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


Questions and Comments
======================

We'd love to hear how you're using mbo, what you find good, and what not
so good.  The preferred way for leaving feedback is via the github issue
tracker for the mbo repository
[https://github.com/d-meiser/mbo/issues](https://github.com/d-meiser/mbo/issues).

You can also contact me directly at dmeiser79@gmail.com.

