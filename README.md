# LatMRG

C++ software tools to analyze the lattice structure of linear generators 
and to search for good ones 

## What this software is about

*LatMRG* is a C++ software library and tool to measure the multivariate uniformity of 
linear congruential and multiple recursive random number generators in terms 
of their lattice structure. It can also search for good parameters in terms of figures of merit
that measure the quality of this lattice structure, e.g., in terms of the spectral test
in up to several dimensions.  
One can analyze the lattice structure of points formed by successive
values in the generator's sequence, or formed by ``leapfrog'' values that can be far apart.
Generators with large moduli and multipliers (e.g. numbers of many
hundreds of bits), as well as combined generators, can also be analyzed.
Multiply-with-carry generators can be studied by analyzing their
corresponding linear congruential generators. 
Matrix multiple recursive generators can also be handled.
*LatMRG* relies on the *Lattice Tester* and *NTL* libraries.

## Documentation

More details on _LatMRG_, its underlying theory, its organization, and examples, can be found in the 
[**LatMRG User's Guide** (in .pdf)](https://www-labs.iro.umontreal.ca/~lecuyer/guides/latmrg-guide.pdf).
[](http://umontreal-simul.github.io/latmrg/)

The interface is specified in the 
[**API documentation**](http://pierrelecuyer.github.io/latmrg/namespaces.html).

## Compiling and Building

### Software Dependencies

Compiling *LatMRG* requires the following software to be installed:

* [LatticeTester](https://github.com/umontreal-simul/lattester)
  bundled and automatically compiled with LatMRG
* [GMP](https://gmplib.org/) compatible version with your NTL installation
* [NTL](http://www.shoup.net/ntl/index.html) 10.4.0 or later, installed with NTL_THREADS=off.
* [yafu](https://sourceforge.net/projects/yafu/) the `yafu` executable must be in the `./data`
  directory for the factorization function in `IntFactorization.h` to work properly
* [Git](http://git-scm.com/) *(optional for downloading the source code)*
* [Python](https://www.python.org/) *(Needed by waf to compile and build the library)*
* [Doxygen](http://www.stack.nl/~dimitri/doxygen/) *(optional for generating
  the API documentation)*

## Configuring the build

LatMRG currently only has a very simple makefile. If NTL is not installed in a
default prefix such as `/usr/local`, or if you use `clang` instead of `gcc` you
will need to modify it manually before building LatMRG. The following commands
will build the library and executables.
```
git clone --recursive https://github.com/umontreal-simul/LatMRG.git
cd LatMRG/latticetester
./waf configure
cd ..
make
```

This will pull and build the LatMRG library in `./LatMRG/lib`, and the executable
programs in `./LatMRG/bin`. There is currently no way to install LatMRG in
standard path to ease the usage of the library or invoke it via the command line.




### Building and Installing

Once everything is configured correctly, the following command will build the
*LatMRG* library and command-line tool:

    ./waf build

If the build process completed without errors, *LatMRG* can be installed to the
directory specified with the `--prefix` option during the configuration step,
with:

    ./waf install

## Authors

François Blouin, Erwan Bourceret, Anna Bragina, Ajmal Chaumun, 
Raymond Couture, Marco Jacques, David Munger, François Paradis, Marc-Antoine Savard, Richard Simard, 
Mamadou Thiongane, Josée Turgeon, and Christian Weiss
have contributed to various versions of this software since around 1986,
under the lead of Pierre L'Ecuyer.

## License

_Lattice Tester_ is free open source software, distributed under the Apache 2.0 License.

