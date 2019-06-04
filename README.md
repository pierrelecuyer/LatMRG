# [LatMRG](https://savamarc.github.io/LatMRG)

*A software package to test and search for new linear congruential random number
generators*

## About this software

The documentation of this software is segmented in multiple locations that each
contain different information:
- This page that contains a list of features and installation instructions.
  Beware, since Github cannot render mathematical equations, you need to be
  a bit familiar with the subject to read it.
- A more comprehensive description of the software, both in terms of number of
  features listed and explanation of the concepts.
- A tutorial containing examples using both the executables and the library.
- The API documentation for the library.

*LatMRG* is a collection of executables and a library available freely to
research and test linear congruential random number generators. The goal of this
software is to provide reliable theoretical tests to the uniformity of the
large family of linear congruential random number generators. For short, in the
documentation of this library, linear congruential random number generators
might simply be called random number generators and even RNG. **Theoretical
tests are not an alternative to statistical testing** but rather a complement,
guarateeing a certain uniformity that is especially desirable in stochastic and
physics simulation.

Linear congruential random number generators possess a very particular **lattice
structure** that can be analysed to obtain a wide range of information about the
generator. *LatMRG* uses the [LatticeTester](https://github.com/umontreal-simul/latticetester)
library to perform operations on this lattice such as the spectral test but also
implements a few more functions to easily build those lattices and test the
period length of the random number generators.

Linear congruential random number generators are a large and widely used family
of random number generator that itself branches in many specific implementations
with different parameters. The library tries to extensively cover the
construction of such generators by implementing classes to represent various
types of RNG: simple and multiple recursive congruential generators, combined
recursive generators, add-with-carry and subtract-with-burrow,
multiply-with-carry and matrix congruential generators.

This software aims at both
- Providing extensive functionnality through executables
- Providing easy and ready to use and interface with representations of RNGs to
  expand the software and build 

## Getting it to work

*LatMRG* is a command line utility that is built and tested to work on Linux
only. This software should also work without much issues on macOS, but we do not
provide instructions for building on something else than Linux. Any user getting
the library to run on Windows or macOS would be welcome to provide us with his
or her process so that it can be added to this guide.

### Dependencies

LatMRG currently depends on
* [LatticeTester](https://github.com/umontreal-simul/latcommon): a utility library
upon which *LatMRG* builds. It currently is necessary to have this library in the
repository to compile *LatMRG*.
* [NTL](http://www.shoup.net/ntl/index.html): *LatMRG* heavily (and shamelesly)
uses the **Number Theory Library** developped by Victor Shoup. Make sure this is
installed with the NTL_THREADS=off option. This has to be in a standard path
as of now because *LatMRG* does not detect it otherwise.
* (*Factoring, free DLC*)[yafu](https://sourceforge.net/projects/yafu/): this factorization
utility unlocks some of the functionnality of *LatMRG*! To do so,
simply download the program and extract the `yafu` executable in `./data`. The
makefile will then include a preprocessor definition that will allow factoring.


### Configuring the build

LatMRG currently only has a very simple makefile that should work on most Linux
distributions. Assuming you have NTL installed in `/usr/local`, you can simply
call
```
git clone --recursive https://github.com/savamarc/LatMRG.git
cd LatMRG
make
```

This will pull and build the LatMRG library in `./LatMRG/lib`, and the executable
programs in `./LatMRG/bin`.

### Current maintainer(s)

[Marc-Antoine Savard](https://github.com/savamarc)
