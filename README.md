# This repository is a work in progress and MOST of the functionnalities probably do not work

# LatMRG

*A front-end software package using 
[LatticeTester](https://github.com/umontreal-simul/latcommon) 
to test and search for new multiple recursive random number generators*

## Dependencies

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

## Compiling

### Configuring the Build

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
