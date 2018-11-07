# This repository is a work in progress and MOST of the functionnalities probably do not work

# LatMRG

*A front-end software package using 
[LatticeTester](https://github.com/umontreal-simul/latcommon) 
to test and search for new multiple recursive random number generators*

## Dependencies

LatMRG currently depends on
* [LatticeTester](https://github.com/umontreal-simul/latcommon)
* [NTL](http://www.shoup.net/ntl/index.html)
* [Boost](https://www.boost.org/)

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
