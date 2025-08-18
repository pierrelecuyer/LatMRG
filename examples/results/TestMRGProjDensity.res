Types: NTL::ZZ, double 

=============================================================
Results for a small MRG example with m=13, k=3, a=(7,0,4). 
===================================================
We build the lattice and look at projections over successive coordinates.

Figure of merit primal succ, with BB.
coordinates      sqlen       merit           minmerit 
{1,2,3,4}       1  0.44285        0.44285
{1,2,3,4,5}     7  0.770305        0.44285
{1,2,3,4,5,6}   11  0.712803        0.44285
{1,2,3,4,5,6,7} 15  0.664498        0.44285
{1,...,8}       25  0.711607        0.44285
FOM value: 0.44285

===================================================
Then we look at other primal projections, over pairs and triples.

Figure of merit primal non-succ, with BB.
coordinates      merit        sqlen    minmerit  
{1,2}           1  1        1
{1,3}           1  1        1
{1,4}           1  1        1
{1,5}           1  1        1
{1,2,3}         1  1        1
{1,2,4}         1  1        1
{1,2,5}         1  1        1
{1,3,4}         5  2.23607        1
{1,3,5}         1  1        1
{1,4,5}         1  1        1
FOM value: 1

===================================================
We now examine the m-dual lattice, first for successive coordinates.

Figure of merit for dual.
coordinates      sqlen       merit           minmerit 
{1,2,3,4}       29  0.66143        0.66143
{1,2,3,4,5}     24  0.853946        0.66143
{1,2,3,4,5,6}   8  0.607881        0.607881
{1,2,3,4,5,6,7} 8  0.700048        0.607881
{1,...,8}       8  0.764366        0.607881
FOM value: 0.607881

===================================================
Then we look at the m-duals of the other projections.

Figure of merit dual non-succ, with BB.
coordinates      sqlen       merit           minmerit 
{1,2}           169  1        1
{1,3}           169  1        1
{1,4}           169  1        1
{1,5}           169  1        1
{1,2,3}         169  1        1
{1,2,4}         169  1        1
{1,2,5}         169  1        1
{1,3,4}         29  0.414243        0.414243
{1,3,5}         169  1        0.414243
{1,4,5}         169  1        0.414243
FOM value: 0.414243

===================================================
A closer look at the projection over coordinates {1,3,4}: 
Full basis B before taking projection {1,3,4}: 
[[1 7 10 9 0 1 4 2]
[0 1 7 10 9 0 1 4]
[0 0 1 7 10 9 0 1]
[0 0 0 13 0 0 0 0]
[0 0 0 0 13 0 0 0]
[0 0 0 0 0 13 0 0]
[0 0 0 0 0 0 13 0]
[0 0 0 0 0 0 0 13]
]
Basis for projection {1,3,4}: 
[[1 -3 -4]
[0 1 -6]
[0 0 13]
]
Basis for m-dual projection {1,3,4}: 
[[13 0 0]
[39 13 0]
[22 6 1]
]

=============================================================
Results for a MRG example with m = 9223372036854773561, a = (1145902849652723, 0, -1184153554609676). 
===================================================
We build the lattice and look at projections over successive coordinates.

Figure of merit primal succ, with BB.
coordinates      sqlen       merit           minmerit 
{1,2,3,4}       1  1.52588e-05        1.52588e-05
{1,2,3,4,5}     1.12465e+08  0.000223483        1.52588e-05
{1,2,3,4,5,6}   1.02613e+15  0.00817336        1.52588e-05
{1,2,3,4,5,6,7} 1.09841e+21  0.358334        1.52588e-05
{1,...,8}       4.67515e+23  0.678151        1.52588e-05
FOM value: 1.52588e-05

===================================================
Then we look at other primal projections, over pairs and triples.

Figure of merit primal non-succ, with BB.
coordinates      merit        sqlen    minmerit  
{1,2}           1  1        1
{1,3}           1  1        1
{1,4}           1  1        1
{1,5}           1  1        1
{1,2,3}         1  1        1
{1,2,4}         1  1        1
{1,2,5}         1  1        1
{1,3,4}         1.02412e+08  10119.9        1
{1,3,5}         1  1        1
{1,4,5}         1  1        1
FOM value: 1

===================================================
We now examine the m-dual lattice, first for successive coordinates.

Figure of merit for dual.
coordinates      sqlen       merit           minmerit 
{1,2,3,4}       9.56894e+15  4.91481e-07        4.91481e-07
{1,2,3,4,5}     9.56894e+15  0.000332039        4.91481e-07
{1,2,3,4,5,6}   9.56894e+15  0.0249593        4.91481e-07
{1,2,3,4,5,6,7} 9.56894e+15  0.541513        4.91481e-07
{1,...,8}       8.92522e+13  0.51637        4.91481e-07
FOM value: 4.91481e-07

===================================================
Then we look at the m-duals of the other projections.

Figure of merit dual non-succ, with BB.
coordinates      sqlen       merit           minmerit 
{1,2}           8.50706e+37  1        1
{1,3}           8.50706e+37  1        1
{1,4}           8.50706e+37  1        1
{1,5}           8.50706e+37  1        1
{1,2,3}         8.50706e+37  1        1
{1,2,4}         8.50706e+37  1        1
{1,2,5}         8.50706e+37  1        1
{1,3,4}         9.56894e+15  1.06058e-11        1.06058e-11
{1,3,5}         8.50706e+37  1        1.06058e-11
{1,4,5}         8.50706e+37  1        1.06058e-11
FOM value: 1.06058e-11

===================================================
A closer look at the projection over coordinates {1,3,4}: 
Full basis B before taking projection {1,3,4}: 
[[1 1145902849652723 6821832294614679088 2092544303070727312 1127730827481453799 1975258172869346907 179500854511123495 5284610717697798191]
[0 1 1145902849652723 6821832294614679088 2092544303070727312 1127730827481453799 1975258172869346907 179500854511123495]
[0 0 1 1145902849652723 6821832294614679088 2092544303070727312 1127730827481453799 1975258172869346907]
[0 0 0 9223372036854773561 0 0 0 0]
[0 0 0 0 9223372036854773561 0 0 0]
[0 0 0 0 0 9223372036854773561 0 0]
[0 0 0 0 0 0 9223372036854773561 0]
[0 0 0 0 0 0 0 9223372036854773561]
]
Basis for projection {1,3,4}: 
[[1 -2401539742240094473 2092544303070727312]
[0 1 1145902849652723]
[0 0 9223372036854773561]
]
Basis for m-dual projection {1,3,4}: 
[[9223372036854773561 0 0]
[22150294503972708037623944697662628353 9223372036854773561 0]
[-2751931234187192216358337832427291 -1145902849652723 1]
]
