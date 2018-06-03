# Stable.jl
A Julia package for distributions related to strictly stable processes. This package includes the following distributions (using Zolotarev's (C) form):
<ul>
  <li>StableUnilateral.jl - Unilateral stable distribution with parameter ![equation](http://latex.codecogs.com/gif.latex?%5Calpha%5Cin%280%2C2%5D),</li>
  <li>Stable.jl - Stable distribution with parameters 
![equation](http://latex.codecogs.com/gif.latex?%5Calpha%5Cin%280%2C2%5D%2C%20%5Cbeta%5Cin%5B-1%2C1%5D),</li>
  <li>StablePositive.jl - Stable distribution conditioned to be positive with parameters 
![equation](http://latex.codecogs.com/gif.latex?%5Calpha%5Cin%280%2C2%5D%2C%20%5Cbeta%5Cin%5B-1%2C1%5D),</li>
  <li>StableSupremum.jl - The supremum of a stable process on [0,1] with parameters 
![equation](http://latex.codecogs.com/gif.latex?%5Calpha%5Cin%280%2C2%5D%2C%20%5Cbeta%5Cin%5B-1%2C1%5D).</li>
</ul>

The following ethods are implemented for every distribution: rand, mean, params, minimum, maximum, insupport
For every distribution but StableSupremum.jl, the following are also implemented: mellin, cf, cdf, pdf, mgf
