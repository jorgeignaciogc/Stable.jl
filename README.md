# StableDistributions.jl
A Julia package for distributions related to strictly stable processes. This package includes the following distributions (using Zolotarev's (C) form):
<ul>
  <li>StableUnilateral.jl - Unilateral stable distribution with parameter alpha in (0,2]</li>
  <li>Stable.jl - Stable distribution with parameters 
alpha in (0,2] and beta in [-1,1]</li>
  <li>StablePositive.jl - Stable distribution conditioned to be positive with parameters alpha in (0,2] and beta in [-1,1]</li>
  <li>StableSupremum.jl - The supremum of a stable process on [0,1] with parameters alpha in (0,2] and beta in [-1,1]</li>
</ul>

## Methods
<ul>
  <li>The following methods are implemented for every distribution: rand, mean, params, minimum, maximum, and insupport</li>
  <li>For every distribution but StableSupremum.jl, the following are also implemented: mellin, cf, cdf, pdf, and mgf</li>
</ul>

### Notes: StableSupremum.jl
This distributions' implementation relies on a recent paper by the authors of the package.
In particular, some additional parameters are used when sampling from it (see the article for details).
  
## Author and Contributor List
Jorge I. González Cázares

Aleksandar Mijatović  
Gerónimo Uribe Bravo
