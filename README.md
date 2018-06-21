# SupStable.jl
A Julia package for exact simulation of the supremum of a stable process. It supports a few methods and an auxiliary distribution (see below for details). Specifically, this package includes the following distributions (using Zolotarev's (C) form of parametrization):
<ul>
  <li>StableSupremum - The supremum of a stable process on [0,1] with parameters alpha in (0,2] and beta in [-1,1]</li>
  <li>StablePositive - Stable random variable conditioned to be positive with parameters alpha in (0,2] and beta in [-1,1]</li>
  <li>StableUnilateral - Stable random variable with parameter alpha in (0,1]</li>
  <li>Stable - Stable random variable with parameters alpha in (0,2] and beta in [-1,1]</li>
</ul>

## Methods
The following methods are implemented
<ul>
  <li>For both distributions: rand, mean, params, minimum, maximum, and insupport</li>
  <li>For StablePositive: mellin, cdf, pdf, cf, mgf</li>
</ul>

### Notes: StableSupremum
This distributions' implementation relies on a recent paper by the authors of the package. See the article for details at:  
Jorge González Cázares and Aleksandar Mijatović and Gerónimo Uribe Bravo, *Exact Simulation of the Extrema of Stable Processes*,  arXiv:1806.01870v1 (2018).  

In particular, some additional parameters are used with a Markov chain when sampling from it (see StableSupremumExact in stablesupremum.jl).

## Examples
    using StatsBase
    using Distributions
    using Gadfly
        
    # Parameters
    a = 1.5
    b = 1
    
    # Distributions
    dSup = StableSupremum(a,b)
    dPos = StablePositive(a,b)
        
    println(mellin(dPos, a/2))
    
    # Sample size
    n = 1000
        
    # Samples
    (sSup,sigma,counter) = rand(dSup,n)
    sPos = rand(dPos,n)
    b = rand(Bernoulli(dPos.rho),n)
    u = rand(Uniform(),n)
    
    # Check perpetuity
    sSup2 = sSup.*u.^(1/a) + sSup.*b.*(1-u).^(1/a)
    plot([ecdf(sSup), ecdf(sSup2)], 0, 5)
        

## Author and Contributor List
Jorge I. González Cázares

Aleksandar Mijatović  
Gerónimo Uribe Bravo
