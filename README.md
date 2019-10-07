# SupStable.jl

A Julia package for exact simulation of the supremum of a stable process. It supports a few methods and an auxiliary distribution (see below for details). Specifically, this package includes the following distributions (using Zolotarev's (C) form of parametrization):
<ul>
<li>Stable - Stable random variable with parameters (α,β)∈(0,2]×[-1,1].</li>
<li>PositiveStable - Stable random variable conditioned to be positive with parameters (α,β)∈(0,2]×[-1,1]-(0,1]×{-1}.</li>
<li>SupremumStable - The supremum of a stable process on [0,1] with parameters (α,β)∈(0,2]×[-1,1].</li>
</ul>

## Table of Contents

1. [Schema](#schema) 

2. [Remarks and References](#references)

3. [Examples](#examples)

4. [Author and Contributor List](#authors)

<a name="schema"/>

## Schema

The distributions included support most of the standard functions as outlined in [Distributions.jl](https://github.com/JuliaStats/Distributions.jl).

### Stable - _Type_

```julia
Stable <: ContinuousUnivariateDistribution
```
This type has a single standard constructor `Stable(α::Real,β::Real)` with parameters (α,β)∈(0,2]×[-1,1]-(0,1]×{-1} and supports the methods `minimum`, `maximum`, `insupport`, `pdf`, `cdf`, `cf`, `mgf`, `mean`, `var`, `mellin`, `params` and `rand`.

#### Remarks

* Method `params(d::Stable)` returns the tuple (α,β,θ,ρ) following Zolotarev's (C) form (i.e.,
`ρ=1-cdf(d,0)` and `θ=2*ρ-1`).
* Method `mellin(d::PositiveStable,X::T)` returns the [Mellin transform](https://en.wikipedia.org/wiki/Mellin_transform), where `T` is either `F` or `AbstractArray{F}` and where `F` is either`Real` or `Complex`.
* Method `rand(d::Stable)` is based on [Chambers-Mellows-Stuck algorithm](https://en.wikipedia.org/wiki/Stable_distribution#Simulation_of_stable_variables).

### PositiveStable - _Type_

```julia
PositiveStable <: ContinuousUnivariateDistribution
```
This type has a single standard constructor `PositiveStable(α::Real,β::Real)` with parameters (α,β)∈(0,2]×[-1,1]-(0,1]×{-1} and supports the methods `minimum`, `maximum`, `insupport`, `pdf`, `cdf`, `cf`, `mgf`, `mean`, `var`, `mellin`, `params` and `rand`.

#### Remarks

* Method `params(d::PositiveStable)` returns the tuple (α,β,θ,ρ) following Zolotarev's (C) form.
* Method `mellin(d::PositiveStable,X::T)` returns the [Mellin transform](https://en.wikipedia.org/wiki/Mellin_transform), where `T` is either `F` or `AbstractArray{F}` and where `F` is either`Real` or `Complex`.

### SupremumStable - _Type_

```julia
SupremumStable <: ContinuousUnivariateDistribution
```
This type has a single standard constructor `SupremumStable(α::Real,β::Real)` with parameters (α,β)∈(0,2]×[-1,1] and supports the methods `minimum`, `maximum`, `insupport`, `mean`, `params`, `rand` and `sampler`.

#### Remarks

* Method `params(d::SupremumStable)` returns the tuple (α,β,θ,ρ) following Zolotarev's (C) form.
* If β=-1, constructor automatically defaults to `PositiveStable(α,β)` since they agree (see Theorem 3.1 in [(Michna, 2013)](https://doi.org/10.1214/ECP.v18-2236)).
* `sampler(d::SupremumStable)` returns a subtype [Sampler](https://juliastats.github.io/Distributions.jl/stable/extends.html) of sub type `PerfectSumpremumStable`. The optional arguments in `sampler(d::SupremumStable,args...)` are as in the constructor below of `PerfectSumpremumStable` below.
* `rand(d::SupremumStable)` calls `rand(sampler(d))[1]`.

### PerfectSupremumStable - _Type_

```
PerfectSupremumStable <: Sampleable{Univariate,Continuous}
```
This is an auxiliary sub type for generating exact samples of `SupremumStable` by means of _perfect simulation_, which has multiple hyperparameters (see the [references](#references) for more details and the conditions they must satisfy). The output of `rand` is a tuple `(x,s)` where `x` is the sample and `s` is the number of steps that the internal process ran for (beyond the user-defined warm up period `Δ`). This sub type supports `params` and the following constructors (for every omitted parameter, the constructor uses a suggested value):

* `PerfectSupremumStable(α::Real,β::Real,d::Real,δ::Real,γ::Real,κ::Real,Δ::Int,mAst::Int)`,
* `PerfectSupremumStable(α::Real,β::Real,d::Real,δ::Real,γ::Real,κ::Real,Δ::Int)`,
* `PerfectSupremumStable(α::Real,β::Real,d::Real,δ::Real,γ::Real,κ::Real)`,
* `PerfectSupremumStable(α::Real,β::Real,d::Real,δ::Real,γ::Real)`,
* `PerfectSupremumStable(α::Real,β::Real)`.

<a name="references"/>

## Remarks and References  

### StableSupremum 

This distribution's implementation relies on a recent paper by the authors of the package. See the article for details at: 
Jorge González Cázares and Aleksandar Mijatović and Gerónimo Uribe Bravo, *Exact Simulation of the Extrema of Stable Processes*, [arXiv:1806.01870v2](https://arxiv.org/abs/1806.01870v2) (2018). Consequently, some additional parameters are used with a Markov chain when sampling from it (see PerfectStableSupremum in stablesupremum.jl). In this reference, the variables `Δ` and `mAst` are denoted <a href="https://www.codecogs.com/eqnedit.php?latex=\Delta(0)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\Delta(0)" title="\Delta(0)" /></a> and <a href="https://www.codecogs.com/eqnedit.php?latex=m^\ast" target="_blank"><img src="https://latex.codecogs.com/gif.latex?m^\ast" title="m^\ast" /></a>, respectively.
Throughout the paper, the authors work on the parameters (α,ρ) where ρ is the positivity parameter and can be computed from (α,β) (see Appendix A in the reference).


<a name="examples"/>

## Examples

### Example 1  

In the case of infinite variation spectrally negative α-stable processes (i.e., when α>1 and β=-1), it is known that the StablePositive and StableSupremum distributions agree according to Theorem 3.1 in [(Michna, 2013)](https://doi.org/10.1214/ECP.v18-2236). We will check this empirically by comparing the empirical distribution function of multiple samples and the real distribution function and applying the [Kolmogorov-Smirnov test](https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test).

```julia
using StatsBase, Distributions, Random, SpecialFunctions

# Fix the random seed
Random.seed!(2019)

# Parameters
(α,β) = (1.5,-1.)

# Sample size
n = Int(1e4)
@time S = sort([rand(PerfectSupremumStable(α,β,2/3,1/3,.95*α,3,30,12))[1] for k = 1:n])
# Get the sample (we force the use of the perfect simulation)
Random.seed!(2019)
@time S = sort([rand(PerfectSupremumStable(α,β,2/3,1/3,.95*α,3,50,12))[1] for k = 1:n])
GC.gc()
# Theoretical CDF
F(x) = cdf(PositiveStable(α,β),x)

# We compute the distance between the empirical CDF and the theoretical CDF
h = 1/n
ks = max(maximum(abs.( F(S) - (h:h:1) )), maximum(abs.( F(S) - (0:h:(1-h)) )) )

# We compute the p-value
print("The p-value is ", 1-cdf(KSDist(n),ks))
# 0.6572

```

### Example 2  

This is the code that outputs (the first plot of) Figure 2 of the reference (the empirical CDF of the sample against the true CDF in the spectrally negative infinite variation case) if run repeatedly for all three parameter choices without resetting the seed.

```julia
using Distributions, StatsBase, SpecialFunctions, Random, FastGaussQuadrature, Gadfly, DataFrames

Random.seed!(2019)

(α,β) = (1.1,-1.)
n = Int(1e4)
df = DataFrame()

# Labels
lab = [i <= 100 ? "Theoretical" : "Practical" for i = 1:200]
df[:Label] = lab

# x-axis grid
x = convert(Array{Float64},0:(1.6/(100-1)):1.6)
df[:x] = vcat(x,x)
# y-axis coordinates of the theoretical CDF
y = cdf(PositiveStable(α,β),x)
# Get the sample and its empirical CDF
S = sort([rand(PerfectSupremumStable(α,β,2/3,1/3,α*.99,5,50,15))[1] for k = 1:n])
G = ecdf(S)
df[:y] = vcat(y, G(x))

# Plot comparing the empirical CDF and CDF
plot(df, x=:x, y=:y, color=:Label, Geom.line, Guide.xlabel("x"),
Guide.xticks(ticks=[0:.4:1.6;]), Guide.ylabel("F(x)"),
Guide.colorkey(title="CDF Type"),
Guide.Theme(panel_fill = nothing, major_label_color = "white",
key_title_color = "white", key_label_color = "white"),
Guide.title(string("ECDF-CDF comparison for α = ",α,", β = ",β)))

# Scaled error
df2 = DataFrame()
df2[:y] = cdf(PositiveStable(α,β),S)
df2[:z] = sqrt(n) .*(G(S) .- df2[:y])
plot(df2, x=:y, y=:z, Geom.line, Guide.xlabel("F(x)"),
Guide.xticks(ticks=[0:.2:1;]), Guide.ylabel("√n(Fn(x)-F(x))"),
Guide.Theme(panel_fill = nothing, major_label_color = "white"),
Guide.title(string("Scaled error (≈ Brownian bridge) for α = ",α,", β = ",β)))

print("\nKolmogorov-Smirnov statistic = ",maximum(abs.(df2[:z])),"\n")
print("\np-value = ",1-cdf(KSDist(n),maximum(abs.(df2[:z]))/sqrt(n)))
# 0.4604

# Confidence Bands
conf = 1.22217561146655
print("\nConfidence level = ",cdf(KSDist(n),conf/sqrt(n)))
# 0.9

```

<a name="authors"/>


## Author and Contributor List
Jorge I. González Cázares


Aleksandar Mijatović  
Gerónimo Uribe Bravo
