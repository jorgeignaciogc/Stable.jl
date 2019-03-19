using Distributions, SpecialFunctions, FastGaussQuadrature

struct PositiveStable <: ContinuousUnivariateDistribution
  α::Float64
  β::Float64
  θ::Float64
  ρ::Float64
  PositiveStable(α::Real,β::Real) = ( β < -1 || β > 1 ||
    (β == -1 && α <= 1) || α <= 0 || α > 2 ) ?
    error("Parameters' requirements unmet:\n (α,β)∈(0,2]×[-1,1]-(0,1]×{-1}") :
    α == 2 ? new(2,0,0,.5) :
    new(α,β,β*(α <= 1 ? 1 :
    (α-2)/α), (1 + β*(α <= 1 ? 1 : (α-2)/α))/2)
end

import Distributions.params
function params(d::PositiveStable)
  return (d.α,d.β,d.θ,d.ρ)
end

import Distributions.minimum
function minimum(d::PositiveStable)
  return β == 1 && α == 1 ? 1 : 0
end

import Distributions.maximum
function maximum(d::PositiveStable)
  return β == 1 && α == 1 ? 1 : Inf
end

import Distributions.insupport
function insupport(d::PositiveStable,x::Real)
  return β == 1 && α == 1 ? x==1 : x >= 0
end

########################
# BASIC SIMULATION TOOLS
########################

import Distributions.rand
function rand(d::PositiveStable)
  ar = d.α*d.ρ
  if d.β == 1 && d.α == 1
    return 1.
  elseif d.β == 1 && d.α < 1
    u1 = rand(Uniform(0,pi))
    e1 = rand(Exponential())
    s1 = sin(u1)
    return (sin(d.α * u1) / s1) *
      (s1 * e1 / sin((1 -d.α) * u1)) ^ (1 - 1 /d.α)
  elseif d.β == -1 && d.α > 1
    u1 = rand(Uniform(0,pi))
    e1 = rand(Exponential())
    s1 = sin(u1)
    return ((sin(d.ρ * u1) / s1) *
      (s1 * e1 / sin((1 - d.ρ) * u1)) ^ (1 - 1 /d.ρ) ) ^ (-d.ρ)
  else
    u1 = rand(Uniform(0,pi))
    u2 = rand(Uniform(0,pi))
    e1 = rand(Exponential())
    e2 = rand(Exponential())
    s1 = sin(u1)
    s2 = sin(u2)
    return ( ((sin(ar * u1) / s1) *
      (s1 * e1 / sin((1 - ar) * u1)) ^ (1 - 1 /ar)) /
      ((sin(d.ρ * u2) / s2) *
      (s2 * e2 / sin.((1 - d.ρ) * u2)) ^ (1 - 1 /d.ρ)) ) ^ d.ρ
  end
end

function rand(d::PositiveStable,n::Int)
  ar = d.α*d.ρ
  if d.β == 1 && d.α == 1
    return ones(n)
  elseif d.β == 1 && d.α < 1
    u1 = rand(Uniform(0,pi),n)
    e1 = rand(Exponential(),n)
    s1 = sin.(u1)
    return (sin.(d.α .* u1) ./ s1) .*
      (s1 .* e1 ./ sin.((1-d.α) .* u1)) .^ (1-1/d.α)
  elseif d.β == -1 && d.α > 1
    u1 = rand(Uniform(0,pi),n)
    e1 = rand(Exponential(),n)
    s1 = sin.(u1)
    return ((sin.(d.ρ .* u1) ./ s1) .*
      (s1 .* e1 ./ sin.((1 - d.ρ) .* u1)) .^ (1 - 1/d.ρ) ) .^ (-d.ρ)
  else
    u1 = rand(Uniform(0,pi),n)
    u2 = rand(Uniform(0,pi),n)
    e1 = rand(Exponential(),n)
    e2 = rand(Exponential(),n)
    s1 = sin.(u1)
    s2 = sin.(u2)
    return ( ((sin.(ar .* u1) ./ s1) .*
      (s1 .* e1 ./ sin.((1 -ar) .* u1)) .^ (1 -1/ar)) ./
      ((sin.(d.ρ .* u2) ./ s2) .*
      (s2 .* e2 ./ sin.((1 -d.ρ) .* u2)) .^ (1 -1/d.ρ)) ) .^ d.ρ
  end
end

########################
# PDF, CDF
########################

function auxV2(x::AbstractArray,a::Real,b::Real)
  y = (pi/2).*x
  t = (pi/2)*a*b
  return ((sin.(a.*y .+ t).^a)./cos.(y)) .^ (1 /(1 -a)) .* cos.((a-1) .* y .+ t)
end

import Distributions.pdf
function pdf(d::PositiveStable,x::Real)
  if d.α <= 1 && β == 1
    return x == 1 ? 1. : 0.
  end
  if x < 0
    return 0.
  elseif x == 0
    return gamma(1 + 1 / d.α)*sin(pi*d.ρ)/(d.ρ*pi)
  end
  m = 10^3
  l = 15
  delta = 1
    #
  pir = pi*d.ρ
  pir1 = pi*(1 - d.ρ)
  #

  if d.α > 1
    a0 = d.α/(abs(d.α-1)*2)
    raa = 1 / (d.α-1)
    a1 = x == 0 ? 0 : abs(x)^raa
    a2 = -a1 .^ d.α
    s1 = d.ρ
    s2 = 1 - d.ρ
    a01 = a0*s1
    a02 = a0*s2
    if x<delta
      v = 1:Int(floor(l*d.α))
      w = (-1) .^ ((v .- 1) .% 2) .* gamma.(v ./ d.α .+ 1) .*
        sin.(pir .* v) ./ (pi .* gamma.(v .+ 1))
      return sum(w .* abs(x) .^ (v .- 1))
    else
      nodes, weights = gausslegendre(m)
      seq1 = auxV2(nodes .* s1 .+ s2,d.α,d.θ)
      return a01*a1*(sum(seq1 .* exp.(a2 .* seq1) .* weights))
    end
  else
    if x>delta
      v = 1:Int(floor(l*d.α))
      w = (-1) .^ ((v .- 1) .% 2) .* gamma.(v .* d.α .+ 1) .*
        sin.((pir*d.α) .* v) ./ gamma.(v .+ 1)
      return sum(w .* abs(x) .^ (-d.α .* v .- 1))/pi
    else
      nodes, weights = gausslegendre(m)
      nodes = nodes ./ 2 .+ .5
      weights = weights ./ 2
      nodes1 = nodes ./ (1 .- nodes)
      nodes2 = nodes ./ (1 .- nodes).^3
      pir = pi*d.ρ
      mat = abs(x) .* nodes1
      mat = pdf(UnilateralStable(d.α),mat)
      fC = nodes2 .* pdf(Cauchy(-cos(pir),sin(pir)),nodes1)
      return sum(fC .* mat .* weights)/d.ρ
    end
  end
end

function pdf(d::PositiveStable,X::AbstractArray{<:Real})
  if d.α == 1 && d.β == 1
    return convert(Array{Float64},X .== 1)
  end
  m = 10^3
  l = 15

  delta = 1
  x = vec(X)
  res = Float64[]
  #
  pir = pi*d.ρ
  pir1 = pi*(1 - d.ρ)
  v = 1:Int(floor(l*d.α))
  w = (-1) .^ ((v .- 1) .% 2) .* gamma.(v ./ d.α .+ 1) .* sin.(pir .* v) ./
    ((d.ρ*pi) .* gamma.(v .+ 1))
  w1 = (-1) .^ ((v .- 1) .% 2) .* gamma.(v ./ d.α .+ 1) .* sin.(pir1 .* v) ./
    ((d.ρ*pi) .* gamma.(v .+ 1))
  #
  aux = gamma(1 + 1 / d.α)*sin(pi*d.ρ)/(d.ρ*pi)
  if d.α > 1
    a0 = d.α/((d.α-1)*2)
    raa = 1 /(d.α-1)
    a1 = abs.(x) .^ raa
    a2 = -a1 .^ d.α
    nodes, weights = gausslegendre(m)
    s1 = d.ρ
    s2 = 1 - d.ρ
    seq1 = auxV2(nodes .* s1 .+ s2,d.α,d.θ)
    weights = (a0*s1) .* weights
    for i = 1:length(x)
      if x[i] > 0
        if x[i] < delta
          push!(res,sum(w .* abs(x[i]) .^ (v .- 1)))
        else
          push!(res,a1[i]*(sum(seq1 .* exp.(a2[i] .* seq1) .* weights)))
        end
      elseif x[i] < 0
        push!(res,0.)
      else
        push!(res,aux)
      end
    end
  else
    nodes, weights = gausslegendre(m)
    nodes = nodes ./2 .+ .5
    weights = weights ./ (d.ρ * 2)
    nodes1 = nodes ./ (1 .- nodes)
    nodes2 = nodes ./ (1 .- nodes).^3
    pir = pi*d.ρ
    fC = nodes2*pdf(Cauchy(-cos(pir),sin(pir)),nodes1)
    mat = abs.(x)*transpose(nodes1)
    mat = pdf(UnilateralStable(d.α),mat)
    w = (-1) .^ ((v .- 1) .% 2) .* gamma.(v .* d.α .+ 1) .*
      sin.((pir*d.α) .* v) ./ gamma.(v .+ 1)
    for i = 1:length(x)
      if x[i] > delta
        v = 1:Int(floor(l*d.α))
        push!(res,sum(w .* abs(x[i]) .^ (-d.α .* v .- 1))/pi)
      elseif x[i] > 0
        push!(res,sum(fC.*mat[i,:] .* weights))
      elseif x[i] < 0
        push!(res,0.)
      else
        push!(res,aux)
      end
    end
  end
  return reshape(res,size(X))
end

import Distributions.cdf
function cdf(d::PositiveStable,x::Real)
  if d.α == 1 && β == 1
    return x <= 1 ? 1. : 0.
  end
  if x <= 0
    return 0.
  end
  m = 10^3
  l = 15
  delta = 1.
  pir = pi*d.ρ
  pir1 = pi*(1-d.ρ)

  if d.α > 1
    raa = 1 / (d.α-1)
    a1 = x == 0 ? 0. : abs(x)^raa
    a2 = -a1 .^ d.α
    s1 = d.ρ
    s2 = 1 - d.ρ
    if x < delta
      v = 1:Int(floor(l*d.α))
      w = (-1) .^ ((v .- 1) .% 2) .* gamma.(v ./ d.α .+ 1) .*
        sin.(pir .* v) ./ (v .* gamma.(v .+ 1))
      return sum(w .* abs(x) .^ v)/(d.ρ*pi)
    else
      nodes, weights = gausslegendre(m)
      seq1 = auxV2(nodes .*s1 .+ s2,d.α,d.θ)
      return 1 - s1*(sum(exp.(a2 .* seq1) .* weights))/(2 *d.ρ)
    end
  else
    if x>delta
      v = 1:Int(floor(l*d.α))
      w = (-1) .^ ((v .- 1) .% 2) .* gamma.(v .* d.α .+ .1) .*
        sin.((pir*d.α) .* v) ./ (v .* gamma.(v .+ 1))
      return 1-sum(w .* abs(x) .^ (-d.α.*v))/(d.α*pi)
    else
      nodes, weights = gausslegendre(m)
      nodes = nodes ./ 2 .+ .5
      weights = weights ./ 2
      nodes1 = nodes ./ (1 .- nodes)
      nodes2 = 1 / (1 .- nodes) .^ 2
      pir = pi*d.ρ
      mat = abs(x)*nodes1
      mat = cdf(UnilateralStable(d.α),mat)
      fC = nodes2*pdf(Cauchy(-cos(pir),sin(pir)),nodes1)
      return sum(fC.*mat .* weights)/d.ρ
    end
  end
end

function cdf(d::PositiveStable,X::AbstractArray{<:Real})
  if d.α == 1 && β == 1
    return convert(Array{Float64},X .<= 1)
  end
  m = 10^3
  l = 15
  delta = 1.
  x = vec(X)
  res = Float64[]
  pir = pi*d.ρ
  pir1 = pi*(1 - d.ρ)
  v = 1:Int(floor(l*d.α))
  w = (-1) .^ ((v .- 1) .% 2) .* gamma.(v ./ d.α .+ 1) .*
    sin.(pir .* v) ./ ((pi*d.ρ) .* v .* gamma.(v .+ 1))
  if d.α > 1
    raa = d.α/(d.α-1)
    a2 = -abs.(x) .^ raa
    nodes, weights = gausslegendre(m)
    s1 = d.ρ
    s2 = 1 - d.ρ
    seq1 = auxV2(nodes .* s1 .+ s2,d.α,d.θ)
    s1 = s1/(2 * d.ρ)
    for i = 1:length(x)
      if x[i] > 0
        if x[i] < delta
          push!(res,sum(w .* (x[i] .^ v)))
        else
          push!(res,1 - s1*(sum(exp.(a2[i] .* seq1) .* weights)))
        end
      else
        push!(res,0.)
      end
    end
  else
    nodes, weights = gausslegendre(m)
    nodes = nodes ./ 2 .+ .5
    weights = weights ./ (2 *d.ρ)
    nodes1 = nodes./(1-nodes)
    nodes2 = 1 ./ (1 .- nodes).^2
    pir = pi*d.ρ
    fC = nodes2*pdf(Cauchy(-cos(pir),sin(pir)),nodes1)
    mat = abs(x)*transpose(nodes1)
    mat = cdf(UnilateralStable(d.α),mat)
    w = (-1) .^ ((v .- 1) .% 2) .* gamma.(v .* d.α .+ 1) .*
      sin.((pir*d.α) .* v) ./ (v .* gamma.(v .+ 1))
    for i = 1:length(x)
      if x[i] > delta
        v = 1:Int(floor(l*d.α))
        push!(res, 1 - sum(w .* abs(x[i]) .^ (-d.α .* v))/(d.α*pi))
      elseif x[i] > 0
        push!(res,aux+sum(fC .* mat[i,:] .* weights))
      else
        push!(res,0.)
      end
    end
  end
  return reshape(res,size(X))
end

import Distributions.mgf
function mgf(d::PositiveStable,x::Real)
  l = 15
  if x == 0
    return 1.
  end
  if d.α == 2
    return 2 * exp(x^2/4.)*cdf(Normal(0,sqrt(2)),x/2)
  end
  if d.β == -1 && x >= -1
    v = 0:Int(floor(l*d.α))
    w = 1 / gamma(v+1)
    return sum(w .* x .^ v)
  end

  nodes, weights = gausslegendre(m)
  nodes = nodes ./ 2 .+ .5
  weights = weights/(2 *d.ρ)
  nodes2 = 1 ./ (1 .- nodes) .^ 2
  nodes = nodes ./ (1 .- nodes)
  pir = pi*d.ρ
  fC = nodes2 .* pdf(Cauchy(-cos(pir),sin(pir)),nodes)
  if d.β == -1
    mat = exp(-(abs(x) .* nodes) .^ d.α)
    return sum(fC.*mat .* weights)
  else
    if x > 0
      return Inf
    else
      mat = exp(-(abs(x) .* nodes) .^ d.α)
      return sum(fC.*mat .* weights)
    end
  end
end

function mgf(d::PositiveStable,X::AbstractArray{<:Real})
  x = vec(X)
  l = 15
  if d.α==2
    return 2 .* exp.(x.^2 ./ 4.) .* cdf(Normal(0,sqrt(2)),x ./ 2)
  end
  nodes, weights = gausslegendre(m)
  nodes = nodes ./ 2 .+ .5
  weights = weights ./ (2 * d.ρ)
  nodes2 = 1 ./ (1 .- nodes) .^ 2
  nodes = nodes ./ (1 .- nodes)
  pir = pi*d.ρ
  fC = nodes2*pdf(Cauchy(-cos(pir),sin(pir)),nodes)
  res = Float64[]
    mat = exp.(-(abs(x)*transpose(nodes)) .^ d.α)
  if d.β == -1
    v = 0:Int(floor(18*d.α))
    w = 1 ./ gamma.(v .+ 1)
    if x[i] == 0
      push!(res,1.)
    elseif x[i] >= -1
      push!(res,sum(w .* x[i] .^ v))
    else
      push!(res,sum(fC.*mat[i,:] .* weights))
    end
  else
    for i = 1:length(x)
      if x[i] > 0
        push!(res,Inf)
      elseif x[i] < 0
        push!(res,sum(fC.*mat[i,:] .* weights))
      else
        push!(res,1.)
      end
    end
  end
  return reshape(res,size(X))
end

import Distributions.mean
function mean(d::PositiveStable)
  if d.α <= 1
    return Inf
  end
  return sin(pi*d.ρ)/(d.α*d.ρ*sin(pi/d.α)*gamma(1 + 1 / d.α))
end

import Distributions.var
function var(d::PositiveStable)
  if d.α < 2 && d.β != -1
    return Inf
  elseif d.α == 1
    if abs(d.β) != 1
      return Inf
    else
      return 0.
    end
  else
    return 2 / gamma(1 + 2 / d.α) - 1 / gamma(1 + 1 / d.α) ^ 2
  end
end

function mellin(d::PositiveStable,x::Complex)
  if (real(x) >= d.α && (d.α <= 1 || (d.α < 2 && d.β != -1))) || real(x) <= -1
    return Inf
  end
  if (d.α > 1 && d.β == -1) || d.α == 2
    return gamma(1 + x) / gamma(1 + x / d.α)
  end
  return (sin(pi * d.ρ * x) * gamma(1 + x)) /
    (d.α * d.ρ * sin(pi * x / d.α) * gamma(1 + x / d.α))
end

function mellin(d::PositiveStable,X::AbstractArray{<:Complex})
  if (d.α > 1 && d.β == -1) || d.α == 2
    return gamma.(1 .+ X) ./ gamma.(1 .+ X ./ d.α)
  end

  res = Complex{Float64}[]
  for x = X
    push!(res,
      (real(x) >= d.α && (d.α <= 1 || (d.α < 2 && d.β != -1))) || real(x) <= -1 ?
      Inf : (sin(pi * d.ρ * x) * gamma(1 + x)) / (d.α * d.ρ * sin(pi * x / d.α) * gamma(1 + x / d.α))
    )
  end
  return reshape(res,size(X))
end

function mellin(d::PositiveStable,x::Real)
  if (real(x) >= d.α && (d.α <= 1 || (d.α < 2 && d.β != -1))) || real(x) <= -1
    return Inf
  end
  if (d.α > 1 && d.β == -1) || d.α == 2
    return gamma(1 + x) / gamma(1 + x / d.α)
  end
  return (sin(pi * d.ρ * x) * gamma(1 + x)) /
  (d.α * d.ρ * sin(pi * x / d.α) * gamma(1 + x / d.α))
end

function mellin(d::PositiveStable,X::AbstractArray{<:Real})
  if (d.α > 1 && d.β == -1) || d.α == 2
    return gamma.(1 .+ X) ./ gamma.(1 .+ X ./ d.α)
  end

  res = Float64[]
  for x = X
    push!(res,
      (real(x) >= d.α && (d.α <= 1 || (d.α < 2 && d.β != -1))) || real(x) <= -1 ?
      Inf : (sin(pi * d.ρ * x) * gamma(1 + x)) / (d.α * d.ρ * sin(pi * x / d.α) * gamma(1 + x / d.α))
    )
  end
  return reshape(res,size(X))
end

export Stable, rand, minimum, maximum, insupport, pdf, cdf, mgf, mean, mellin, params
