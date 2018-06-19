using Distributions
using FastGaussQuadrature
importall Distributions

struct StableUnilateral <: ContinuousUnivariateDistribution
  a::Float64
  s::Float64
  StableUnilateral(a) = ( a<=0 || a>1 ) ? error("Parameter requirement unmet:\n α∈(0,2]") : new(a,1.0)
  StableUnilateral(a,s) = ( a<=0 || a>1 || abs(b)!=1) ? error("Parameter requirement unmet:\n α∈(0,2]") : new(a,b)
end

import Distributions.minimum
function minimum(d::StableUnilateral)
  return d.a==1 ? d.s : 0.0
end

import Distributions.params
function params(d::StableUnilateral)
  return d.a
end

import Distributions.maximum
function maximum(d::StableUnilateral)
  return d.a==1 ? d.s : Inf
end

import Distributions.insupport
function insupport(d::StableUnilateral,x::Real)
  return d.a==1 ? x==d.s : x*d.s>=0
end

########################
# BASIC SIMULATION TOOLS
########################

import Distributions.rand
function rand(d::StableUnilateral)
  if d.a==1
    return d.s
  end
  c = 1/d.a
  u = rand(Uniform(0,pi))
  return d.s*(sin(u*d.a)/sin(u)^c)*(sin(u*(1-d.a))/rand(Exponential()))^(c-1)
end

function rand(d::StableUnilateral,n::Int)
  if d.a==1
    return ones(n).*d.s
  end
  c = 1/d.a
  ra = 1-d.a
  c1 = c-1
  u = rand(Uniform(0,pi),n)
  return d.s.*(sin(u.*d.a)./sin(u).^c).*(sin(u.*ra)./rand(Exponential(),n)).^c1
end

import Distributions.pdf
function pdf(d::StableUnilateral,x::Real)
  m=10^3
  if d.a==1
    return x==d.s ? 1.0 : 0.0
  elseif d.s*x<=0
    return 0.0
  end
  c = 1/d.a-1
  C = d.a/(1-d.a)
  C1 = 1/(d.a-1)
  nodes, weights = gausslegendre(m)
  k = pi/2
  nodes = nodes*k + k
  nodes = sin(d.a.*nodes).^(d.a/(1-d.a)).*sin((1-d.a).*nodes)./sin(nodes).^(1/(1-d.a))
  return abs(x)^C1*dot(nodes.*exp(-nodes./(abs(x)^C)),weights)*C/2
end

function pdf(d::StableUnilateral,X::AbstractArray)
  m=10^3
  x=vec(X)
  if d.a==1
    return reshape([x[i]==d.s ? 1.0 : 0.0 for i=1:length(x)], size(X))
  end
  x = x.*d.s
  c = 1/d.a-1
  C = d.a/(1-d.a)
  C1 = 1/(d.a-1)
  nodes, weights = gausslegendre(m)
  k = pi/2
  weights = weights*(C/2)
  nodes = nodes*k + k
  nodes = sin(d.a.*nodes).^(d.a/(1-d.a)).*sin((1-d.a).*nodes)./sin(nodes).^(1/(1-d.a))
  return reshape([x[i]<=0 ? 0.0 : x[i]^C1*dot(nodes.*exp(-nodes./(x[i]^C)),weights) for i=1:length(x)],size(X))
end

import Distributions.cd
function cdf(d::StableUnilateral,x::Real)
  m=10^3
  if d.a==1
    return x>=d.s ? 1.0 : 0.0
  elseif d.s*x<=0
    return (1-d.s)/2
  end
  c = 1/d.a-1
  C = d.a/(1-d.a)
  nodes, weights = gausslegendre(m)
  k = pi/2
  nodes = nodes*k + k
  nodes = sin(d.a.*nodes).^(d.a/(1-d.a)).*sin((1-d.a).*nodes)./sin(nodes).^(1/(1-d.a))
  return (1-d.s)/2+d.s*dot(exp(-nodes./(abs(x)^C)),weights)/2
end

function cdf(d::StableUnilateral,X::AbstractArray)
  m=10^3
  x=vec(X)
  if d.a==1
    return reshape( [x[i]>=d.s ? 1.0 : 0.0 for i=1:length(x)],size(X))
  end
  x = d.s.*x
  c = 1/d.a-1
  C = d.a/(1-d.a)
  nodes, weights = gausslegendre(m)
  k = pi/2
  weights = d.s*weights/2
  nodes = nodes*k + k
  nodes = sin(d.a.*nodes).^C.*sin((1-d.a).*nodes)./sin(nodes).^(1/(1-d.a))
  c1 = (1-d.s)/2
  return reshape([x[i]<=0 ? c1 : c1+dot(exp(-nodes./(x[i]^C)),weights) for i=1:length(x)],size(X))
end

import Distributions.cf
function cf(d::StableUnilateral,x::Real)
  if x==0
    return 1.0
  end
  if d.a==1
    return exp(im*d.s*x)
  end
  return exp(-abs(x)^d.a*exp(-pi*d.a*sign(d.s*x)*im/2))
end

function cf(d::StableUnilateral,X::AbstractArray)
  x=vec(X).*d.s
  if x==0
    return reshape(ones(length(x)),size(X))
  end
  if d.a==1
    return reshape(exp(im.*x),size(X))
  end
  return reshape(exp(-abs(x)^d.a*exp(-pi*d.a*sign(x)*im/2)),size(X))
end

import Distributions.mgf
function mgf(d::StableUnilateral,x::Real)
  if d.a==1
    return exp(d.s*x)
  end
  if x==0
    return 1.0
  elseif d.s*x<0
    return exp(-(-d.s*x)^d.a)
  else
    return Inf
  end
end

function mgf(d::StableUnilateral,X::AbstractArray)
  x = vec(X).*d.s
  if d.a==1
    return reshape(exp(x),size(X))
  end
  res = zeros(0)
  for i=1:length(x)
    if x[i]==0
      push!(res,1.0)
    elseif x[i]<0
      push!(res,exp(-x[i]^d.a))
    else
      push!(res,Inf)
    end
  end
  return reshape(res,size(X))
end

function mellin(d::StableUnilateral,x::Complex)
  if real(x)>=d.a
    return Inf
  end
  return gamma(1-x/d.a)/gamma(1-x)
end

function mellin(d::StableUnilateral,x::Real)
  if x>=d.a
    return Inf
  end
  return gamma(1-x/d.a)/gamma(1-x)
end

export Stable, rand, minimum, maximum, insupport, pdf, cdf, cf, mgf, mellin, params
